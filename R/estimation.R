#' Internal estimation functions
#'
#' @name estimation-internal
#' @keywords internal
#' @importFrom stats as.formula coef fitted lm.fit na.omit residuals var
NULL


#' Parse the xtfifevd formula
#' @keywords internal
#' @noRd
.parse_formula <- function(formula, data, id, time, instruments = NULL,
                           na.action = na.omit) {
  
  # Extract formula components
  formula_str <- as.character(formula)
  if (length(formula_str) != 3) {
    stop("Formula must be of the form: y ~ x1 + x2 | z1 + z2")
  }
  
  # Parse y ~ x | z structure
  rhs <- formula_str[3]
  lhs <- formula_str[2]
  
  if (!grepl("\\|", rhs)) {
    stop("Formula must include '|' to separate time-varying (x) and ",
         "time-invariant (z) variables.\n",
         "Example: y ~ x1 + x2 | z1 + z2")
  }
  
  parts <- strsplit(rhs, "\\|")[[1]]
  x_part <- trimws(parts[1])
  z_part <- trimws(parts[2])
  
  # Build formulas for extraction
  y_formula <- as.formula(paste("~", lhs))
  x_formula <- as.formula(paste("~", x_part))
  z_formula <- as.formula(paste("~", z_part))
  
  # Get variable names
  y_name <- all.vars(y_formula)
  x_names <- all.vars(x_formula)
  z_names <- all.vars(z_formula)
  iv_names <- NULL
  
  if (!is.null(instruments)) {
    iv_names <- all.vars(instruments)
  }
  
  # Check all variables exist
  all_vars <- c(y_name, x_names, z_names, iv_names, id, time)
  missing_vars <- setdiff(all_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }
  
  # Subset and handle NAs
  vars_needed <- unique(c(y_name, x_names, z_names, iv_names, id, time))
  data_sub <- data[, vars_needed, drop = FALSE]
  data_sub <- na.action(data_sub)
  
  # Sort by id, time
  data_sub <- data_sub[order(data_sub[[id]], data_sub[[time]]), ]
  
  # Extract matrices
  y <- data_sub[[y_name]]
  X <- as.matrix(data_sub[, x_names, drop = FALSE])
  Z <- as.matrix(data_sub[, z_names, drop = FALSE])
  
  panel_id <- data_sub[[id]]
  time_id <- data_sub[[time]]
  
  IV <- NULL
  if (!is.null(iv_names)) {
    IV <- as.matrix(data_sub[, iv_names, drop = FALSE])
  }
  
  # Panel info
  panels <- unique(panel_id)
  N_g <- length(panels)
  N <- length(y)
  
  # Compute panel sizes
  panel_sizes <- table(panel_id)
  T_bar <- mean(panel_sizes)
  T_min <- min(panel_sizes)
  T_max <- max(panel_sizes)
  
  # Check balanced
  balanced <- (T_min == T_max)
  
  # Validate time-invariant variables
  for (z_name in z_names) {
    z_var <- data_sub[[z_name]]
    by_panel <- tapply(z_var, panel_id, function(x) length(unique(x)))
    if (any(by_panel > 1)) {
      warning("Variable '", z_name, "' is not strictly time-invariant. ",
              "Results may be unreliable.")
    }
  }
  
  list(
    y = y,
    X = X,
    Z = Z,
    IV = IV,
    panel_id = panel_id,
    time_id = time_id,
    y_name = y_name,
    x_names = x_names,
    z_names = z_names,
    iv_names = iv_names,
    id_name = id,
    time_name = time,
    panels = panels,
    N = N,
    N_g = N_g,
    T_bar = T_bar,
    T_min = T_min,
    T_max = T_max,
    balanced = balanced,
    data = data_sub
  )
}


#' Fixed Effects estimation (Stage 1)
#' @keywords internal
#' @noRd
.fe_stage1 <- function(y, X, panel_id) {
  
  N <- length(y)
  k_x <- ncol(X)
  panels <- unique(panel_id)
  N_g <- length(panels)
  
  # Compute panel means
  panel_means <- function(v) {
    tapply(v, panel_id, mean)[as.character(panel_id)]
  }
  
  y_bar <- panel_means(y)
  X_bar <- apply(X, 2, panel_means)
  if (!is.matrix(X_bar)) X_bar <- matrix(X_bar, ncol = k_x)
  
  # Within transformation
  y_within <- y - y_bar
  X_within <- X - X_bar
  
  # FE regression (demeaned)
  fit_fe <- lm.fit(X_within, y_within)
  beta <- coef(fit_fe)
  names(beta) <- colnames(X)
  
  # Residuals
  e_it <- residuals(fit_fe)  # idiosyncratic residuals
  u_it <- y - X %*% beta     # combined: alpha_i + e_it
  
  # Time-averaged residual (= alpha_i + mean(e_it))
  u_bar <- tapply(u_it, panel_id, mean)
  
  # Variance components
  sigma2_e <- sum(e_it^2) / (N - N_g - k_x)
  
  # Between variance from panel means
  overall_mean <- mean(y - X %*% beta)
  alpha_hat <- u_bar - overall_mean
  sigma2_u <- var(alpha_hat)
  
  # Robust variance for beta (clustered by panel)
  # Bread: (X'X)^{-1}
  XtX <- crossprod(X_within)
  XtX_inv <- solve(XtX)
  
  # Meat: sum of cluster-specific outer products
  # For clustered (panel-robust) standard errors
  meat <- matrix(0, k_x, k_x)
  for (i in panels) {
    idx <- which(panel_id == i)
    Xi <- X_within[idx, , drop = FALSE]
    ei <- e_it[idx]
    # Compute cluster score: X_i' e_i (k_x x 1 vector)
    score_i <- as.vector(crossprod(Xi, ei))
    meat <- meat + outer(score_i, score_i)
  }
  
  V_beta_robust <- XtX_inv %*% meat %*% XtX_inv
  
  # Standard (non-robust) variance
  V_beta <- sigma2_e * XtX_inv
  
  # Compute X_bar_i as a matrix (N_g x k_x)
  X_bar_i <- t(sapply(panels, function(i) {
    colMeans(X[panel_id == i, , drop = FALSE])
  }))
  if (k_x == 1) {
    X_bar_i <- matrix(X_bar_i, ncol = 1)
    colnames(X_bar_i) <- colnames(X)
  }
  
  list(
    beta = beta,
    e_it = e_it,
    u_it = as.vector(u_it),
    u_bar = u_bar,
    y_bar_i = tapply(y, panel_id, mean),
    X_bar_i = X_bar_i,
    sigma2_e = sigma2_e,
    sigma2_u = sigma2_u,
    V_beta = V_beta,
    V_beta_robust = V_beta_robust,
    df_residual = N - N_g - k_x
  )
}


#' FEF Stage 2 estimation
#' @keywords internal
#' @noRd
.fef_stage2 <- function(u_bar, Z, panel_id, panels) {
  
  N_g <- length(panels)
  k_z <- ncol(Z)
  
  # Get one observation per panel (first row of each panel)
  first_row <- !duplicated(panel_id)
  Z_panel <- Z[first_row, , drop = FALSE]
  
  # OLS of u_bar on Z with intercept
  Z_aug <- cbind(1, Z_panel)
  fit <- lm.fit(Z_aug, u_bar)
  
  coefs <- coef(fit)
  alpha <- coefs[1]
  gamma <- coefs[-1]
  names(gamma) <- colnames(Z)
  
  # Residuals: chat_i = ubar_i - alpha - Z_i' gamma
  chat <- residuals(fit)
  
  list(
    gamma = gamma,
    alpha = alpha,
    chat = chat,
    Z_panel = Z_panel
  )
}


#' Pesaran-Zhou variance for gamma (Eq. 17)
#' @keywords internal
#' @noRd
.pz_variance_fef <- function(Z_panel, X_bar_i, chat, V_beta_robust, N_g) {
  
  k_z <- ncol(Z_panel)
  k_x <- ncol(X_bar_i)
  
  # Ensure X_bar_i is a matrix
  if (!is.matrix(X_bar_i)) {
    X_bar_i <- matrix(X_bar_i, ncol = 1)
  }
  
  # Center
  Z_mean <- colMeans(Z_panel)
  X_mean <- colMeans(X_bar_i)
  
  Zc <- sweep(Z_panel, 2, Z_mean)
  Xc <- sweep(X_bar_i, 2, X_mean)
  
  # Q_zz = (1/N) * Z_c' Z_c  (Eq. 8)
  Qzz <- crossprod(Zc) / N_g
  
  # V_zz = (1/N) * sum(chat_i^2 * (z_i - zbar)(z_i - zbar)')  (Eq. 19)
  Vzz <- matrix(0, k_z, k_z)
  for (i in seq_len(N_g)) {
    zi <- Zc[i, ]
    Vzz <- Vzz + chat[i]^2 * outer(zi, zi)
  }
  Vzz <- Vzz / N_g
  
  # Q_zxbar = (1/N) * Z_c' X_c  (Eq. 9)
  Qzxbar <- crossprod(Zc, Xc) / N_g
  
  # V_gamma = (1/N) * Qzz^{-1} * [Vzz + Qzxbar * N*V(beta) * Qzxbar'] * Qzz^{-1}
  Qzz_inv <- solve(Qzz)
  mid <- Vzz + Qzxbar %*% (N_g * V_beta_robust) %*% t(Qzxbar)
  V_gamma <- (1 / N_g) * Qzz_inv %*% mid %*% Qzz_inv
  
  # Clean up names
  rownames(V_gamma) <- colnames(Z_panel)
  colnames(V_gamma) <- colnames(Z_panel)
  
  V_gamma
}


#' FEF-IV Stage 2 estimation with instruments
#' @keywords internal
#' @noRd
.fef_iv_stage2 <- function(u_bar, Z, IV, panel_id, panels) {
  
  N_g <- length(panels)
  k_z <- ncol(Z)
  k_iv <- ncol(IV)
  
  if (k_iv < k_z) {
    stop("Number of instruments (", k_iv, ") must be >= number of ",
         "endogenous z-variables (", k_z, ")")
  }
  
  # Get one observation per panel (first row of each panel)
  first_row <- !duplicated(panel_id)
  Z_panel <- Z[first_row, , drop = FALSE]
  IV_panel <- IV[first_row, , drop = FALSE]
  
  # 2SLS: First stage - regress Z on IV
  IV_aug <- cbind(1, IV_panel)
  Z_hat <- matrix(0, N_g, k_z)
  for (j in seq_len(k_z)) {
    fit_j <- lm.fit(IV_aug, Z_panel[, j])
    Z_hat[, j] <- fitted(fit_j)
  }
  
  # Second stage - regress u_bar on Z_hat
  Z_hat_aug <- cbind(1, Z_hat)
  fit_iv <- lm.fit(Z_hat_aug, u_bar)
  
  coefs <- coef(fit_iv)
  alpha <- coefs[1]
  gamma <- coefs[-1]
  names(gamma) <- colnames(Z)
  
  # IV residuals: upsilon_i = ubar_i - alpha - Z_i' gamma_iv
  upsilon <- as.vector(u_bar) - alpha - as.vector(Z_panel %*% gamma)
  
  list(
    gamma = gamma,
    alpha = alpha,
    upsilon = as.vector(upsilon),
    Z_panel = Z_panel,
    IV_panel = IV_panel
  )
}


#' Pesaran-Zhou variance for gamma_IV (Eq. 51)
#' @keywords internal
#' @noRd
.pz_variance_fef_iv <- function(Z_panel, IV_panel, X_bar_i, upsilon,
                                V_beta_robust, N_g) {
  
  k_z <- ncol(Z_panel)
  k_iv <- ncol(IV_panel)
  k_x <- ncol(X_bar_i)
  
  if (!is.matrix(X_bar_i)) {
    X_bar_i <- matrix(X_bar_i, ncol = 1)
  }
  
  # Center all
  Z_mean <- colMeans(Z_panel)
  R_mean <- colMeans(IV_panel)
  X_mean <- colMeans(X_bar_i)
  
  Zc <- sweep(Z_panel, 2, Z_mean)
  Rc <- sweep(IV_panel, 2, R_mean)
  Xc <- sweep(X_bar_i, 2, X_mean)
  
  # Q_zr = (1/N) * Z_c' R_c
  Qzr <- crossprod(Zc, Rc) / N_g
  
  # Q_rr = (1/N) * R_c' R_c
  Qrr <- crossprod(Rc) / N_g
  
  # Q_rxbar = (1/N) * R_c' X_c
  Qrxbar <- crossprod(Rc, Xc) / N_g
  
  # V_rr = (1/N) * sum(upsilon_i^2 * (r_i - rbar)(r_i - rbar)')
  Vrr <- matrix(0, k_iv, k_iv)
  for (i in seq_len(N_g)) {
    ri <- Rc[i, ]
    Vrr <- Vrr + upsilon[i]^2 * outer(ri, ri)
  }
  Vrr <- Vrr / N_g
  
  # H_zr = (Q_zr * Q_rr^{-1} * Q_zr')^{-1} * Q_zr * Q_rr^{-1}
  Qrr_inv <- solve(Qrr)
  H_mid <- Qzr %*% Qrr_inv %*% t(Qzr)
  Hzr <- solve(H_mid) %*% Qzr %*% Qrr_inv
  
  # V_gamma = (1/N) * H_zr * [V_rr + Q_rxbar * N*V(beta) * Q_rxbar'] * H_zr'
  mid <- Vrr + Qrxbar %*% (N_g * V_beta_robust) %*% t(Qrxbar)
  V_gamma <- (1 / N_g) * Hzr %*% mid %*% t(Hzr)
  
  rownames(V_gamma) <- colnames(Z_panel)
  colnames(V_gamma) <- colnames(Z_panel)
  
  V_gamma
}


#' FEVD estimator
#' @keywords internal
#' @noRd
.estimate_fevd <- function(parsed) {
  
  # Stage 1: Fixed effects
  stage1 <- .fe_stage1(parsed$y, parsed$X, parsed$panel_id)
  
  # Stage 2: Decompose unit effects
  stage2 <- .fef_stage2(stage1$u_bar, parsed$Z, parsed$panel_id, parsed$panels)
  
  # Stage 3: Pooled OLS (point estimates equal FEF, SEs wrong - we don't use them)
  # We compute Pesaran-Zhou corrected SEs instead
  
  # Variance for gamma using Pesaran-Zhou Eq. 17
  V_gamma <- .pz_variance_fef(
    stage2$Z_panel, stage1$X_bar_i, stage2$chat,
    stage1$V_beta_robust, parsed$N_g
  )
  
  # Build combined coefficient vector
  k_x <- length(stage1$beta)
  k_z <- length(stage2$gamma)
  
  coefficients <- c(stage1$beta, stage2$gamma, "_cons" = stage2$alpha)
  
  # Build combined variance matrix (block diagonal)
  k_total <- k_x + k_z + 1
  V <- matrix(0, k_total, k_total)
  
  # Beta block
  V[1:k_x, 1:k_x] <- stage1$V_beta
  
  # Gamma block
  V[(k_x + 1):(k_x + k_z), (k_x + 1):(k_x + k_z)] <- V_gamma
  
  # Alpha variance (from Stage 2 OLS)
  Z_aug <- cbind(1, stage2$Z_panel)
  se_alpha <- sqrt(sum(stage2$chat^2) / (parsed$N_g - k_z - 1) *
                     solve(crossprod(Z_aug))[1, 1])
  V[k_total, k_total] <- se_alpha^2
  
  rownames(V) <- colnames(V) <- names(coefficients)
  
  # Fitted values
  fitted <- parsed$X %*% stage1$beta +
    parsed$Z %*% stage2$gamma +
    stage2$alpha
  
  list(
    coefficients = coefficients,
    vcov = V,
    beta = stage1$beta,
    gamma = stage2$gamma,
    intercept = stage2$alpha,
    residuals = stage1$e_it,
    fitted.values = as.vector(fitted),
    sigma2_e = stage1$sigma2_e,
    sigma2_u = stage1$sigma2_u,
    N = parsed$N,
    N_g = parsed$N_g,
    T_bar = parsed$T_bar,
    k_x = k_x,
    k_z = k_z,
    x_names = parsed$x_names,
    z_names = parsed$z_names,
    y_name = parsed$y_name,
    V_gamma_pz = V_gamma
  )
}


#' FEF estimator
#' @keywords internal
#' @noRd
.estimate_fef <- function(parsed) {
  
  # Stage 1: Fixed effects
  stage1 <- .fe_stage1(parsed$y, parsed$X, parsed$panel_id)
  
  # Stage 2: FEF regression
  stage2 <- .fef_stage2(stage1$u_bar, parsed$Z, parsed$panel_id, parsed$panels)
  
  # Pesaran-Zhou variance
  V_gamma <- .pz_variance_fef(
    stage2$Z_panel, stage1$X_bar_i, stage2$chat,
    stage1$V_beta_robust, parsed$N_g
  )
  
  # Build results (same structure as FEVD)
  k_x <- length(stage1$beta)
  k_z <- length(stage2$gamma)
  
  coefficients <- c(stage1$beta, stage2$gamma, "_cons" = stage2$alpha)
  
  k_total <- k_x + k_z + 1
  V <- matrix(0, k_total, k_total)
  V[1:k_x, 1:k_x] <- stage1$V_beta
  V[(k_x + 1):(k_x + k_z), (k_x + 1):(k_x + k_z)] <- V_gamma
  
  Z_aug <- cbind(1, stage2$Z_panel)
  se_alpha <- sqrt(sum(stage2$chat^2) / (parsed$N_g - k_z - 1) *
                     solve(crossprod(Z_aug))[1, 1])
  V[k_total, k_total] <- se_alpha^2
  
  rownames(V) <- colnames(V) <- names(coefficients)
  
  fitted <- parsed$X %*% stage1$beta +
    parsed$Z %*% stage2$gamma +
    stage2$alpha
  
  list(
    coefficients = coefficients,
    vcov = V,
    beta = stage1$beta,
    gamma = stage2$gamma,
    intercept = stage2$alpha,
    residuals = stage1$e_it,
    fitted.values = as.vector(fitted),
    sigma2_e = stage1$sigma2_e,
    sigma2_u = stage1$sigma2_u,
    N = parsed$N,
    N_g = parsed$N_g,
    T_bar = parsed$T_bar,
    k_x = k_x,
    k_z = k_z,
    x_names = parsed$x_names,
    z_names = parsed$z_names,
    y_name = parsed$y_name,
    V_gamma_pz = V_gamma
  )
}


#' FEF-IV estimator
#' @keywords internal
#' @noRd
.estimate_fef_iv <- function(parsed) {
  
  if (is.null(parsed$IV)) {
    stop("FEF-IV requires instruments. Provide 'instruments' argument.")
  }
  
  # Stage 1: Fixed effects
  stage1 <- .fe_stage1(parsed$y, parsed$X, parsed$panel_id)
  
  # Stage 2: FEF-IV regression
  stage2 <- .fef_iv_stage2(stage1$u_bar, parsed$Z, parsed$IV, parsed$panel_id, parsed$panels)
  
  # Pesaran-Zhou variance for IV (Eq. 51)
  V_gamma <- .pz_variance_fef_iv(
    stage2$Z_panel, stage2$IV_panel, stage1$X_bar_i, stage2$upsilon,
    stage1$V_beta_robust, parsed$N_g
  )
  
  k_x <- length(stage1$beta)
  k_z <- length(stage2$gamma)
  k_iv <- ncol(parsed$IV)
  
  coefficients <- c(stage1$beta, stage2$gamma, "_cons" = stage2$alpha)
  
  k_total <- k_x + k_z + 1
  V <- matrix(0, k_total, k_total)
  V[1:k_x, 1:k_x] <- stage1$V_beta
  V[(k_x + 1):(k_x + k_z), (k_x + 1):(k_x + k_z)] <- V_gamma
  
  # IV variance for alpha
  IV_aug <- cbind(1, stage2$IV_panel)
  se_alpha <- sqrt(sum(stage2$upsilon^2) / (parsed$N_g - k_iv - 1) *
                     MASS::ginv(crossprod(IV_aug))[1, 1])
  V[k_total, k_total] <- se_alpha^2
  
  rownames(V) <- colnames(V) <- names(coefficients)
  
  fitted <- parsed$X %*% stage1$beta +
    parsed$Z %*% stage2$gamma +
    stage2$alpha
  
  list(
    coefficients = coefficients,
    vcov = V,
    beta = stage1$beta,
    gamma = stage2$gamma,
    intercept = stage2$alpha,
    residuals = stage1$e_it,
    fitted.values = as.vector(fitted),
    sigma2_e = stage1$sigma2_e,
    sigma2_u = stage1$sigma2_u,
    N = parsed$N,
    N_g = parsed$N_g,
    T_bar = parsed$T_bar,
    k_x = k_x,
    k_z = k_z,
    k_iv = k_iv,
    x_names = parsed$x_names,
    z_names = parsed$z_names,
    iv_names = parsed$iv_names,
    y_name = parsed$y_name,
    V_gamma_pz = V_gamma
  )
}
