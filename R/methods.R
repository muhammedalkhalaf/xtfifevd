#' Methods for xtfifevd objects
#'
#' @name xtfifevd-methods
#' @keywords internal
#' @importFrom stats coef vcov pnorm qnorm printCoefmat
NULL


#' @export
print.xtfifevd <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  
  cat("\n")
  cat(x$method, "Estimation Results\n")
  cat(rep("-", 50), sep = "")
  cat("\n\n")
  
  cat("Dep. variable:", x$y_name, "\n")
  cat("Method:       ", x$method, "\n")
  cat("Observations: ", x$N, "\n")
  cat("Groups:       ", x$N_g, "\n")
  cat("T (average):  ", round(x$T_bar, 2), "\n")
  cat("\n")
  
  cat("Coefficients:\n")
  print(round(x$coefficients, digits))
  cat("\n")
  
  invisible(x)
}


#' @export
summary.xtfifevd <- function(object, ...) {
  
  coefs <- object$coefficients
  se <- sqrt(diag(object$vcov))
  z_val <- coefs / se
  p_val <- 2 * pnorm(-abs(z_val))
  
  coef_table <- cbind(
    Estimate = coefs,
    `Std. Error` = se,
    `z value` = z_val,
    `Pr(>|z|)` = p_val
  )
  
  # Add significance stars
  stars <- ifelse(p_val < 0.001, "***",
                  ifelse(p_val < 0.01, "**",
                         ifelse(p_val < 0.05, "*",
                                ifelse(p_val < 0.1, ".", ""))))
  
  result <- list(
    call = object$call,
    method = object$method,
    coefficients = coef_table,
    stars = stars,
    N = object$N,
    N_g = object$N_g,
    T_bar = object$T_bar,
    sigma2_e = object$sigma2_e,
    sigma2_u = object$sigma2_u,
    y_name = object$y_name,
    x_names = object$x_names,
    z_names = object$z_names,
    k_x = object$k_x,
    k_z = object$k_z
  )
  
  class(result) <- "summary.xtfifevd"
  result
}


#' @export
print.summary.xtfifevd <- function(x, digits = max(3L, getOption("digits") - 3L),
                                   signif.stars = getOption("show.signif.stars"),
                                   ...) {
  
  cat("\n")
  cat(rep("=", 70), sep = "")
  cat("\n")
  cat(x$method, "Estimation Results", "                      xtfifevd 1.0.0\n")
  cat(rep("=", 70), sep = "")
  cat("\n")
  
  cat("Dep. variable:  ", x$y_name, "\n")
  cat("Method:         ", x$method, "\n")
  cat("Variance:        Pesaran-Zhou (2016)\n")
  cat(sprintf("Observations:    %-10d Groups:     %d\n", x$N, x$N_g))
  cat(sprintf("T (average):     %-10.2f\n", x$T_bar))
  cat(rep("-", 70), sep = "")
  cat("\n\n")
  
  # Print coefficients
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
               P.values = TRUE, has.Pvalue = TRUE, ...)
  
  cat("\n")
  cat(rep("-", 70), sep = "")
  cat("\n")
  cat("Time-varying (FE):     ", paste(x$x_names, collapse = ", "), "\n")
  cat("Time-invariant:        ", paste(x$z_names, collapse = ", "), "\n")
  cat(sprintf("sigma_e: %.4f    sigma_u: %.4f\n",
              sqrt(x$sigma2_e), sqrt(x$sigma2_u)))
  cat(rep("=", 70), sep = "")
  cat("\n")
  
  invisible(x)
}


#' @export
coef.xtfifevd <- function(object, ...) {
  object$coefficients
}


#' @export
vcov.xtfifevd <- function(object, ...) {
  object$vcov
}


#' @export
confint.xtfifevd <- function(object, parm, level = 0.95, ...) {
  
  cf <- coef(object)
  se <- sqrt(diag(vcov(object)))
  
  if (missing(parm)) {
    parm <- names(cf)
  } else if (is.numeric(parm)) {
    parm <- names(cf)[parm]
  }
  
  a <- (1 - level) / 2
  fac <- qnorm(1 - a)
  
  ci <- cbind(cf[parm] - fac * se[parm],
              cf[parm] + fac * se[parm])
  
  pct <- paste0(format(100 * c(a, 1 - a), digits = 3), "%")
  colnames(ci) <- pct
  rownames(ci) <- parm
  
  ci
}
