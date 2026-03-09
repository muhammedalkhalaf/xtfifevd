#' Panel Fixed Effects Estimation for Time-Invariant Variables
#'
#' @description
#' Estimates panel models with time-invariant regressors using FEVD, FEF, or
#' FEF-IV methods. Standard fixed effects estimation cannot identify
#' coefficients on time-invariant variables; these methods decompose or filter
#' the unit effects to recover these coefficients.
#'
#' @param formula A formula of the form `y ~ x1 + x2 | z1 + z2` where variables
#'   before `|` are time-varying and variables after `|` are time-invariant.
#' @param data A data frame containing the variables.
#' @param id Character string naming the panel (individual) identifier variable.
#' @param time Character string naming the time identifier variable.
#' @param method Estimation method: `"fevd"` (default), `"fef"`, or `"fef_iv"`.
#' @param instruments For `method = "fef_iv"`, a one-sided formula specifying

#'   instrumental variables, e.g., `~ iv1 + iv2`.
#' @param na.action How to handle missing values. Default is `na.omit`.
#'
#' @return An object of class `"xtfifevd"` containing:
#' \describe{
#'   \item{coefficients}{Named vector of all coefficients (beta, gamma, intercept)}
#'   \item{vcov}{Variance-covariance matrix using Pesaran-Zhou estimator}
#'   \item{beta}{Coefficients on time-varying variables (from FE stage)}
#'   \item{gamma}{Coefficients on time-invariant variables}
#'   \item{intercept}{Overall intercept}
#'   \item{residuals}{Idiosyncratic residuals from FE stage}
#'   \item{fitted.values}{Fitted values}
#'   \item{sigma2_e}{Variance of idiosyncratic error}
#'   \item{sigma2_u}{Variance of unit effects}
#'   \item{N}{Total number of observations}
#'   \item{N_g}{Number of groups (panels)}
#'   \item{T_bar}{Average time periods per panel}
#'   \item{method}{Estimation method used}
#'   \item{call}{The matched call}
#' }
#'
#' @details
#' ## Model
#' The panel model is:
#' \deqn{y_{it} = \alpha_i + x_{it}'\beta + z_i'\gamma + \varepsilon_{it}}
#'
#' where \eqn{x_{it}} are time-varying regressors, \eqn{z_i} are time-invariant
#' regressors, and \eqn{\alpha_i} are individual fixed effects.
#'
#' ## Stage 1 (All Methods)
#' Fixed effects regression of \eqn{y_{it}} on \eqn{x_{it}} yields consistent
#' \eqn{\hat{\beta}} and combined residuals \eqn{\hat{u}_{it} = \hat{\alpha}_i + \hat{\varepsilon}_{it}}.
#'
#' ## Stage 2
#' Time-averaged residuals \eqn{\bar{u}_i} are regressed on \eqn{z_i}:
#' \itemize{
#'   \item **FEF**: OLS of \eqn{\bar{u}_i} on \eqn{z_i} with intercept
#'   \item **FEF-IV**: 2SLS using instruments \eqn{r_i}
#'   \item **FEVD**: Same as FEF, then Stage 3 pooled OLS (point estimates identical)
#' }
#'
#' ## Variance Estimation
#' Uses Pesaran-Zhou (2016) Equation 17/51 which properly accounts for
#' estimation uncertainty from Stage 1. The naive pooled OLS standard errors
#' from FEVD Stage 3 are **inconsistent** and can understate true SEs by
#' factors of 2-5x or more.
#'
#' @examples
#' # Simulate panel data
#' set.seed(123)
#' N <- 100  # panels
#' T <- 10   # time periods
#' n <- N * T
#'
#' # Generate data
#' id <- rep(1:N, each = T)
#' time <- rep(1:T, N)
#' alpha_i <- rep(rnorm(N), each = T)  # Fixed effects
#' z <- rep(rnorm(N), each = T)        # Time-invariant
#' x <- rnorm(n)                        # Time-varying
#' y <- 1 + 2 * x + 0.5 * z + alpha_i + rnorm(n, sd = 0.5)
#'
#' data <- data.frame(id = id, time = time, y = y, x = x, z = z)
#'
#' # Estimate with different methods
#' fit_fevd <- xtfifevd(y ~ x | z, data = data, id = "id", time = "time")
#' summary(fit_fevd)
#'
#' fit_fef <- xtfifevd(y ~ x | z, data = data, id = "id", time = "time",
#'                     method = "fef")
#' summary(fit_fef)
#'
#' @references
#' Plumper, T., & Troeger, V. E. (2007). Efficient Estimation of Time-Invariant
#' and Rarely Changing Variables in Finite Sample Panel Analyses with Unit Fixed
#' Effects. \emph{Political Analysis}, 15(2), 124-139.
#' \doi{10.1093/pan/mpm002}
#'
#' Pesaran, M. H., & Zhou, Q. (2018). Estimation of time-invariant effects in
#' static panel data models. \emph{Econometric Reviews}, 37(10), 1137-1171.
#' \doi{10.1080/07474938.2016.1222225}
#'
#' @seealso [fevd()], [fef()], [fef_iv()], [bw_ratio()]
#'
#' @export
xtfifevd <- function(formula, data, id, time,
                     method = c("fevd", "fef", "fef_iv"),
                     instruments = NULL,
                     na.action = na.omit) {
  
  call <- match.call()
  method <- match.arg(method)
  
  # Parse the formula: y ~ x1 + x2 | z1 + z2
  parsed <- .parse_formula(formula, data, id, time, instruments, na.action)
  
  # Dispatch to appropriate estimator
  result <- switch(method,
                   "fevd" = .estimate_fevd(parsed),
                   "fef" = .estimate_fef(parsed),
                   "fef_iv" = .estimate_fef_iv(parsed))
  
  result$call <- call
  result$formula <- formula
  result$method <- toupper(method)
  result$method <- gsub("_", "-", result$method)
  
  class(result) <- "xtfifevd"
  result
}


#' @describeIn xtfifevd FEVD estimation (3-stage, Plümper-Troeger)
#' @export
fevd <- function(formula, data, id, time, na.action = na.omit) {
  xtfifevd(formula, data, id, time, method = "fevd", na.action = na.action)
}


#' @describeIn xtfifevd FEF estimation (2-stage, Pesaran-Zhou)
#' @export
fef <- function(formula, data, id, time, na.action = na.omit) {
  xtfifevd(formula, data, id, time, method = "fef", na.action = na.action)
}


#' @describeIn xtfifevd FEF-IV estimation with instruments
#' @export
fef_iv <- function(formula, data, id, time, instruments,
                   na.action = na.omit) {
  xtfifevd(formula, data, id, time, method = "fef_iv",
           instruments = instruments, na.action = na.action)
}
