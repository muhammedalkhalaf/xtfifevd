#' Between/Within Variance Ratio for Time-Invariant Variables
#'
#' @importFrom stats sd
#' @description
#' Computes the between-panel and within-panel standard deviations for
#' specified variables, along with their ratio. This diagnostic helps assess
#' whether FEVD/FEF methods may improve upon standard FE estimation.
#'
#' @param data A data frame containing the panel data.
#' @param variables A character vector of variable names to analyze.
#' @param id Character string naming the panel identifier variable.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{variable}{Variable name}
#'   \item{sd_between}{Between-panel standard deviation}
#'   \item{sd_within}{Within-panel standard deviation}
#'   \item{bw_ratio}{Ratio of between to within SD}
#' }
#'
#' @details
#' For truly time-invariant variables, the within-panel SD should be zero
#' (or near-zero due to numerical precision), giving an infinite ratio.
#'
#' According to Plümper and Troeger (2007), FEVD/FEF methods tend to improve
#' upon standard FE when:
#' \itemize{
#'   \item The between/within ratio exceeds approximately 1.7
#'   \item The correlation between z and the unobserved unit effect is not too high
#' }
#'
#' @examples
#' # Create example data
#' set.seed(42)
#' N <- 50
#' T <- 5
#' id <- rep(1:N, each = T)
#' z_invariant <- rep(rnorm(N), each = T)  # Truly time-invariant
#' z_slow <- rep(rnorm(N), each = T) + rnorm(N * T, sd = 0.1)  # Slowly varying
#' x_varying <- rnorm(N * T)  # Time-varying
#'
#' data <- data.frame(id = id, z_inv = z_invariant,
#'                    z_slow = z_slow, x = x_varying)
#'
#' bw_ratio(data, c("z_inv", "z_slow", "x"), id = "id")
#'
#' @references
#' Plumper, T., & Troeger, V. E. (2007). Efficient Estimation of Time-Invariant
#' and Rarely Changing Variables in Finite Sample Panel Analyses with Unit Fixed
#' Effects. \emph{Political Analysis}, 15(2), 124-139.
#' \doi{10.1093/pan/mpm002}
#'
#' @export
bw_ratio <- function(data, variables, id) {
  
  # Validate inputs
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  
  missing_vars <- setdiff(c(variables, id), names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }
  
  panel_id <- data[[id]]
  
  results <- data.frame(
    variable = character(),
    sd_between = numeric(),
    sd_within = numeric(),
    bw_ratio = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (v in variables) {
    x <- data[[v]]
    
    # Between SD: SD of panel means
    panel_means <- tapply(x, panel_id, mean, na.rm = TRUE)
    sd_between <- sd(panel_means, na.rm = TRUE)
    
    # Within SD: SD of deviations from panel means
    mean_expanded <- panel_means[as.character(panel_id)]
    deviations <- x - mean_expanded
    sd_within <- sd(deviations, na.rm = TRUE)
    
    # Ratio
    if (sd_within > 1e-10) {
      ratio <- sd_between / sd_within
    } else {
      ratio <- Inf
    }
    
    results <- rbind(results, data.frame(
      variable = v,
      sd_between = sd_between,
      sd_within = sd_within,
      bw_ratio = ratio,
      stringsAsFactors = FALSE
    ))
  }
  
  # Print nicely via message() so output can be suppressed
  lines <- c(
    "",
    "Between/Within Variance Ratios",
    paste(rep("-", 60), collapse = ""),
    sprintf("%-15s %12s %12s %12s",
            "Variable", "Between SD", "Within SD", "B/W Ratio"),
    paste(rep("-", 60), collapse = "")
  )
  for (i in seq_len(nrow(results))) {
    if (is.finite(results$bw_ratio[i])) {
      lines <- c(lines, sprintf("%-15s %12.4f %12.4f %12.2f",
                  results$variable[i],
                  results$sd_between[i],
                  results$sd_within[i],
                  results$bw_ratio[i]))
    } else {
      lines <- c(lines, sprintf("%-15s %12.4f %12.4f %12s",
                  results$variable[i],
                  results$sd_between[i],
                  results$sd_within[i],
                  "Inf"))
    }
  }
  lines <- c(
    lines,
    paste(rep("-", 60), collapse = ""),
    "Note: B/W > 1.7 suggests FEVD/FEF may improve on FE",
    "      (threshold depends on corr(z,u); see Plumper & Troeger 2007)",
    ""
  )
  message(paste(lines, collapse = "\n"))
  
  invisible(results)
}
