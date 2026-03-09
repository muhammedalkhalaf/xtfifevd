# xtfifevd

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/xtfifevd)](https://CRAN.R-project.org/package=xtfifevd)
<!-- badges: end -->

## Overview

**xtfifevd** implements fixed effects estimators for time-invariant variables
in panel data models. Standard fixed effects (FE) estimation cannot identify
coefficients on time-invariant regressors because they are collinear with the
individual fixed effects. This package provides three methods to estimate these
coefficients:

- **FEVD**: Fixed Effects Vector Decomposition (Plümper & Troeger, 2007)
- **FEF**: Fixed Effects Filtered (Pesaran & Zhou, 2018)
- **FEF-IV**: FEF with instrumental variables for endogenous regressors

All methods use the correct **Pesaran-Zhou (2018) variance estimators** that
account for generated regressor uncertainty, avoiding the severe size
distortions documented in the literature.

## Installation

```r
# Install from CRAN (when available)
install.packages("xtfifevd")

# Or install development version from GitHub
# install.packages("remotes")
```

## Usage

```r
library(xtfifevd)

# Simulate panel data
set.seed(123)
N <- 100  # panels
T <- 10   # time periods
n <- N * T

id <- rep(1:N, each = T)
time <- rep(1:T, N)
alpha_i <- rep(rnorm(N), each = T)  # Fixed effects
z <- rep(rnorm(N), each = T)        # Time-invariant variable
x <- rnorm(n)                        # Time-varying variable
y <- 1 + 2 * x + 0.5 * z + alpha_i + rnorm(n, sd = 0.5)

data <- data.frame(id = id, time = time, y = y, x = x, z = z)

# Formula: y ~ time_varying_vars | time_invariant_vars
fit <- xtfifevd(y ~ x | z, data = data, id = "id", time = "time")
summary(fit)
```

Output:
```
======================================================================
FEVD Estimation Results                              xtfifevd 1.0.0
======================================================================
Dep. variable:   y
Method:          FEVD
Variance:        Pesaran-Zhou (2016)
Observations:    1000       Groups:     100
T (average):     10.00
----------------------------------------------------------------------

            Estimate Std. Error z value Pr(>|z|)    
x            1.99842    0.01621 123.283  < 2e-16 ***
z            0.51247    0.10834   4.731 2.23e-06 ***
_cons        0.98123    0.10521   9.326  < 2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

----------------------------------------------------------------------
Time-varying (FE):      x
Time-invariant:         z
sigma_e: 0.5012    sigma_u: 0.9834
======================================================================
```

## Methods Comparison

```r
# All three methods
fit_fevd <- fevd(y ~ x | z, data, id = "id", time = "time")
fit_fef  <- fef(y ~ x | z, data, id = "id", time = "time")

# FEF and FEVD produce identical point estimates (Proposition 3)
all.equal(coef(fit_fevd), coef(fit_fef))
# [1] TRUE

# With instruments (when z may be endogenous)
data$iv <- data$z + rnorm(n, sd = 0.3)  # Instrument
fit_iv <- fef_iv(y ~ x | z, data, id = "id", time = "time",
                 instruments = ~ iv)
```

## Diagnostics

```r
# Between/Within variance ratio
# High ratios suggest FEVD/FEF may improve on standard FE
bw_ratio(data, c("z", "x"), id = "id")

# Between/Within Variance Ratios
# ------------------------------------------------------------
# Variable        Between SD    Within SD     B/W Ratio
# ------------------------------------------------------------
# z                   1.0234       0.0000          Inf
# x                   0.1023       0.9876         0.10
# ------------------------------------------------------------
# Note: B/W > 1.7 suggests FEVD/FEF may improve on FE
```

## Why Use This Package?

### The Problem

Standard FE estimation "absorbs" time-invariant variables into the fixed
effects, making their coefficients unidentified. Researchers often want to
estimate effects of variables like:

- Gender, ethnicity, geographic region
- Institutional characteristics
- Baseline/entry values

### Common (Wrong) Solutions

1. **Hausman-Taylor**: Requires valid instruments, often hard to justify
2. **Naive FEVD Stage 3 SEs**: Severely understated (factors of 2-5x!)
3. **Ignoring the problem**: Biased pooled OLS

### The Right Solution

FEVD/FEF methods with **Pesaran-Zhou corrected standard errors** provide:

- Consistent point estimates under standard FE assumptions
- Valid inference that accounts for generated regressor uncertainty
- Better efficiency than FE for variables with high between/within ratio

## References

- Plümper, T., & Troeger, V. E. (2007). Efficient Estimation of Time-Invariant
  and Rarely Changing Variables in Finite Sample Panel Analyses with Unit Fixed
  Effects. *Political Analysis*, 15(2), 124-139.
  [doi:10.1017/S1047198700006362](https://doi.org/10.1017/S1047198700006362)

- Pesaran, M. H., & Zhou, Q. (2018). Estimation of time-invariant effects in
  static panel data models. *Econometric Reviews*, 37(10), 1137-1171.
  [doi:10.1080/07474938.2016.1243552](https://doi.org/10.1080/07474938.2016.1243552)

- Greene, W. H. (2011). Fixed Effects Vector Decomposition: A Magical Solution
  to the Problem of Time-Invariant Variables in Fixed Effects Models?
  *Political Analysis*, 19(2), 135-146.
  [doi:10.1093/pan/mpq034](https://doi.org/10.1093/pan/mpq034)

## License

GPL-3
