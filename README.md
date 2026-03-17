# xtfifevd

**Fixed Effects Filtered and Vector Decomposition for Panel Data**

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/xtfifevd)](https://CRAN.R-project.org/package=xtfifevd)
<!-- badges: end -->

## Overview

`xtfifevd` estimates panel data models with time-invariant regressors. It
implements three methods:

- **FEVD**: Fixed Effects Vector Decomposition (Plumper & Troeger, 2007)
- **FEF**: Fixed Effects Filtered (Pesaran & Zhou, 2016)
- **FEF-IV**: FEF with instruments for endogenous time-invariant regressors

All methods use the **Pesaran-Zhou corrected variance estimator** (2016),
which is consistent. The raw FEVD Stage-3 standard errors are inconsistent
(Remark 4 in Pesaran & Zhou 2016) and are reported for diagnostic comparison only.

## Installation

```r
install.packages("xtfifevd")
```

## Usage

```r
library(xtfifevd)

set.seed(42)
n <- 10; tt <- 15
dat <- data.frame(
  id   = rep(1:n, each = tt),
  time = rep(1:tt, times = n),
  y    = rnorm(n * tt),
  x    = rnorm(n * tt),
  z    = rep(rnorm(n), each = tt)   # time-invariant
)

res <- xtfifevd(y ~ x, data = dat, index = c("id", "time"),
                zinvariants = "z", method = "fevd")
print(res)
summary(res)
```

## References

Plumper, T. and Troeger, V.E. (2007). Efficient Estimation of Time-Invariant
and Rarely Changing Variables in Finite Sample Panel Analyses with Unit Fixed
Effects. *Political Analysis*, 15(2), 124–139.
<https://doi.org/10.1093/pan/mpm002>

Pesaran, M.H. and Zhou, Q. (2016). Estimation of Time-Invariant Effects in
Static Panel Data Models. *Econometric Reviews*, 37(10), 1137–1171.
<https://doi.org/10.1080/07474938.2015.1032164>
