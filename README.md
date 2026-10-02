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

- **FEVD**: Fixed Effects Vector Decomposition (Plumper and Troeger, 2007),
  three stages with an intercept in stage 2
- **FEF**: Fixed Effects Filtered (Pesaran and Zhou, 2018)
- **FEF-IV**: FEF with instrumental variables for endogenous regressors

All methods use the **Pesaran and Zhou (2018) variance estimators**, which
account for generated regressor uncertainty (the naive FEVD stage 3 standard
errors are too small for the time-invariant coefficients), and report the full covariance matrix of the
time-varying coefficients, the time-invariant coefficients and the intercept.

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
# Transformations and factors are allowed, e.g. log(y) ~ x + I(x^2) | z
fit <- xtfifevd(y ~ x | z, data = data, id = "id", time = "time")
summary(fit)
fit$delta   # FEVD stage 3 coefficient on h_i, equal to 1 by construction
```

Output:
```
======================================================================
FEVD Estimation Results                            xtfifevd 1.1.0
======================================================================
Dep. variable:   y
Method:          FEVD
Variance:        Pesaran and Zhou (2018), beta vcov: robust
Observations:    1000       Groups:     100
T (average):     10.00
----------------------------------------------------------------------

      Estimate Std. Error z value Pr(>|z|)
x      2.03092    0.01619 125.469  < 2e-16 ***
z      0.42660    0.08702   4.902 9.48e-07 ***
_cons  1.08528    0.09019  12.033  < 2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

----------------------------------------------------------------------
Time-varying (FE):      x
Time-invariant:         z
sigma_e: 0.4881    sigma_u (unexplained unit effect): 0.9046
FEVD stage 3 coefficient on h_i (delta): 1.000000  [equals 1 by construction]
Naive stage 3 OLS SEs are too small for the time-invariant
coefficients; see ?xtfifevd.
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
# Between/Within SD ratio (Plumper and Troeger 2007 define it as the
# between SD divided by the within SD)
bw_ratio(data, c("z", "x"), id = "id")
```

Plumper and Troeger (2007, Fig. 4; N = 30, T = 20) find that FEVD has lower
RMSE than FE for a rarely changing variable when its b/w ratio exceeds about
0.2 if corr(z, u) = 0, 1.7 at corr 0.3, 2.8 at 0.5 and 3.8 at 0.8. The
correlation with the unit effects is not observable or testable, the FEVD
and FEF coefficients are biased whenever it is non-zero, and Plumper and
Troeger state that they "cannot offer a simple rule of thumb".

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
2. **Naive FEVD Stage 3 SEs**: too small for the time-invariant coefficients
   (Breusch et al. 2010, Theorem 3; Greene 2011); in a Monte Carlo check
   (400 replications, N = 200, T = 8, AR(1) errors with coefficient 0.8,
   heteroskedastic across units) their 95 percent coverage was 35 to 40
   percent against 92 to 95 percent for the Pesaran and Zhou SEs
3. **Ignoring the problem**: Biased pooled OLS

### The Right Solution

FEVD/FEF methods with **Pesaran-Zhou corrected standard errors** provide:

- Consistent point estimates under standard FE assumptions
- Valid inference that accounts for generated regressor uncertainty
- Better efficiency than FE for variables with high between/within ratio

## References

- Plumper, T. and Troeger, V. E. (2007). Efficient Estimation of Time-Invariant
  and Rarely Changing Variables in Finite Sample Panel Analyses with Unit Fixed
  Effects. *Political Analysis*, 15(2), 124-139.
  [doi:10.1093/pan/mpm002](https://doi.org/10.1093/pan/mpm002)

- Pesaran, M. H. and Zhou, Q. (2018). Estimation of time-invariant effects in
  static panel data models. *Econometric Reviews*, 37(10), 1137-1171.
  [doi:10.1080/07474938.2016.1222225](https://doi.org/10.1080/07474938.2016.1222225)

- Breusch, T., Ward, M. B., Nguyen, H. and Kompas, T. (2010). On the
  fixed-effects vector decomposition. MPRA Paper No. 21452.
  https://mpra.ub.uni-muenchen.de/21452/

- Greene, W. H. (2011). Fixed Effects Vector Decomposition: A Magical Solution
  to the Problem of Time-Invariant Variables in Fixed Effects Models?
  *Political Analysis*, 19(2), 135-146.
  [doi:10.1093/pan/mpq034](https://doi.org/10.1093/pan/mpq034)

## License

GPL-3
