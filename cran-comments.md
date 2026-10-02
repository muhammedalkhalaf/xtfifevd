## xtfifevd 1.1.0

This release corrects the computations below; the version on CRAN is 1.0.2.

* Formula handling: the right-hand side is now split at the top-level `|`
  and each part is built with `model.frame()` and `model.matrix()`, so
  transformations such as `log(x)`, `I(x^2)`, `x:w` and `log(y)` and factors
  (with their contrasts) are handled as in `lm()`. Previously `all.vars()`
  was used and the raw columns were read, so transformations were silently
  dropped or split. Intercept columns are removed from both parts. A factor
  `id` with unused levels now works (levels are dropped). Duplicated
  `(id, time)` pairs are now an error.
* Variance of the time-varying coefficients `beta`: the reported block of
  `vcov()` was the homoskedastic `sigma2_e (X'MX)^-1`, while the variance of
  `gamma` used the panel-robust matrix of Pesaran and Zhou (2018, eq. 18).
  The robust matrix is now reported for `beta` everywhere. New argument
  `vcov_beta = c("robust", "classical")` in `xtfifevd()`, `fevd()`, `fef()`
  and `fef_iv()`; the chosen matrix is also used inside the Pesaran and Zhou
  variance of `gamma` and of the intercept. Both matrices are returned in
  `V_beta_robust` and `V_beta_classical`.
* Full covariance matrix: `vcov()` was block diagonal. It now contains the
  covariance between `gamma` and `beta` from Pesaran and Zhou (2018, eq.
  A.11), `Cov(gamma, beta) = -Qzz^-1 Qzxbar Var(beta)`, and the variance of
  and covariances with the intercept `alpha = ubar - zbar' gamma` (eq. 5) by
  the delta method. The same derivation is applied to FEF-IV with the
  instrument projections (eq. 48 and 51). The previous intercept standard
  errors (stage 2 OLS for FEF, `s^2 ginv([1, R]'[1, R])[1, 1]` for FEF-IV)
  are no longer used.
* `fevd()` now runs the three stages of Plumper and Troeger (2007)
  literally: within FE regression, unit-level regression of the
  time-averaged FE residuals on an intercept and z, and pooled OLS of y on
  an intercept, x, z and the unexplained unit effect h_i. Previously no
  stage 3 was run and `fevd()` returned the FEF results. The stage 3
  coefficients are returned in `stage3$coefficients`, the coefficient on
  h_i in `delta` (equal to 1 by construction, Pesaran and Zhou 2018,
  Proposition 3, also in unbalanced panels), and the naive stage 3 standard
  errors in `stage3$se_naive` for reference only. Inference continues to use
  the Pesaran and Zhou standard errors; the documentation states that the
  naive stage 3 standard errors are too small for the time-invariant
  coefficients (Breusch, Ward, Nguyen and Kompas 2010, Theorem 3; Greene
  2011) and why an intercept is included in
  stage 2. Stage 2 is unit-level and unweighted in unbalanced panels;
  Plumper and Troeger are silent on unbalanced panels.
* Rarely changing variables after `|`: the unit (panel) mean is now used in
  all stages (previously the first-period value was used silently); the
  warning is kept and documented.
* `bw_ratio()`: the print header now says "SD Ratios" (Plumper and Troeger
  2007 define the b/w ratio as the between SD over the within SD). The
  guidance no longer states a single 1.7 threshold: the thresholds of
  Plumper and Troeger (2007, Fig. 4; N = 30, T = 20) are about 0.2, 1.7,
  2.8 and 3.8 for corr(z, u) of 0, 0.3, 0.5 and 0.8, the correlation is not
  observable or testable, the coefficient is biased whenever it is non-zero,
  and Plumper and Troeger offer no simple rule of thumb.
* `sigma2_u` is now the variance of the unexplained unit effect (stage 2
  residual variance minus `sigma2_e * mean(1 / T_i)`, truncated at zero).
  Previously it was the variance of the time-averaged FE residuals, which
  includes `z'gamma` and `sigma2_e / T`. The summary label says so.
* `print(summary())` now shows the installed package version (previously a
  hard-coded "1.0.0") and cites Pesaran and Zhou (2018) (previously 2016).
* `MASS` is no longer imported; `utils` is.
* New tests pin the FEF and FEF-IV coefficients and the full covariance
  matrix to hand computations of the Pesaran and Zhou formulas, the formula
  parser to `lm()` on the transformed data, the robust `beta` covariance to
  `plm::vcovHC(method = "arellano", type = "HC0")` (skipped when `plm` is not
  available), and `delta = 1`.

## Test environments

* Ubuntu 24.04, R 4.3.3 and R-devel, R CMD check --as-cran

## R CMD check results

0 errors | 0 warnings | 0 notes
