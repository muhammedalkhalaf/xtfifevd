# xtfifevd 1.0.0

* Initial CRAN release.

* Implements three estimation methods for time-invariant variables in panel FE models:
  - `fevd()`: Fixed Effects Vector Decomposition (Plümper & Troeger, 2007)
  - `fef()`: Fixed Effects Filtered (Pesaran & Zhou, 2018)
  - `fef_iv()`: FEF with instrumental variables (Pesaran & Zhou, 2018)

* Uses correct Pesaran-Zhou (2018) variance estimators that account for generated
  regressor uncertainty.

* Provides `bw_ratio()` diagnostic for between/within variance analysis.

* Full S3 methods: `print()`, `summary()`, `coef()`, `vcov()`, `confint()`.
