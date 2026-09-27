# xtfifevd 1.0.2

* Resubmission to CRAN.
* Removed UTF-8 BOM (Byte Order Mark) from `inst/CITATION`. The BOM was
  causing the CRAN incoming auto-check to fail with:
  `Invalid citation information in 'inst/CITATION': 1:1: unexpected input`.
* Synchronised the version string in `inst/CITATION` with `DESCRIPTION`.
* No code or behavior changes.

# xtfifevd 1.0.1

* Documentation and metadata corrections (DOI fixes, `cat()` -> `message()`).

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
