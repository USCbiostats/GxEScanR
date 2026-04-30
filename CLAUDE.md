# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Package Does

GxEScanR is an R package for high-efficiency genome-wide association studies (GWAS) and genome-wide by environment interaction studies (GWEIS) on **BinaryDosage**-formatted imputed genetic data. It supports both continuous (linear regression) and binary (logistic regression) traits.

Main exported functions: `gwas()` and `gweis()`.

## Build and Development Commands

```r
devtools::document()   # Regenerate NAMESPACE and Rd files from roxygen2 comments
devtools::build()      # Build the package tarball
devtools::check()      # Run R CMD check (includes tests)
devtools::test()       # Run tests only
devtools::install()    # Install locally
```

From the shell:
```bash
R CMD build GxEScanR
R CMD check GxEScanR_*.tar.gz
R CMD INSTALL GxEScanR
```

Run a single test file:
```r
testthat::test_file("tests/testthat/test-results.R")
```

## Architecture

### Data Flow

1. User supplies a covariate/phenotype `data.frame` + a `bdinfo` object (BinaryDosage file metadata).
2. `gwas()`/`gweis()` validates inputs (`validateinput()`), subsets to complete cases (`subsetdata()`), and optionally filters SNPs (`subsetsnps()`).
3. A baseline regression model is fit to the null (no genetic effect).
4. `assignblocks()` calculates adaptive block sizes (smaller blocks for larger sample counts) to manage memory.
5. Each block: binary dosage data is read from disk (`bdread.cpp` → `readblock()`/`getdosages()`), then regression models are updated via QR decomposition rank-1 updates.
6. Results are returned as a `data.frame` or written to a tab-delimited file.

### Source Files

**R/**
- `GxEScanR.R` — Main entry points (`gwas`, `gweis`), input validation, and block management helpers.
- `LinReg.R` — R wrappers for large-scale linear regression (`linreggwas`, `linreggweis`).
- `LogReg.R` — R wrappers for large-scale logistic regression with Newton-Raphson scoring (`logreggwas`, `logreggweis`).
- `RcppExports.R` — Auto-generated Rcpp bindings; do not edit manually.

**src/**
- `LSLinReg.cpp` — C++ LS linear regression via QR decomposition (`initlslinreg`, `lslinreg`).
- `LSLogReg.cpp` — C++ logistic regression via iterative scoring (`initlslogreg`, `lslogreg`).
- `bdread.cpp` — Binary dosage file I/O, supports formats 1.1, 1.2, and extended (`readblock`, `getdosages`).
- `Makevars` / `Makevars.win` — Require **C++11**, **RcppArmadillo** (LAPACK/BLAS), and **OpenMP**.

### Key Dependencies

- **RcppArmadillo** — C++ linear algebra; requires LAPACK/BLAS at link time.
- **prodlim** — Used for `row.match()` in subject ID alignment.
- **BinaryDosage** — Soft dependency; provides the `bdinfo` objects consumed by this package.

### Output Columns

- **GWAS:** `snpid`, `betag` (genetic effect), `lrtg` (likelihood ratio test statistic)
- **GWEIS:** `betadg`, `lrtdg`, `betagxe` (GxE interaction coefficient), `lrtgxe`, `lrt2df` (2-degree-of-freedom test)

## Testing

Framework: **testthat** (>= 2.1.0). Test data lives in `inst/extdata/`.

- `tests/testthat/test-gwasinput.R` — Input validation coverage (bad arguments, ID handling, parameter edge cases).
- `tests/testthat/test-results.R` — Numerical correctness: results are compared against reference `.rds` files (`lingwas.rds`, `lingweis.rds`, `loggwas.rds`, `loggweis.rds`) with tolerance `< 1e-6`.

When adding C++ code, regenerate `RcppExports.R` and `RcppExports.cpp` with:
```r
Rcpp::compileAttributes()
```
