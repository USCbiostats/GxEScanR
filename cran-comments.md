## Resubmission

This is a resubmission. The package was previously archived because of the
following NOTE:

> Specified C++11: please drop specification unless essential.

All C++ source code has been moved to the `lsReg` package. GxEScanR now depends
on `lsReg` and `BinaryDosage` and contains no compiled code of its own,
eliminating the C++11 specification entirely.

Additional changes since the archived version (2.0.2):

* New API: `gweis.mem()` and `rungweis()` replace the former `gwas()` and
  `gweis()` functions, supporting a wider set of GWEIS tests (gene-only,
  GxE interaction, joint, E|G, case-only, control-only).
* Test suite added (testthat).
* Vignette added.

## Test environments

* Windows 11, R 4.5.3 (local)
* GitHub Actions, ubuntu-latest, R release
* GitHub Actions, ubuntu-latest, R devel
* GitHub Actions, macOS-latest, R release
* GitHub Actions, windows-latest, R release

## R CMD check results

0 ERRORs | 0 WARNINGs | 0 NOTEs

(Local check produces 1 NOTE — "unable to verify current time" — due to no
internet access on the local machine. This does not appear in any CI
environment.)

## Reverse dependencies

There are no reverse dependencies on the archived version of this package.
