CRAN comments
================

## Test environments
* local OS X install, R 4.0.2
* local Windows 10 install, R 4.0.2
* ubuntu 16.04.6 (on travis-ci), R-release
* win-builder (devel, release)

## R CMD check results
1 note on all systems
* checking for future file timestamps ... NOTE
unable to verify current time

1 warning on local Windows 10
* checking compiled code ... OK
   WARNING
  'qpdf' is needed for checks on size reduction of PDFs
This appears to be a local system setup issue.

## rhub results
2 Notes
* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'John Morrison <jmorr@usc.edu>'
  
  
  Possibly mis-spelled words in DESCRIPTION:
    BinaryDosage (26:3)
    GWAS (3:14, 23:31)
    GWIS (3:19, 23:36)
    getbdinfo (26:17)
  New submission

Words are correct. GWAS and GWIS are acronyms.
BinaryDosage is an R package and getbdinfo is
a routine in that package.

* checking for future file timestamps ... NOTE
  unable to verify current time

## check_win_xxx results
1 Notes
* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'John Morrison <jmorr@usc.edu>'
  
  
  Possibly mis-spelled words in DESCRIPTION:
    BinaryDosage (26:3)
    GWAS (3:14, 23:31)
    GWIS (3:19, 23:36)
    getbdinfo (26:17)
  New submission

Words are correct. GWAS and GWIS are acronyms.
BinaryDosage is an R package and getbdinfo is
a routine in that package.

