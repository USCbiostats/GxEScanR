CRAN comments
================

## Test environments
* local OS X install, R 4.0.2
* local Windows 10 install, R 4.0.2
* ubuntu 16.04.6 (on travis-ci), R-release
* win-builder (devel, release, and oldrelease)

## R CMD check results
One note on all systems
* checking for future file timestamps ... NOTE
unable to verify current time
One warning on local Windows 10
> checking compiled code ... OK
   WARNING
  'qpdf' is needed for checks on size reduction of PDFs
  
## rhub results
2 Notes
   Possibly mis-spelled words in DESCRIPTION:
     VCF (21:56)
N  checking for non-standard things in the check directory
   Found the following files/directories:
     'BinaryDosage-Ex_i386.Rout' 'BinaryDosage-Ex_x64.Rout'
     'examples_i386' 'examples_x64' 'tests_i386' 'tests_x64'
VCF is correct
Second note appears to be issue with rhub

## check_win_xxx results
1 Notes
   Possibly mis-spelled words in DESCRIPTION:
     VCF (21:56)
VCF is correct.

