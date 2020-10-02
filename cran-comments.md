CRAN comments
================

## Test environments
* local OS X install, R 4.0.2
* local Windows 10 install, R 4.0.2
* ubuntu 16.04.6 (on travis-ci), R-release
* win-builder (oldrelease, devel, release)

## R CMD check results
1 warning on local Windows 10
* checking compiled code ... OK
   WARNING
  'qpdf' is needed for checks on size reduction of PDFs
This appears to be a local system setup issue.

## rhub results
1 Note
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'John Morrison <jmorr@usc.edu>'

Days since last update: 3

This is a requested update to fix to an error identified in testing

## check_win_xxx results
1 Notes
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'John Morrison <jmorr@usc.edu>'

Days since last update: 3

This is a requested update to fix to an error identified in testing

## Comments from second submission

The second submission was accepted. After acceptance I received an email indicating there was an installation error on an r-patched-solaris-x86 system and that the error needed to be fixed.

## Response to comments

There was a type-casting issue in the C++ code that only appears
when using the solaris compiler. A type-cast has been added that will
address the issue. No other code was modified. The NEWS.md file
was modified to indicate the update.

## Comments from first submission

Please do not start the description with "This package", package name,
title or similar.

Please always explain all acronyms in the description text.

Please add () behind all function names in the description texts
(DESCRIPTION file). e.g: --> getbdinfo()

## Responses to comments

Updated description.

Updated all documentation to have meanings of acronyms used.

Added () behind function name getbdinfo.

