GxEScanR: Performs GWAS/GWIS analyses using BinaryDosage files
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

# GxEScanR

<!-- badges: start -->

[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/USCbiostats/GxEScanR?branch=master&svg=true)](https://ci.appveyor.com/project/USCbiostats/GxEScanR)
[![Travis build
status](https://travis-ci.com/USCbiostats/GxEScanR.svg?branch=master)](https://travis-ci.com/USCbiostats/GxEScanR)
[![Codecov test
coverage](https://codecov.io/gh/USCbiostats/GxEScanR/branch/master/graph/badge.svg)](https://codecov.io/gh/USCbiostats/GxEScanR?branch=master)
<!-- badges: end -->

GxEScanR is designed to efficiently run GWAS/GWIS scans using imputed
genotypes stored in the BinaryDosage format. The phenotype to be
analyzed can either be a continuous or binary trait. The GWIS scan
performs multiple tests that can be used in two-step methods.

## Installation

You can install the released version of GxEScanR from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("GxEScanR")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("USCbiostats/GxEScanR")
```

## Example

The following is a step by step example on how to run a GWAS.

The first step is to load the subject phenotype and covariate data. In
this example the data is included with the package. The data is stored
in a RDS formatted file. The first five lines of the data frame are
shown. The first column is the subject ID, sid. The second column is the
phenotype, y. The last column is a covariate,e.

``` r
library(GxEScanR)

covdatafile <- system.file("extdata", "covdata.rds", package = "GxEScanR")
covdata <- readRDS(covdatafile)
covdata[1:5,]
#>   sid y e
#> 1  I1 0 0
#> 2  I2 0 0
#> 3  I3 0 0
#> 4  I4 0 0
#> 5  I5 0 0
```

The second step is to load the information about the binary dosage file.
This is obtained by running
BinaryDosage::getbdinfo(<binary dosage file name>). More information
about this can be found in the BinaryDosage package.

A binary dosage file and information file about it are included with
this package. The getbdinfo routine stores the complete file path to the
binary dosage file. The installation routine moved the binary dosage
from its original location. The third line of code corrects this. The
user will not need to run the third line in normal usage.

``` r
  bdinfofile <- system.file("extdata", "pdata_4_1.bdinfo", package = "GxEScanR")
  bdinfo <- readRDS(bdinfofile)
  # Not normally run - This is needed only for the example data file
  bdinfo$filename <- system.file("extdata", "pdata_4_1.bdose", package = "GxEScanR")
```

Everything is now ready to run a GWAS. The number of subjects used in
the analysis is displayed at the start of the run. There are a lot of
options not being shown in this README file. Information on these
options can be found in the documentation and vignettes.

``` r
results <- gwas(data = covdata, bdinfo = bdinfo)
#> [1] "200 subjects have complete data"
results
#>       snp      betag      lrtg
#> 1 1:10001  0.7129412 9.0230000
#> 2 1:10002 -0.0943161 0.1219431
#> 3 1:10003 -0.3143876 1.3776592
#> 4 1:10004  0.1397320 0.2356375
#> 5 1:10005  0.2680398 1.2002662
```

In this example the output was stored in a data frame and displayed. An
option exists to output the results to a text file that can easily be
read into R using read.table.

The columns in the output are the SNP ID, the coefficient estimate for
the genetic effect (betag), and the likelihood ratio test for the
estimate (lrtg). The first SNP was simulated with a log odds ratio of
0.75, and the others were simulated to have no effect. The results are
consistent with the modelling.
