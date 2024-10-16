GxEScanR: Performs GWAS/GWEIS scans using BinaryDosage files
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

# GxEScanR

GxEScanR is designed to efficiently run genome-wide association study
(GWAS) and genome-wide by environmental interaction study (GWEIS) scans
using imputed genotypes stored in the BinaryDosage format. The phenotype
to be analyzed can either be a continuous or binary trait. The GWEIS
scan performs multiple tests that can be used in two-step methods.

### Development version

The development version has an improved convergence routine. It is not
required but it is recommended to use the latest development version of
the BinaryDosage package. Click
[here](https://github.com/USCbiostats/BinaryDosage/tree/htslib) for
instructions

## Installation

To install the development version from [GitHub](https://github.com/)
run the following code.

``` r
# Uncomment the following line if devtools is not installed.
# install.packages("devtools")
devtools::install_github("USCbiostats/GxEScanR@develop")
```

## Example

### Notation and models

There are five models that can be fit using GxEScanR when running a
GWEIS and one that has to be fit before running the GWEIS.

$Y$ represents the outcome or phenotype

$X$ represents the confounding covariates

$E$ represents the covariate that will have its genetic interaction
tested

$G$ represents the genetic data

$\beta$ represents the effect estimate that will be fitted in the model

#### Model 0 - Environment-only Model

This model contains the phenotype, confounding covariates, and the
interaction covariate. This must be fit prior to running the GWEIS using
the glm() routine.

$Y = \beta_X X + \beta_E E$

#### Model 1 - Gene-only Model

This models contains all the values as the environment only model and
adds the genetic data

$Y = \beta_X X + \beta_E E + \beta_G G$

$\beta_G$ is the only term tested in this model and is represented by
bg_ge in the output.

#### Model 2 - GxE Interaction Model

This models contains all the values as the gene-only model and adds the
gene-environment interaction term.

$Y = \beta_X X + \beta_E E + \beta_G G + \beta_{GxE} GE$

$\beta_G$ and $\beta_{GxE}$ are both tested for this model as well as a
joint test for both terms. They are represented by bg_gxe, bgxe, and
joint, respectively in the output.

#### Model 3 - EG model

This models tests for a relationship between the environmental value of
interest and the gene. If $E$ is coded 0,1, a logistic model is fitted,
otherwise a linear model is fitted.

$E = \beta_X X + \beta_G G$

$\beta_G$ is tested in this model and is represented by bg_eg in the
output.

#### Model 4 - Case-only Model

This model is the same as model 3 except it is run on the cases only.

$\beta_G$ is tested and is represented by bg_case in the output.

This model can only be run on a dichotomous phenotype.

#### Model 5 - Control-only Model

This model is the same as model 3 except it is run on the controls only.

$\beta_G$ is tested and is represented by bg_ctrl in the output.

This model can only be run on a dichotomous phenotype.

### Steps needed to run a GWEIS

There are 5 steps to running a GWEIS using GxEScanR

1.  Preparing the Data
    - Genetic Data
    - Subject Data
2.  Fitting the Base Model
3.  Allocating Memory for the GWEIS
4.  Running the GWEIS
5.  Processing the results

### Preparing the Data

#### Genetic Data

GxEScanR using data that is stored in the BinaryDosage format. The
BinaryDosage package converts files from the VCF format to the
BinaryDosage format. It is recommended that the latest version of the
BinaryDosage package be used. It is available on github
[here](https://github.com/USCbiostats/BinaryDosage/tree/htslib).

The following code is for the latest version of BinaryDosage.

This example uses a VCF file included with the GxEScanR package. The
complete file name is obtained using the system.file() routine. The file
is then converted into a BinaryDosage formatted file using the vcftobd()
routine. This requires two file names, one for the BinaryDosage
formatted file and a BinaryDosage information file. The first file
contains the genetic data. The second file contains information on how
to read the file.

The BinaryDosage information file is saved as an rds file, R dataset,
and can be read using the readRDS() routine.

``` r
library(GxEScanR)

# Get the file name for the VCF file.
vcffile <- system.file("extdata", "gendata.vcf.gz", package = "GxEScanR")

exampledir <- tempdir()
bdosefile <- paste(exampledir, "\\gendata.bdose", sep = "")
bdinfofile <- paste(exampledir, "\\gendata.bdinfo", sep = "")

bdinfo <- vcftobd(vcffiles = vcffile, TRUE, c(bdosefile, bdinfofile))
```

#### Subject Data

The subject data consists of the subject IDs that will be used to link
the genetic data, and the phenotype and the covariates that will be used
in the GWEIS. The data included with the GxEScanR package has two
phenotypes, one linear and one dichotomous. The data also contains two
covariates, x1 and x2.

**Important** All values for phenotypes and covariates must be numeric,
factors are not allowed.

**Important** All subject IDs must be character values.

**Important** When using a dichotomous phenotype, it must be coded 0,1.

``` r
library(GxEScanR)

subjectfile <- system.file("extdata", "subdata.rds", package = "GxEScanR")

subjectdata <- readRDS(subjectfile)

head(subjectdata)
#>      subid y_linear y_logistic x1 x2
#> 1 subject1 7.939120          1  0  0
#> 2 subject2 5.965786          0  1  0
#> 3 subject3 3.047969          0  1  0
#> 4 subject4 1.173026          0  0  0
#> 5 subject5 7.691592          1  1  0
#> 6 subject6 3.931234          0  1  0
```

### Fitting the Base Model

The base model is fit using the glm routine using data that contains
only the subject IDs, phenotype, and covariates of interest. No genetic
data will be in the model.

It is important to remove subjects that have missing data and subjects
that don’t have genetic data. The subjects that have genetic are listed
in the BinaryDosage information file. It is also important to keep the
subject IDs in the data used by the glm() routine because this is used
to match the genetic data to the correct phenotype and covariates.

The formula used in the glm() routine must have the covariate that will
have its interaction with the genetic data tested be listed last. In
this example we will be testing for a genetic interaction with the x1
phenotype. There the formula value will be

y ~ x2 + x1

Note: It is recommended to look at the results from the glm() routine to
verify the model was fitted correctly.

#### Continuous Phenotype

``` r
library(GxEScanR)

subjectfile <- system.file("extdata", "subdata.rds", package = "GxEScanR")

subjectdata <- readRDS(subjectfile)

# Remove y_logistic from the subject data
lineardata <- subjectdata[, c(1, 2, 4, 5)]

# Keep only subjects with complete data
lineardata <- lineardata[complete.cases(lineardata), ]

# Keep only subjects with genetic data
lineardata <- lineardata[!is.na(match(bdinfo$samples$sid, lineardata$subid)), ]

# Fit the linear model using the glm routine
linearmodel <- glm(y_linear ~ x2 + x1, data = lineardata)
```

#### Dichotomous Phenotype

``` r
library(GxEScanR)

subjectfile <- system.file("extdata", "subdata.rds", package = "GxEScanR")

subjectdata <- readRDS(subjectfile)

# Remove y_linear from the subject data
logisticdata <- subjectdata[, c(1, 3, 4, 5)]

# Keep only subjects with complete data
logisticdata <- logisticdata[complete.cases(logisticdata), ]

# Keep only subjects with genetic data
logisticdata <- logisticdata[!is.na(match(bdinfo$samples$sid, logisticdata$subid)),
    ]

# Fit the logistic model using the glm routine
logisticmodel <- glm(y_logistic ~ x2 + x1, "binomial", data = logisticdata)
```

### Allocating Memory for the GWEIS

Memory needs to be allocated for the GWEIS. This makes the GWEIS run
much quicker. Allocating memory for the GWEIS is done using the
gweis.mem() in GxEScanR passing it the model fitted without the genetic
data, the subject IDs, and the list of GWEIS tests to perform. The
subject IDs should come from the dataset that was used in the glm()
routine.

The following tests are the ones normally performed in a GWEIS with a
continuous phenotype, “bg_ge”, “bg_gxe”, “bgxe”, and “joint”

The following tests are the ones normally performed in a GWEIS with a
dichotomous phenotype, “bg_ge”, “bg_gxe”, “bgxe”, “joint”, “bg_eg”,
“bg_case”, and “bg_ctrl”

``` r
library(GxEScanR)

# GWEIS with continuous outcome
lineartests <- c("bg_ge", "bg_gxe", "bgxe", "joint")
linearmem <- gweis.mem(gemdl = linearmodel, subids = lineardata$subid, tests = lineartests)

# GWEIS with dichotomous outcome
logistictests <- c("bg_ge", "bg_gxe", "bgxe", "joint", "bg_eg", "bg_case", "bg_ctrl")
logisticmem <- gweis.mem(gemdl = logisticmodel, subids = logisticdata$subid, tests = logistictests)
```

### Running the GWEIS

To run the GWEIS you will need the BinaryDosage information file, the
list returned from the gweis.mem() routine, a list of SNPs to test, and
the name of the file to save the results.

The list of SNPs to perform the GWEIS on can be either the SNP IDs or
the indices in the BinaryDosage information. In this example all the
SNPs in the genetic data will be tested and they will be identified my
index.

``` r
library(GxEScanR)

# GWEIS with continuous outcome
snpindex <- 1:nrow(bdinfo$snps)
linearresults <- paste(exampledir, "linear.txt", sep = "//")
res <- rungweis(gweismem = linearmem, bdinfo = bdinfo, snps = snpindex, outfilename = linearresults)

# GWEIS with dichotomous outcome
snpindex <- 1:nrow(bdinfo$snps)
logisticresults <- paste(exampledir, "logistic.txt", sep = "//")
res <- rungweis(gweismem = logisticmem, bdinfo = bdinfo, snps = snpindex, outfilename = logisticresults)
```

### Processing the results

Processing the results is done mostly with other tools. The following
talks about what is in the output files. Only the first 3 lines of the
results will be displayed for this example.

``` r
# GWEIS with continuous outcome
lineardf <- read.table(linearresults, header = TRUE, sep = "\t")
# GWEIS with continuous outcome
logisticdf <- read.table(logisticresults, header = TRUE, sep = "\t")
# print(logisticdf)
```

For a continuous outcome only the first 2 models listed above are
fitted. The SNP data consists of SNP ID (snpid), chromosome number
(chr), location (loc), reference allele (ref), alternate allele (alt),
and alternate allele frequency (aaf).

``` r
knitr::kable(lineardf[1:3, 1:6])
```

| snpid              | chr   |      loc | ref | alt |       aaf |
|:-------------------|:------|---------:|:----|:----|----------:|
| chr22:10527916:T:G | chr22 | 10527916 | T   | G   | 0.6494005 |
| chr22:10648779:G:T | chr22 | 10648779 | G   | T   | 0.0854975 |
| chr22:10648794:G:A | chr22 | 10648794 | G   | A   | 0.0536727 |

For a dichotomous outcome, two additional columns are added, the
alternate allele frequency in cases (aaf_case) and the alternate allele
frequency in controls (aaf_ctrl).

``` r
knitr::kable(logisticdf[1:3, 1:8])
```

| snpid              | chr   |      loc | ref | alt |       aaf |  aaf_case |  aaf_ctrl |
|:-------------------|:------|---------:|:----|:----|----------:|----------:|----------:|
| chr22:10527916:T:G | chr22 | 10527916 | T   | G   | 0.6494005 | 0.6264272 | 0.6540861 |
| chr22:10648779:G:T | chr22 | 10648779 | G   | T   | 0.0854975 | 0.0882864 | 0.0849287 |
| chr22:10648794:G:A | chr22 | 10648794 | G   | A   | 0.0536727 | 0.0457573 | 0.0552871 |

For model 1 there are two value output, bg_ge and bg_ge_lrt. The first
value is the $\beta_G$ estimate in model 1. The second value is the
likelihood ratio test (LRT) 1 df chi-squared value for $\beta_G = 0$.
These values are output when “bg_ge” is included in the tests value
passed to the gweis.mem() routine.

``` r
knitr::kable(lineardf[1:3, 7:8])
```

|      bg_ge | bg_ge_lrt |
|-----------:|----------:|
| -0.8900621 | 5.2558824 |
| -0.4766408 | 0.6737040 |
| -0.6481754 | 0.8070534 |

For model 2 there are five possible values that can be output.

bg_gxe and bg_gxe_lrt are the estimate of $\beta_G$ and the LRT 1 df
test that the valueof the estimate is zero. These values are output when
“bg_gxe” is included in the tests value passed to the gweis.mem()
routine.

bgxe and bgxe_lrt are the estimate of $\beta_{GxE}$ and the LRT 1 df
test that the value of the estimate is zero. These values are output
when “bgxe” is included in the tests value passed to the gweis.mem()
routine.

joint_lrt is the LRT 2 df chi-squared value for the test that
$\beta_G = 0$ and $\beta_{GxE} = 0$. This value is output when “joint”
is included in the tests value passed to the gweis.mem() routine.

``` r
knitr::kable(lineardf[1:3, 9:13])
```

|     bg_gxe | bg_gxe_lrt |       bgxe |  bgxe_lrt | joint_lrt |
|-----------:|-----------:|-----------:|----------:|----------:|
| -1.4864649 |  8.4363261 |  1.4016489 | 3.2015268 | 8.4574092 |
| -0.4591966 |  0.3719000 | -0.0430555 | 0.0013243 | 0.6750284 |
| -0.6509083 |  0.5512644 |  0.0084781 | 0.0000301 | 0.8070835 |

For models 3 through 5 the estimates for $\beta_G$ along with the LRT 1
df chi-squared value for $\beta_G = 0$. The estimates for $\beta_G$ are
labeled bg_eg, bg_case, and bg_ctrl for the 3 tests, respectively. The
test statistics are labeled bg_eg_lrt, bg_case_lrt, and bg_ctrl,
respectively.

``` r
knitr::kable(logisticdf[1:3, 16:21])
```

|      bg_eg | bg_eg_lrt |    bg_case | bg_case_lrt |    bg_ctrl | bg_ctrl_lrt |
|-----------:|----------:|-----------:|------------:|-----------:|------------:|
| -0.3460221 | 0.7862570 |  0.8264053 |   0.8541469 | -0.4673483 |   1.1215720 |
| -0.0867130 | 0.0219140 | -0.7545895 |   0.4064202 |  0.0815513 |   0.0144966 |
| -0.2255839 | 0.0949399 | -0.4952261 |   0.0361789 | -0.0220965 |   0.0008210 |
