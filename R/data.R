#' Simulated subject data set
#' 
#' A data set with subject IDs, an outcome, and two covariates
#' 
#' @format A data frame with 100 rows and 4 variables
#' \describe{
#'   \item{iid}{subject id}
#'   \item{y}{outcome variable}
#'   \item{x1}{first covariate}
#'   \item{x2}{second covariate}
#' }
#' @source Writers of GxEScanR
"simSubjectData"

#' Simulated genetic data set
#' 
#' A data set containing information about a binary dosage file
#' 
#' @format A list with 13 entries
#' \describe{
#'   \item{filetype}{type of genetic file}
#'   \item{filename}{name of file with genetic data}
#'   \item{format}{format of file}
#'   \item{versin}{version of the format of the file}
#'   \item{Groups}{number of groups in the data set}
#'   \item{GroupSizes}{vector of sizes of each group}
#'   \item{NumSamples}{number of samples in data set}
#'   \item{Samples}{data frame with family and subject id of each sample}
#'   \item{NumSNPs}{number of SNPs in data set}
#'   \item{SNPs}{data frame with SNP ID, chromosome, location, reference allele, and alternate allele}
#'   \item{SNPInfo}{data frame with alternate allele frequency, minor allele frequency, average call, and r-squared}
#'   \item{Indices}{vector of integers used to locate SNP data in file}
#' }
#' @source Writers of GxEScanR
"simGeneticData"
