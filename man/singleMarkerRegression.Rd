% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polyqtlR_functions.R
\name{singleMarkerRegression}
\alias{singleMarkerRegression}
\title{Run a single marker regression using marker dosages}
\usage{
singleMarkerRegression(
  dosage_matrix,
  Phenotype.df,
  genotype.ID,
  trait.ID,
  maplist = NULL,
  perm_test = FALSE,
  N_perm = 1000,
  alpha = 0.05,
  ncores = 1,
  return_R2 = FALSE,
  log = NULL
)
}
\arguments{
\item{dosage_matrix}{An integer matrix with markers in rows and individuals in columns.
All markers in this matrix will be tested for association with the trait.}

\item{Phenotype.df}{A data.frame containing phenotypic values}

\item{genotype.ID}{The colname of \code{Phenotype.df} that contains the population identifiers (F1 names) (must be a colname of \code{Phenotype.df})}

\item{trait.ID}{The colname of \code{Phenotype.df} that contains the response variable to use in the model (must be a colname of \code{Phenotype.df})}

\item{maplist}{Option to include linkage map in the format returned by \code{MDSMap_from_list} from \code{polymapR}. If maplist is
not specified (by default \code{NULL}) then no ordering of markers from dosage-matrix is performed. Note that all markers in dosage_matrix are tested;
markers with dosages that were not on the maplist will be assigned unordered to linkage group 0 with dummy cM positions 1,2,3 etc.}

\item{perm_test}{Logical, by default \code{FALSE}. If \code{TRUE}, a permutation test will be performed to determine a
genome-wide significance threshold.}

\item{N_perm}{Integer. The number of permutations to run if \code{perm_test} is \code{TRUE}; by default this is 1000.}

\item{alpha}{Numeric. The P-value to be used in the selection of a threshold if \code{perm_test} is \code{TRUE};
by default 0.05 (i.e. the 0.95 quantile).}

\item{ncores}{Number of cores to use if parallel processing required. Works both for Windows and UNIX (using \code{doParallel}).
Use \code{parallel::detectCores()} to find out how many cores you have available.}

\item{return_R2}{Should the (adjusted) R2 of the model fit also be determined?}

\item{log}{Character string specifying the log filename to which standard output should be written. If \code{NULL} log is send to stdout.}
}
\value{
A list containing the following components:
\describe{
\item{QTL.res}{ The -log(p) of the model fit per marker are returned as "LOD" scores, although "LOP" would have been a better description.
If requested, R2 values are also returned in column "R2adj"}
\item{Perm.res}{ The results of the permutation test if performed, otherwise \code{NULL}}
\item{Map}{ The linkage map if provided, otherwise \code{NULL}}
\item{LG_names}{Names of the linkage groups, if a map was provided, otherwise \code{NULL}}
}
}
\description{
Function to run a single marker regression using marker dosages
}
\examples{
data("SNP_dosages.4x","BLUEs.pheno")
Trait_1.smr <- singleMarkerRegression(dosage_matrix = SNP_dosages.4x,
Phenotype.df = BLUEs.pheno,genotype.ID = "Geno",trait.ID = "BLUE")
}
