% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polyqtlR_functions.R
\name{exploreQTL}
\alias{exploreQTL}
\title{Explore the possible segregation type of a QTL peak using Schwarz Information Criterion}
\usage{
exploreQTL(
  IBD_list,
  Phenotype.df,
  genotype.ID,
  trait.ID,
  linkage_group,
  LOD_data,
  cM = NULL,
  QTLconfig = NULL,
  plotBIC = TRUE,
  deltaBIC = 6,
  testAllele_Effects = TRUE,
  log = NULL
)
}
\arguments{
\item{IBD_list}{List of IBD probabilities}

\item{Phenotype.df}{A data.frame containing phenotypic values}

\item{genotype.ID}{The colname of \code{Phenotype.df} that contains the population identifiers (F1 names) (must be a colname of \code{Phenotype.df})}

\item{trait.ID}{The colname of \code{Phenotype.df} that contains the response variable to use in the model (must be a colname of \code{Phenotype.df})}

\item{linkage_group}{Numeric identifier of the linkage group being tested, based on the order of \code{IBD_list}.
Only a single linkage group is allowed.}

\item{LOD_data}{Output of \code{\link{QTLscan}} function}

\item{cM}{By default \code{NULL}, in which case the position of maximum LOD score is taken as the position of interest. 
Otherwise, the cM position to be explored.}

\item{QTLconfig}{Nested list of homologue configurations and modes of action of QTL to be explored and compared, the output of \code{\link{segMaker}}.
Note that a default List is available of all possible bi-allelic QTL if none is provided.
Each list element is itself a list with components
\describe{
 \item{homs : }{ a vector of length at least 1, describing the proposed homologues the functional allele Q is on}
 \item{mode : }{ Vector of same length as \code{homs} with codes "a" for additive and "d" for dominant.}
}}

\item{plotBIC}{Logical, with default \code{TRUE} - should the calculated BIC values be plotted?}

\item{deltaBIC}{Numeric, by default 6. Configurations within this distance of the minimum BIC are considered plausible.}

\item{testAllele_Effects}{Logical, with default \code{TRUE} - should the effects of the different alleles be tested
using the most likely QTL configuration?}

\item{log}{Character string specifying the log filename to which standard output should be written. If \code{NULL} log is send to stdout.}
}
\value{
List with the following items:
\describe{
\item{linkage_group}{Linkage group of the QTL peak being explored}
\item{cM}{CentiMorgan position of the locus being explored}
\item{BIC}{Vector of BIC values corresponding to elements of \code{QTLconfig} provided for testing}
\item{Allele.effects}{Summary of the means and standard errors of groups with (+)
and without(-) the specified allele combinations for the most likely QTLconfig
if \code{testAllele_Effects} = \code{TRUE} (\code{NULL} otherwise).}
\item{genotype.means}{A one-column matrix of mean phenotype values of offspring classes, with rownames
corresponding to the genotype class. If the probability of certain genotype classes is 0 (e.g. double reduction
classes where no double reduction occurred), then the genotype mean for that class will be \code{NA}}
}
}
\description{
Function to explore the possible segregation type at a QTL position using the Schwarz Information Criterion
}
\examples{
data("IBD_4x","BLUEs.pheno","qtl_LODs.4x")
exploreQTL(IBD_list = IBD_4x,
           Phenotype.df = BLUEs.pheno,
           genotype.ID = "Geno",
           trait.ID = "BLUE",
           linkage_group = 1,
           LOD_data = qtl_LODs.4x)
}
