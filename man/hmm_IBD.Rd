% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/accessory_functions.R
\name{hmm_IBD}
\alias{hmm_IBD}
\title{Generate IBD probabilities from marker genotypes and a phased linkage map using HMM}
\usage{
hmm_IBD(
  input_type = "discrete",
  genotypes,
  phased_maplist,
  remove_markers = NULL,
  ploidy,
  ploidy2 = NULL,
  parent1 = "P1",
  parent2 = "P2",
  individuals = "all",
  log = NULL,
  map_function = "haldane",
  bivalent_decoding = TRUE,
  error = 0.01,
  full_multivalent_hexa = FALSE,
  verbose = FALSE,
  ncores = 1
)
}
\arguments{
\item{input_type}{Can be either one of 'discrete' or 'probabilistic'. For the former (default), \code{dosage_matrix} must be supplied,
while for the latter \code{probgeno_df} must be supplied}

\item{genotypes}{Marker genotypes, either a 2d matrix of integer marker scores or a data.frame of dosage probabilities. 
Details are as follows:
\itemize{
\item{discrete : }{
 If \code{input_type} is 'discrete', \code{genotypes} is a matrix of marker dosage scores with markers in rows and individuals in columns.
 Both (marker) rownames and (individual or sample) colnames are needed.
}
\item{probabilistic : }{
If \code{input_type} is 'probabilistic', \code{genotypes} is a data frame as read from the scores file produced by function \code{saveMarkerModels} of R package 
\code{fitPoly}, or alternatively, a data frame containing at least the following columns:
\itemize{
\item{SampleName : }{
Name of the sample (individual)
}
\item{MarkerName : }{
Name of the marker
}
\item{P0 : }{
Probabilities of dosage score '0'
}
\item{P1, P2... etc. : }{
Probabilities of dosage score '1' etc. (up to max offspring dosage, e.g. P4 for tetraploid population)
}
}
}
}}

\item{phased_maplist}{A list of phased linkage maps, the output of \code{polymapR::create_phased_maplist}}

\item{remove_markers}{Optional vector of marker names to remove from the maps. Default is \code{NULL}.}

\item{ploidy}{Integer. Ploidy of the organism.}

\item{ploidy2}{Optional integer, by default \code{NULL}. Ploidy of parent 2, if different from parent 1.}

\item{parent1}{Identifier of parent 1, by default assumed to be \code{"P1"}}

\item{parent2}{Identifier of parent 2, by default assumed to be \code{"P2"}}

\item{individuals}{By default "all" offspring are included, but otherwise a subset can be selected, using a vector of offspring indexing numbers (1,2, etc.)
according to their order in \code{dosage_matrix}}

\item{log}{Character string specifying the log filename to which standard output should be written. If \code{NULL} log is send to stdout.}

\item{map_function}{Mapping function to use when converting map distances to recombination frequencies.
Currently only \code{"haldane"} or \code{"kosambi"} are allowed.}

\item{bivalent_decoding}{Option to consider only bivalent pairing during formation of gametes (ignored for diploid populations, as only bivalents possible there), by default \code{TRUE}}

\item{error}{The (prior) probability of errors in the offspring dosages, usually assumed to be small but non-zero}

\item{full_multivalent_hexa}{Option to allow multivalent pairing in both parents at the hexaploid level, by default \code{FALSE}. Note that if \code{TRUE},
a very large available RAM may be required (>= 32Gb) to process the data.}

\item{verbose}{Logical, by default \code{TRUE}. Should progress messages be written?}

\item{ncores}{How many CPU cores should be used in the evaluation? By default 1 core is used.}
}
\value{
A list of IBD probabilities, organised by linkage group (as given in the input \code{phased_maplist}). Each
list item is itself a list containing the following:
\itemize{
\item{IBDtype: The type of IBD; for this function only "genotypeIBD" are calculated.}
\item{IBDarray: A 3d array of IBD probabilities, with dimensions marker, genotype-class and F1 individual.}
\item{map: A 3-column data-frame specifying chromosome, marker and position (in cM)}
\item{parental_phase: Phasing of the markers in the parents, as given in the input \code{phased_maplist}}
\item{biv_dec: Logical, whether bivalent decoding was used in the estimation of the F1 IBD probabilities.}
\item{gap: The size of the gap (in cM) used when interpolating the IBD probabilities, if performed.}
\item{genocodes: Ordered list of genotype codes used to represent different genotype classes.}
\item{pairing: log likelihoods of each of the different pairing scenarios considered (can be used e.g. for post-mapping check of preferential pairing)}
\item{ploidy: ploidy of parent 1}
\item{ploidy2: ploidy of parent 2}
\item{method : The method used, here "hmm" (Hidden Markov Model)}
\item{error : The error prior used}
}
}
\description{
\code{hmm_IBD} is a function for creating identity-by-descent (IBD) probabilities using hidden Markov models,
from marker genotypes (either discrete marker dosages (ie scores 0, 1, ..., ploidy representing the number of copies of the marker allele),
or the probabilities of these dosages) and a phased linkage map. Unlike the original TetraOrigin software, it does not re-estimate parental linkage phase,
and has been generalised for use in diploid, triploid, tetraploid and hexaploid populations.
}
\examples{
data("phased_maplist.4x", "SNP_dosages.4x")
hmm_IBD(phased_maplist=phased_maplist.4x,genotypes=SNP_dosages.4x,ploidy=4)
}
\references{
\itemize{
\item{Durbin R, Eddy S, Krogh A, Mitchison G (1998) Biological sequence analysis: Probabilistic models of proteins and nucleic acids. Cambridge: Cambridge University Press.}
\item{Hackett et al. (2013) Linkage analysis and QTL mapping using SNP dosage data in a tetraploid potato mapping population. PLoS One 8(5): e63939}
\item{Zheng et al. (2016) Probabilistic multilocus haplotype reconstruction in outcrossing tetraploids. Genetics 203: 119-131}
}
}
