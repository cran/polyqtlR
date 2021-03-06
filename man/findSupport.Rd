% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polyqtlR_functions.R
\name{findSupport}
\alias{findSupport}
\title{Function to find a LOD - x support interval around a QTL position}
\usage{
findSupport(LOD_data, linkage_group, LOD_support = 2)
}
\arguments{
\item{LOD_data}{Output of \code{\link{QTLscan}} function.}

\item{linkage_group}{Numeric identifier of the linkage group being tested, based on the order of \code{IBD_list}.
Only a single linkage group is allowed.}

\item{LOD_support}{The level of support around a QTL peak, by default 2 (giving a LOD - 2 support interval, the 
range of positions with a LOD score within 2 LOD units of the maximum LOD on that linkage group).}
}
\description{
Given QTL output, this function returns the LOD - x support for a specified linkage group, taking the 
maximum LOD position as the desired QTL peak.
}
\examples{
data("qtl_LODs.4x")
findSupport(LOD_data=qtl_LODs.4x,linkage_group=1)
}
