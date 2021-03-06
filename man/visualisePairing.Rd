% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polyqtlR_functions.R
\name{visualisePairing}
\alias{visualisePairing}
\title{Visualise pairing of parental homologues}
\usage{
visualisePairing(
  meiosis_report.ls,
  pos.col = "red",
  neg.col = "blue",
  parent,
  max.lwd = 20,
  datawidemax,
  add.label = TRUE,
  return.data = FALSE,
  ...
)
}
\arguments{
\item{meiosis_report.ls}{List output of function \code{\link{meiosis_report}}}

\item{pos.col}{Colour corresponding to excess of pairing associations predicted (positive deviations), by default red}

\item{neg.col}{Colour corresponding to lack of pairing associations predicted (negative deviations), by default blue}

\item{parent}{The parent, either "P1" (mother) or "P2 (father)}

\item{max.lwd}{Maximum line width, by default 20}

\item{datawidemax}{This argument is currently a work-around to allow multiple plots to have the same scale (line thicknesses consistent).
No default is provided. To estimate this value, simply set argument \code{return.data = TRUE}, and record the
maximum absolute value over columns 'count', which are the deviations from random expectations. This should be done
over multiple function calls if e.g. comparing both P1 and P2 values. When a global maximum (absolute) deviation is known,
re-run the function with this value for \code{datawidemax}. The line width specified by \code{max.lwd} will then be
 used for this, and all other line widths re-scaled accordingly.}

\item{add.label}{Should a label be applied, giving the maximum deviation in the plot? By default \code{TRUE}}

\item{return.data}{Should plot data be returned? By default \code{FALSE}}

\item{\dots}{Optional arguments passed to \code{\link[igraph]{plot.igraph}}}
}
\value{
If \code{return.data = TRUE}, the values for pairwise deviations from the expected numbers are 
returned, useful for determining the value \code{datawidemax} to provide consistent scaling across multiple plots
}
\description{
Function to visualise the pairing of parental homologues across the population using graph, with nodes 
to denote parental homologues and edges to denote deviations from expected proportions under a polysomic model of inheritance
}
\examples{
data("mr.ls")
visualisePairing(meiosis_report.ls = mr.ls,
                 parent = "P1",
                 datawidemax = 3)
}
