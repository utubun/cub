% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Methods.R
\name{codons,DNAStringSet-method}
\alias{codons,DNAStringSet-method}
\title{Split DNAString / DNAString set into codons}
\usage{
\S4method{codons}{DNAStringSet}(object, ...)
}
\arguments{
\item{character}{DNAString or DNAStringSet representing the Gene or Genome}
}
\value{
character vector of codones, or list of characters for DNAStringSet
}
\description{
Split DNAString / DNAString set into codons
}
\examples{
dat <- PA14CDS[[1]]
codons(dat)
}
