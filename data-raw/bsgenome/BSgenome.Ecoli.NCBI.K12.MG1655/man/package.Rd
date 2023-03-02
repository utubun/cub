\name{BSgenome.Ecoli.NCBI.K12.MG1655}
\docType{package}

\alias{BSgenome.Ecoli.NCBI.K12.MG1655-package}
\alias{BSgenome.Ecoli.NCBI.K12.MG1655}
\alias{Ecoli}

\title{Full genome of Escherichia coli str. K-12 substr. MG1655}

\description{
  Full genome sequence Escherichia coli str. K-12 substr. MG1655 provided by NCBI
}

\details{
  
}

\note{
  This BSgenome data package was made from the following source data files:
  \preformatted{
-- information not available --
  }

  See \code{?\link[BSgenome]{BSgenomeForge}} and the BSgenomeForge
  vignette (\code{vignette("BSgenomeForge")}) in the \pkg{BSgenome}
  software package for how to make a BSgenome data package.
}

\author{The Bioconductor Dev Team}

\seealso{
  \itemize{
    \item \link[BSgenome]{BSgenome} objects and the
          \code{\link[BSgenome]{available.genomes}} function
          in the \pkg{BSgenome} software package.
    \item \link[Biostrings]{DNAString} objects in the \pkg{Biostrings}
          package.
    \item The BSgenomeForge vignette (\code{vignette("BSgenomeForge")})
          in the \pkg{BSgenome} software package for how to make a BSgenome
          data package.
  }
}

\examples{
BSgenome.Ecoli.NCBI.K12.MG1655
genome <- BSgenome.Ecoli.NCBI.K12.MG1655
head(seqlengths(genome))


## ---------------------------------------------------------------------
## Genome-wide motif searching
## ---------------------------------------------------------------------
## See the GenomeSearching vignette in the BSgenome software
## package for some examples of genome-wide motif searching using
## Biostrings and the BSgenome data packages:
if (interactive())
    vignette("GenomeSearching", package="BSgenome")
}

\keyword{package}
\keyword{data}
