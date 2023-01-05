\name{BSgenome.Paeruginosa.NCBI.ASM1462v1}
\docType{package}

\alias{BSgenome.Paeruginosa.NCBI.ASM1462v1-package}
\alias{BSgenome.Paeruginosa.NCBI.ASM1462v1}
\alias{Paeruginosa}

\title{Full genome of Pseudomonas aeruginosa UCBPP-PA14}

\description{
  Complete Genome
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

\author{Dmytro Strunin}

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
BSgenome.Paeruginosa.NCBI.ASM1462v1
genome <- BSgenome.Paeruginosa.NCBI.ASM1462v1
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
