pkgname <- "BSgenome.Paeruginosa.NCBI.ASM676v1"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('BSgenome.Paeruginosa.NCBI.ASM676v1')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("package")
### * package

flush(stderr()); flush(stdout())

### Name: BSgenome.Paeruginosa.NCBI.ASM676v1
### Title: Full genome of Pseudomonas aeruginosa PAO1 (g-proteobacteria)
### Aliases: BSgenome.Paeruginosa.NCBI.ASM676v1-package
###   BSgenome.Paeruginosa.NCBI.ASM676v1 Paeruginosa
### Keywords: package data

### ** Examples

BSgenome.Paeruginosa.NCBI.ASM676v1
genome <- BSgenome.Paeruginosa.NCBI.ASM676v1
head(seqlengths(genome))


## ---------------------------------------------------------------------
## Genome-wide motif searching
## ---------------------------------------------------------------------
## See the GenomeSearching vignette in the BSgenome software
## package for some examples of genome-wide motif searching using
## Biostrings and the BSgenome data packages:
if (interactive())
    vignette("GenomeSearching", package="BSgenome")



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
