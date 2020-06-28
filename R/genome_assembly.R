# forging bsgenome for P.aeruginosa PAO1 and PA14 strains
library(Biostrings)
library(BSgenome)

wd <- getwd()

setwd('./dat/gnm/asm/out')

# P.aeruginosa PA01

forgeBSgenomeDataPkg('../src/PAO1/seed.txt')

# installation from cl (tested on Ubuntu 18.04)
system('R CMD build BSgenome.Paeruginosa.NCBI.ASM676v1')
system('R CMD check BSgenome.Paeruginosa.NCBI.ASM676v1_1.0.0.tar.gz')
system('R CMD INSTALL BSgenome.Paeruginosa.NCBI.ASM676v1_1.0.0.tar.gz')

# P.aeruginosa PA14

forgeBSgenomeDataPkg('../src/PA14/seed.txt')

# installation from cl (tested on Ubuntu 18.04)
system('R CMD build BSgenome.Paeruginosa.NCBI.ASM1462v1')
system('R CMD check BSgenome.Paeruginosa.NCBI.ASM1462v1_1.0.0.tar.gz')
system('R CMD INSTALL BSgenome.Paeruginosa.NCBI.ASM1462v1_1.0.0.tar.gz')

# check installation

# library(BSgenome.Paeruginosa.NCBI.ASM676v1)
# library(BSgenome.Paeruginosa.NCBI.ASM1462v1)

# go back to working directory
setwd(wd)
