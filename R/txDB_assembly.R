library('GenomicFeatures')
library(BSgenome.Paeruginosa.NCBI.ASM676v1)
library(BSgenome.Paeruginosa.NCBI.ASM1462v1)
# set location to the project directory!

setwd('~/Documents/prj/bac/ribotracker')
# P.aeruginosa strain P01

pa01_txdb <- makeTxDbFromGFF(
  file       = './dat/gnm/asm/src/PA01/chr.gff',
  dataSource = 'NCBI, Pseudomonas aeruginosa PAO1, assembly ASM676v1',
  organism   = 'Pseudomonas aeruginosa'
)

seqlevels(pa01_txdb)  <- 'chr'

dir.create('./dat/gnm/sql', F, T)

saveDb(pa01_txdb, file = './dat/gnm/sql/TxDb.Paeruginosa.ASM676v1.sqlite')

rm(pa01_txdb)

# P.aeruginosa strain PA14

pa14_txdb <- makeTxDbFromGFF(
  file = './dat/gnm/asm/src/PA14/chr.gff',
  dataSource = 'NCBI, Pseudomonas aeruginosa PA14, assembly ASM1462v1',
  organism = 'Pseudomonas aeruginosa'
)

# Warning message:
#   In makeTxDbFromGRanges(gr, metadata = metadata) :
#   The following transcripts were dropped because their exon ranks could not be inferred
# (either because the exons are not on the same chromosome/strand or because they are
#   not separated by introns): gene-PA14_RS01280, gene-PA14_RS01345, gene-PA14_RS21005,
# gene-PA14_RS22430

seqlevels(pa14_txdb)  <- 'chr'

saveDb(pa14_txdb, file = './dat/gnm/sql/TxDb.Paeruginosa.ASM1462v1.sqlite')

rm(pa14_txdb)
