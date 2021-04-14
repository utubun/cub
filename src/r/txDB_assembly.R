library('GenomicFeatures')
library(BSgenome.Paeruginosa.NCBI.ASM676v1)
library(BSgenome.Paeruginosa.NCBI.ASM1462v1)

pao1_gff <- 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/765/GCA_000006765.1_ASM676v1/GCA_000006765.1_ASM676v1_genomic.gff.gz'
pa14_gff <- 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/014/625/GCA_000014625.1_ASM1462v1/GCA_000014625.1_ASM1462v1_genomic.gff.gz'

# P.aeruginosa strain P01

pa01_txdb <- makeTxDbFromGFF(
  file       = pao1_gff,
  organism   = 'Pseudomonas aeruginosa',
  dataSource = 'NCBI, Pseudomonas aeruginosa PAO1, assembly ASM676v1',
  circ_seqs  = 'AE004091.2'
)

dir.create('./dat/gnm/sql', F, T)

saveDb(pa01_txdb, file = './dat/gnm/sql/TxDb.Paeruginosa.ASM676v1.sqlite')

rm(pa01_txdb)

# P.aeruginosa strain PA14

pa14_txdb <- makeTxDbFromGFF(
  file       = pa14_gff,
  organism   = 'Pseudomonas aeruginosa',
  dataSource = 'NCBI, Pseudomonas aeruginosa PA14, assembly ASM1462v1',
  circ_seqs  = 'CP000438.1'
)

saveDb(pa14_txdb, file = './dat/gnm/sql/TxDb.Paeruginosa.ASM1462v1.sqlite')

rm(pa14_txdb)
