# DECIFER TEST
library(tidyverse)
library(magrittr)
library(DECIPHER)
#library(RSQLite)

db = dbConnect(SQLite(), './dat/gnm/sql/genomes.db')

genomes <- c(
  PA01 = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/765/GCA_000006765.1_ASM676v1/GCA_000006765.1_ASM676v1_genomic.fna.gz',
  PA14 = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/014/625/GCA_000014625.1_ASM1462v1/GCA_000014625.1_ASM1462v1_genomic.fna.gz'
)

for(i in names(genomes)) {
  Seqs2DB(genomes[i], 'FASTA', db, i)
}

(synteny <- FindSynteny(db, 'Seqs', minScore = 50))

plot(synteny)
pairs(synteny)
head(synteny[[1, 2]])
head(synteny[[2, 1]])

(psae <- AlignSynteny(synteny, db, verbose = TRUE))

library(SynExtend)

GeneCalls <- vector(mode = 'list', length = ncol(synteny))

GeneCalls[[1L]] <- gffToDataFrame(
  GFF = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/765/GCA_000006765.1_ASM676v1/GCA_000006765.1_ASM676v1_genomic.gff.gz',
  Verbose = TRUE
)

GeneCalls[[2L]] <- gffToDataFrame(
  GFF = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/014/625/GCA_000014625.1_ASM1462v1/GCA_000014625.1_ASM1462v1_genomic.gff.gz',
  Verbose = TRUE
)

names(GeneCalls) <- names(genomes)

head(GeneCalls[[1]])
head(GeneCalls[[2]])

class(GeneCalls[[1]])

links <- NucleotideOverlap(synteny, GeneCalls, LimitIndex = FALSE, Verbose = TRUE)

linked_pairs <- PairSummaries(
  SyntenyLinks           = links,
  GeneCalls              = GeneCalls,
  DBPATH                 =  './dat/gnm/sql/genomes.db',
  PIDs                   = FALSE,
  IgnoreDefaultStringSet = FALSE,
  Verbose                = TRUE,
  Model                  = 'Global',
  Correction             = 'none'
)

# indices are in the form \\d_\\d_\\d e.g. 1_1_1 2_1_1 first pair of digits is constant
# for a given genome. The last is the index of gene, and the only important one
# it indicates gene index in GFF objects above. Use the indices to get gene info

linked_pairs %<>%
  rownames_to_column('idx') %>%
  filter(ModelSelect)  %>%
  select(-ModelSelect) %>%
  mutate(
    gn_01 = GeneCalls[[1]][as.numeric(gsub('.*_(\\d+)\\s.+', '\\1', idx)), ]$ID,
    gn_02 = GeneCalls[[2]][as.numeric(gsub('^.*_(\\d+)$', '\\1',     idx)), ]$ID,
    idx   = NULL
    )

# identical genes as a negative control

GeneCalls[[1]][GeneCalls[[1]]$ID %in% linked_pairs[linked_pairs$TotalCoverage == 1, ]$gn_01, ]
GeneCalls[[2]][GeneCalls[[2]]$ID %in% linked_pairs[linked_pairs$TotalCoverage == 1, ]$gn_02, ]
