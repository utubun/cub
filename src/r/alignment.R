# DECIFER TEST
library(tidyverse)
library(magrittr)
library(DECIPHER)
library(SynExtend)
#library(RSQLite)

db <- dbConnect(SQLite(), './dat/gnm/sql/genomes.db')

genomes <- c(
  'AE004091.2' = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/765/GCA_000006765.1_ASM676v1/GCA_000006765.1_ASM676v1_genomic.fna.gz',
  'CP000438.1' = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/014/625/GCA_000014625.1_ASM1462v1/GCA_000014625.1_ASM1462v1_genomic.fna.gz'
)

dbRemoveTable(db, 'Seqs')
dbRemoveTable(db, '_Seqs')

for(i in names(genomes)) {
  Seqs2DB(genomes[i], 'FASTA', db, i)
}

dbListTables(db)
dbGetTable(db, 'Seqs')

rna <- SearchDB(db, nameBy = 'identifier', type = "RNAStringSet")

(synteny <- FindSynteny(db, 'Seqs', minScore = 50))
write_rds(synteny, './dat/gnm/syn/synteny.Rds')

dir.create('./dat/img')

png('./dat/img/synteny.png')
plot(synteny, main = 'Synteny')
dev.off()

png('./dat/img/pairs.png')
pairs(synteny)
dev.off()

head(synteny[[1, 2]])
head(synteny[[2, 1]])

(psae <- AlignSynteny(synteny, db, verbose = TRUE))

readr::write_rds(psae, './dat/gnm/syn/syn_aligned.Rds')
# Synthextend part

GeneCalls <- vector(mode = 'list', length = ncol(synteny))

GeneCalls[[1L]] <- gffToDataFrame(
  GFF = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/765/GCA_000006765.1_ASM676v1/GCA_000006765.1_ASM676v1_genomic.gff.gz',
  AdditionalTypes = 'CDS',
  Verbose         = TRUE
)

GeneCalls[[2L]] <- gffToDataFrame(
  GFF = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/014/625/GCA_000014625.1_ASM1462v1/GCA_000014625.1_ASM1462v1_genomic.gff.gz',
  AdditionalTypes = 'CDS',
  Verbose         = TRUE
)

names(GeneCalls) <- names(genomes)

GeneCalls <- lapply(
  GeneCalls,
  function(i) {
    i <- i[i$Type == 'CDS', ]
    i$ID <- gsub('^cds-', '', i$ID)
    i
  }
)

readr::write_rds(GeneCalls, './dat/gnm/syn/gene_calls.Rds')

head(GeneCalls[[1]])
head(GeneCalls[[2]])

class(GeneCalls[[1]])

links <- NucleotideOverlap(synteny, GeneCalls, LimitIndex = FALSE, Verbose = TRUE)

readr::write_rds(links, './dat/gnm/syn/gene_links.Rds')

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

readr::write_rds(linked_pairs, './dat/gnm/syn/linked_pairs.Rds')

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

readr::write_rds(linked_pairs, './dat/gnm/syn/linked_pairs.Rds')

# identical genes as a negative control

GeneCalls[[1]][GeneCalls[[1]]$ID %in% linked_pairs[linked_pairs$TotalCoverage == 1, ]$gn_01, ]
GeneCalls[[2]][GeneCalls[[2]]$ID %in% linked_pairs[linked_pairs$TotalCoverage == 1, ]$gn_02, ]

# Disconnect database!

dbDisconnect(db)
