# Load libraried
library(parallel)
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(BSgenome.Paeruginosa.NCBI.ASM676v1))
suppressPackageStartupMessages(library(BSgenome.Paeruginosa.NCBI.ASM1462v1))
library(xlsx)

# Load genomes

gnm_pao1 <- BSgenome.Paeruginosa.NCBI.ASM676v1
gnm_pa14 <- BSgenome.Paeruginosa.NCBI.ASM1462v1

# Load TxDb for both strains

txdb_pao1 <- loadDb('./dat/gnm/sql/TxDb.Paeruginosa.ASM676v1.sqlite')
txdb_pa14 <- loadDb('./dat/gnm/sql/TxDb.Paeruginosa.ASM1462v1.sqlite')

# extract CDS as granges

cds_pao1 <- cds(txdb_pao1, use.names = T)
cds_pa14 <- cds(txdb_pa14, use.names = T)

# # make views on genomes --------------------------------------------------------
#
# cds_pao1_views <- Views(gnm_pao1, cds_pao1)
# params <- new('BSParams', X = gnm_pao1, FUN = codons)
# cdns <- bsapply(params)

# extract CDS as a DNAStringSets

dna_pao1 <- getSeq(gnm_pao1, cds_pao1)
dna_pa14 <- getSeq(gnm_pa14, cds_pa14)

# # check that it is correctly extracted -----------------------------------------
#
# strand(cds_pao1[1])
# strand(cds_pa01[5036])
# (cds_01_coords <- c(start(cds_pa01[1]), end(cds_pa01[1])))
# (cds_5036_coords <- c(start(cds_pa01[5036]), end(cds_pa01[5036])))
# gnm_pa01$chr[cds_01_coords[1]:cds_01_coords[2]]
# dna_pa01[1]
# gnm_pa01$chr[cds_5036_coords[1]:cds_5036_coords[2]]
# dna_pa01 <- gnm_pa01$chr
# dna_pa01[5036]

triplets <- function(x) {
  #unlist(strsplit(toString(x), '(?<=.{3})', perl = TRUE), use.names = F)
  res <- tryCatch(
    {
      as.character(codons(x))
    },
    error = function(cond) {
      return(NA)
    }
  )

  return(res)
}

ncores <- detectCores()

cl     <- makeCluster(ncores)

cdn_pao1 <- parLapply(cl, dna_pao1, triplets)
cdn_pa14 <- parLapply(cl, dna_pa14, triplets)

stopCluster(cl)

dir.create('./dat/ftb', F, T)

readr::write_rds(cdn_pao1, './dat/ftb/ASM676v1.Rds')
readr::write_rds(cdn_pa14, './dat/ftb/ASM1462v1.Rds')

cdn_pao1 <- readr::read_rds('./dat/ftb/ASM676v1.Rds')
cdn_pa14 <- readr::read_rds('./dat/ftb/ASM1462v1.Rds')

dir.create('./dat/grn', F, T)

# count codons for each gene
cdnCount <- function(x) {
  ncores <- detectCores()
  cl     <- makeCluster(ncores)

  res <- parLapply(cl, x, function(datum) {
    tb <- table(datum)
    tb <- tibble::tibble(codon = names(tb), count = tb)
    return(tb)
    }
  )

  stopCluster(cl)
  rm(cl)

  return(res)
}

cdn_cnt_pao1 <- data.table::rbindlist(cdnCount(cdn_pao1), idcol = 'name') %>%
  group_by(codon) %>%
  summarise(count = sum(count)) %>%
  mutate(
    strain        = 'PAO1',
    usage_total   = count / sum(count) * 100,
    aa            = map_chr(codon, ~GENETIC_CODE[.x])
  ) %>%
  group_by(aa) %>%
  mutate(
    usage_aa = count / sum(count) * 100
  ) %>%
  ungroup() %>%
  dplyr::select(strain, codon, aa, usage_total, usage_aa) %>%
  arrange(aa, desc(usage_total), desc(usage_aa), codon)

cdn_cnt_pa14 <- data.table::rbindlist(cdnCount(cdn_pa14), idcol = 'name', fill = TRUE) %>%
  group_by(codon) %>%
  summarise(count = sum(count)) %>%
  mutate(
    strain        = 'PA14',
    usage_total   = count / sum(count) * 100,
    aa            = GENETIC_CODE[codon]
  ) %>%
  group_by(aa) %>%
  mutate(
    usage_aa = count / sum(count) * 100
  ) %>%
  ungroup() %>%
  dplyr::select(strain, codon, aa, usage_total, usage_aa) %>%
  arrange(aa, desc(usage_total), desc(usage_aa), codon)

# Combine data for total codon usage
bind_rows(cdn_cnt_pao1, cdn_cnt_pa14) %>%
  dplyr::select(strain, aa, codon, usage = usage_total) %>%
  pivot_wider(names_from = strain, values_from = usage) %>%
  rename(
    `Amino acid` = aa, Codon = codon
  ) %>%
  filter(!is.na(`Amino acid`)) %>%
  write.xlsx('./dat/cdn/cdn_count_total.xlsx')

# Combine data for per aa codon usage
bind_rows(cdn_cnt_pao1, cdn_cnt_pa14) %>%
  dplyr::select(strain, aa, codon, usage = usage_aa) %>%
  pivot_wider(names_from = strain, values_from = usage) %>%
  rename(
    `Amino acid` = aa, Codon = codon
  ) %>%
  filter(!is.na(`Amino acid`)) %>%
  write.xlsx('./dat/cdn/cdn_count_byaa.xlsx')
