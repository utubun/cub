## code to prepare `DATASET` dataset goes here
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(BSgenome.Paeruginosa.NCBI.ASM676v1))
suppressPackageStartupMessages(library(BSgenome.Paeruginosa.NCBI.ASM1462v1))

dat <- loadDb("./data/TxDb.Paeruginosa.ASM676v1.sqlite")

usethis::use_data(DATASET, overwrite = TRUE)
