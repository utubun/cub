ANN <- TxDb.Paeruginosa.PA14::TxDb.Paeruginosa.PA14
GNM <- BSgenome.Paeruginosa.NCBI.PA14::BSgenome.Paeruginosa.NCBI.PA14
GRN <- GenomicFeatures::cds(ANN, use.names = TRUE)
PA14CDS <- Biostrings::getSeq(GNM, GRN)
PA14CDS <- PA14CDS[names(PA14CDS) != '']
usethis::use_data(PA14CDS, overwrite = TRUE, compress = 'xz')
