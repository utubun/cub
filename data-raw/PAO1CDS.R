ANN <- TxDb.Paeruginosa.PAO1::TxDb.Paeruginosa.PAO1
GNM <- BSgenome.Paeruginosa.NCBI.PAO1::BSgenome.Paeruginosa.NCBI.PAO1
GRN <- GenomicFeatures::cds(ANN, use.names = TRUE)

PAO1CDS <- Biostrings::getSeq(GNM, GRN)

usethis::use_data(PAO1CDS, overwrite = TRUE, compress = 'xz')
