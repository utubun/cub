ANN <- AnnotationDbi::loadDb("./data-raw/txdb/TxDb.Paeruginosa.ASM676v1.sqlite")
GNM <- BSgenome.Paeruginosa.NCBI.ASM676v1::BSgenome.Paeruginosa.NCBI.ASM676v1

GRN <- GenomicFeatures::cds(ANN, use.names = TRUE)

PAO1CDS <- Biostrings::getSeq(GNM, GRN)

usethis::use_data(PAO1CDS, overwrite = TRUE, compress = 'xz')
