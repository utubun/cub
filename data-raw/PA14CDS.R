ANN <- AnnotationDbi::loadDb("./data-raw/txdb/TxDb.Paeruginosa.ASM1462v1.sqlite")
GNM <- BSgenome.Paeruginosa.NCBI.ASM1462v1::BSgenome.Paeruginosa.NCBI.ASM1462v1

GRN <- GenomicFeatures::cds(ANN, use.names = TRUE)

PA14CDS <- Biostrings::getSeq(GNM, GRN)

usethis::use_data(PA14CDS, overwrite = TRUE, compress = 'xz')
