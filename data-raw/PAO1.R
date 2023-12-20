PAO1 <- cub::CubSet(cub::PAO1CDS)
ANN <- TxDb.Paeruginosa.PAO1::TxDb.Paeruginosa.PAO1
lookup <- AnnotationDbi::mapIds(ANN, keys(ANN), 'CDSNAME', 'GENEID')
lookup <- setNames(names(lookup), lookup)
PAO1@data <- setNames(PAO1@data, names(cub::PAO1CDS@ranges))
PAO1@data <- setNames(PAO1@data, lookup[names(PAO1@data)])
usethis::use_data(PAO1, overwrite = TRUE, compress = 'xz')
