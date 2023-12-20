PA14 <- cub::CubSet(cub::PA14CDS)
ANN <- TxDb.Paeruginosa.PA14::TxDb.Paeruginosa.PA14
lookup <- AnnotationDbi::mapIds(ANN, keys(ANN), 'CDSNAME', 'GENEID')
lookup <- setNames(names(lookup), lookup)
PA14@data <- setNames(PA14@data, names(cub::PA14CDS@ranges))
PA14@data <- setNames(PA14@data, lookup[names(PA14@data)])
usethis::use_data(PA14, overwrite = TRUE, compress = 'xz')
