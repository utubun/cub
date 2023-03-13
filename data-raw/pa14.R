library(BSgenome.Paeruginosa.NCBI.PA14)
library(TxDb.Paeruginosa.PA14)

GNM <- BSgenome.Paeruginosa.NCBI.PA14
ANN <- TxDb.Paeruginosa.PA14

GRN <- GenomicFeatures::cds(ANN, columns = 'tx_name', use.names = TRUE)

pa14 <- Biostrings::getSeq(GNM, GRN)

#nm <- biomaRt::select(ANN, names(pa14), 'TXNAME', 'CDSNAME')

#names(pa14) <- nm[match(nm[, 1], names(k12)), ][, 2]

#res <- cub::CubSet(k12)

# nm <- Sharp166$gene[Sharp166$expression == 'high']
#
# subs <- eck12[txNames$GENEID[txNames$TXNAME %in% nm]]
#
# hsubs <- cub::CubSet(subs)
#
# k12 <- cub::rscu(hsubs)

usethis::use_data(pa14, overwrite = TRUE, compress = 'xz')
