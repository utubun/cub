library(BSgenome.Ecoli.NCBI.K12.MG1655)
library(TxDb.Ecoli.K12.MG1655)

#devtools::load_all()

GNM <- BSgenome.Ecoli.NCBI.K12.MG1655
ANN <- TxDb.Ecoli.K12.MG1655

GRN <- GenomicFeatures::cds(ANN, columns = 'tx_name', use.names = TRUE)

k12 <- Biostrings::getSeq(GNM, GRN)

# nm <- biomaRt::select(ANN, names(k12), 'TXNAME', 'CDSNAME')
#
# names(k12) <- nm[match(nm[, 1], names(k12)), ][, 2]
#
# res <- cub::CubSet(k12)

# nm <- Sharp166$gene[Sharp166$expression == 'high']
#
# subs <- eck12[txNames$GENEID[txNames$TXNAME %in% nm]]
#
# hsubs <- cub::CubSet(subs)
#
# k12 <- cub::rscu(hsubs)

usethis::use_data(k12, overwrite = TRUE, compress = 'xz')
