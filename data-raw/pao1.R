library(BSgenome.Paeruginosa.NCBI.PAO1)
library(TxDb.Paeruginosa.PAO1)

GNM <- BSgenome.Paeruginosa.NCBI.PAO1
ANN <- TxDb.Paeruginosa.PAO1

GRN <- GenomicFeatures::cds(ANN, columns = 'tx_name', use.names = TRUE)

pao1 <- Biostrings::getSeq(GNM, GRN)

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

usethis::use_data(pao1, overwrite = TRUE, compress = 'xz')
