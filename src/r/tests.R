pao1_txdb <- loadDb('./data/txdb/sql/TxDb.Paeruginosa.pao1.sqlite')
seqlevels(pao1_txdb) <- seqnames(pao1)

# split gene to codones, may be done by SlidingWindow with width = 3 on granges
# but can't be dnastrings, only granges

#cds_range <- cdsBy(pao1_txdb, by = 'tx', use.names = TRUE)
cds_rg <- cds(pao1_txdb, use.names = T)
cds_sq <- getSeq(pao1, cds_rg)
rna <- RNAStringSet(reverseComplement(tx_sq), use.names = T)
names(rna) <- names(cds_sq)
# see that extractTranscriptSeqs return transcript from negative strand
#BSgenome.Paeruginosa.EMBL.PAO1$AE004091.2[6259671:6260054] # genomic
#sds_seqs$PA5566

library(parallel)
ncores <- detectCores()
cl     <- makeCluster(ncores)
cdns <- parLapply(cl, rna, function(x) as.character(codons(x)))
stopCluster(cl)

names(cdns) <- names(rna)

codon_table <- unique(unlist(cdns))
set.seed(11235813)
codon_table_list <-abs(rnorm(length(codon_table)))
codon_table_list <- codon_table_list / sum(codon_table_list)
names(codon_table_list) <- codon_table
codon_table_list
cdns_freq_tab <- sapply(cdns, function(i) codon_table_list[i])

cds_rg$freqs <- cdns_freq_tab
# table E.coli from here https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/codon-usage-tables/master/codon_usage_data/tables/e_coli_316407.csv
cdn_table <- readr::read_csv('./data/codonusage/e_coli_316407.csv')
cdn_table_list <- as.list(cdn_table$relative_frequency)
names(cdn_table_list) <- cdn_table$codon
cdns_freq_tab <- sapply(cdns, function(i) codon_table_list[i])
cds_rg$freqs <- cdns_freq_tab

x <- sapply(cdns_freq_tab, function(x) cumsum(1/x))
