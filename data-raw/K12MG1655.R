## code to prepare `K12MG1655` dataset goes here
ANN <- GenomicFeatures::makeTxDbFromGFF(
  './data-raw/txdb/GCF_000005845.2_ASM584v2_genomic.gff.gz',
  'gff',
  'NCBI',
  'Escherichia coli str. K-12 substr. MG1655',
  511145
)

GNM <- BSgenome.Ecoli.NCBI.K12.MG1655::BSgenome.Ecoli.NCBI.K12.MG1655

GRN <- GenomicFeatures::genes(ANN)

K12MG1655GN <- getSeq(GNM, GRN)

# map between b-ids and names
txNames <- select(ANN, names(K12MG1655GN), 'TXNAME', 'GENEID')

# very highly expressed
vhgenes <- c('rpsU', 'rpsJ', 'rpsT', 'rpsL', 'rpsA', 'rpsB', 'rpsO', 'rpsG', 'rpmB',
            'rpmG', 'rpmH', 'rplK', 'rplJ', 'rplA', 'rplL', 'rplQ', 'rplC', 'lpp',
            'ompA', 'ompC', 'ompF', 'tufA', 'tufB', 'tsf', 'fusA', 'recA', 'dnaK')
mean(vhgenes %in% txNames[, 2])

# highly expressed
hgenes <- c(
  'rpoA', 'rpoB', 'rpoC', 'rpoD', 'uncA', 'uncD', 'uncE', 'alaS', 'metG', 'glnS',
  'glyS', 'gly2', 'thrS', 'trpS', 'tyrS'
)

hgenes %in% txNames[, 2] # uncA-E renamed

# K12 <- cub::CubSet(K12MG1655GN)

subs <- K12MG1655GN[txNames$GENEID[txNames$TXNAME %in% vhgenes]]
hsubs <- cub::CubSet(subs)

usethis::use_data(K12MG1655GN, overwrite = TRUE)
