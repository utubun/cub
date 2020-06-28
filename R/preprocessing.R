# Load libraried
library(parallel)
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(BSgenome.Paeruginosa.NCBI.ASM676v1))
suppressPackageStartupMessages(library(BSgenome.Paeruginosa.NCBI.ASM1462v1))

# Load genomes

#gnm_pao1 <- BSgenome.Paeruginosa.NCBI.ASM676v1
gnm_pa14 <- BSgenome.Paeruginosa.NCBI.ASM1462v1

# Load TxDb for both strains

#txdb_pao1 <- loadDb('./dat/gnm/sql/TxDb.Paeruginosa.ASM676v1.sqlite')
txdb_pa14 <- loadDb('./dat/gnm/sql/TxDb.Paeruginosa.ASM1462v1.sqlite')

# extract CDS as granges

#cds_pao1 <- cds(txdb_pao1, use.names = T)
cds_pa14 <- cds(txdb_pa14, use.names = T)

# # make views on genomes --------------------------------------------------------
#
# cds_pao1_views <- Views(gnm_pao1, cds_pao1)
# params <- new('BSParams', X = gnm_pao1, FUN = codons)
# cdns <- bsapply(params)

# extract CDS as a DNAStringSets

#dna_pao1 <- getSeq(gnm_pao1, cds_pao1)
dna_pa14 <- getSeq(gnm_pa14, cds_pa14)

# # check that it is correctly extracted -----------------------------------------
#
# strand(cds_pa01[1])
# strand(cds_pa01[5036])
# (cds_01_coords <- c(start(cds_pa01[1]), end(cds_pa01[1])))
# (cds_5036_coords <- c(start(cds_pa01[5036]), end(cds_pa01[5036])))
# gnm_pa01$chr[cds_01_coords[1]:cds_01_coords[2]]
# dna_pa01[1]
# gnm_pa01$chr[cds_5036_coords[1]:cds_5036_coords[2]]
# dna_pa01 <- gnm_pa01$chr
# dna_pa01[5036]

triplets <- function(x) {
  unlist(strsplit(toString(x), '(?<=.{3})', perl = TRUE), use.names = F)
}

ncores <- detectCores()

cl     <- makeCluster(ncores)

#cdn_pao1 <- parLapply(cl, dna_pao1, triplets)
cdn_pa14 <- parLapply(cl, dna_pa14, triplets)

stopCluster(cl)

dir.create('./dat/ftb', F, T)

#readr::write_rds(cdn_pao1, './dat/ftb/ASM676v1.Rds')
readr::write_rds(cdn_pa14, './dat/ftb/ASM1462v1.Rds')

#cdn_pao1 <- readr::read_rds('./dat/ftb/ASM676v1.Rds')
cdn_pa14 <- readr::read_rds('./dat/ftb/ASM1462v1.Rds')

dir.create('./dat/grn', F, T)

# calculate relative tRNA concentration for each codon -------------------------

freq <- readr::read_rds('./dat/ftb/cdn_freq_lookup.Rds')

cdn_freq_list_pa14  <- NumericList(sapply(cdn_pa14, function(codon) unlist(freq[codon], use.names = F)))

readr::write_rds(cdn_freq_list_pa14, './dat/ftb/cdn_freq_list_pa14.Rds')

cdn_usage_list_pa14 <- lapply(
  cdn_pa14,
  function(codons) {
    cdn_tab <- table(codons)
    setNames(as.vector(cdn_tab) / sum(cdn_tab), names(cdn_tab))
  }
)

max_lenght = max(sapply(cdn_freq_list_pa14, length))

cdn_usage_pa14_mat <- data.table::rbindlist(lapply(cdn_usage_list_pa14, as.list), fill = TRUE)
cdn_usage_pa14_mat <- as.matrix(cdn_usage_pa14_mat)
cdn_usage_pa14_mat[is.na(cdn_usage_pa14_mat)] <- 0
rownames(cdn_usage_pa14_mat) <- names(cdn_freq_list_pa14)

library(kohonen)
x <- as.matrix(pa14_spline@freq)
cdn_som <- som(x)
plot(cdn_som)
plot(cdn_som, type = "dist.neighbours", palette.name = terrain.colors)
pca <- prcomp(t(cdn_usage_pa14_mat), center = FALSE, scale = FALSE)

library(factoextra)
library(ggthemes)

fviz_screeplot(
  pca,
  choice    = 'variance',
  barfill   = 'gray50',
  barcolor  = 'gray50',
  linecolor = 'gray10',
  ggtheme   = theme_few(),
  main      = '',
  caption   = 'Gen III, PCA: Variances of dimensions'
) +
  theme(
    text = element_text(family = 'Alegreya', color = 'gray10')
  )

main_components <- as.matrix(pca$rotation[, 1])

fviz_nbclust(
  main_components,
  FUNcluster = kmeans,
  method     = 'wss',
  linecolor = 'red'
) +
  labs(title = element_blank(), caption = 'Gen III, PCA: Optimal number of clusters') +
  theme_few() +
  theme(
    text = element_text(family = 'Alegreya', color = 'gray10')
  )

set.seed(11235813)
kmean_clus <- kmeans(main_components, 2, nstart = 50, iter.max = 1000)

fviz_cluster(
  kmean_clus,
  data = main_components,
  legend = FALSE,
  palette = RColorBrewer::brewer.pal(3, 'Set2')
) +
  labs(
    title   = element_blank(),
    caption = 'Gen III, PCA: Three clusters of compounds'
  ) +
  theme_few() +
  theme(
    text         = element_text(family = 'Alegreya'),
    legend.title = element_blank()
  )


# assign relative tRNA frequencies to granges objects for both strains ---------

cds_pa14$freq <- cdn_freq_list_pa14
