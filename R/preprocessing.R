# Load libraried
library(parallel)
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(GenomicFeatures))
library(BSgenome.Paeruginosa.NCBI.ASM676v1)
library(BSgenome.Paeruginosa.NCBI.ASM1462v1)

# Load genomes

gnm_pa01 <- BSgenome.Paeruginosa.NCBI.ASM676v1
gnm_pa14 <- BSgenome.Paeruginosa.NCBI.ASM1462v1

# Load TxDb for both strains

txdb_pa01 <- loadDb('./dat/gnm/sql/TxDb.Paeruginosa.ASM676v1.sqlite')
txdb_pa14 <- loadDb('./dat/gnm/sql/TxDb.Paeruginosa.ASM1462v1.sqlite')

# Make seqlevels for each strain equal to genome seqnames

seqlevels(txdb_pa01) <- seqnames(gnm_pa01)
seqlevels(txdb_pa14) <- seqnames(gnm_pa14)

# extract CDS as granges

cds_pa01 <- cds(txdb_pa01, use.names = T)
cds_pa14 <- cds(txdb_pa14, use.names = T)

# make views on genomes --------------------------------------------------------

cds_pa01_view <- Views(gnm_pa01, cds_pa01)
params <- new('BSParams', X = gnm_pa01, FUN = codons)
cdns <- bsapply(params)

# extract CDS as a DNAStringSets

dna_pa01 <- getSeq(gnm_pa01, cds_pa01)
dna_pa14 <- getSeq(gnm_pa14, cds_pa14)

# check that it is correctly extracted -----------------------------------------

strand(cds_pa01[1])
strand(cds_pa01[5036])
(cds_01_coords <- c(start(cds_pa01[1]), end(cds_pa01[1])))
(cds_5036_coords <- c(start(cds_pa01[5036]), end(cds_pa01[5036])))
gnm_pa01$chr[cds_01_coords[1]:cds_01_coords[2]]
dna_pa01[1]
gnm_pa01$chr[cds_5036_coords[1]:cds_5036_coords[2]]
dna_pa01 <- gnm_pa01$chr
dna_pa01[5036]

# triplets <- function(x) {
#   unlist(strsplit(toString(x), '(?<=.{3})', perl = TRUE), use.names = F)
# }
#
# ncores <- detectCores()
#
# cl     <- makeCluster(ncores)
#
# cdn_pa01 <- parLapply(cl, dna_pa01, triplets)
# cdn_pa14 <- parLapply(cl,  dna_pa14, triplets)
#
# stopCluster(cl)
#
# dir.create('./dat/ftb', F, T)
#
# readr::write_rds(cdn_pa01, './dat/ftb/ASM676v1.Rds')
# readr::write_rds(cdn_pa14, './dat/ftb/ASM1462v1.Rds')

cdn_pa01 <- readr::read_rds('./dat/ftb/ASM676v1.Rds')
cdn_pa14 <- readr::read_rds('./dat/ftb/ASM1462v1.Rds')

dir.create('./dat/grn', F, T)

# calculate relative tRNA concentration for each codon -------------------------

cdn_freq_list_pa01 <- NumericList(sapply(cdn_pa01, function(codon) 1 / unlist(freq[codon], use.names = F)))
cdn_freq_list_pa14 <- NumericList(sapply(cdn_pa14, function(codon) 1 / unlist(freq[codon], use.names = F)))

# assign relative tRNA frequencies to granges objects for both strains ---------

cds_pa01$rel_tRNA_conc <- cdn_freq_list_pa01
cds_pa14$rel_tRNA_conc <- cdn_freq_list_pa14


require(KEGGREST)

orgs <- keggList('organism')
orgs <- orgs[grepl('Pseudomonas aeruginosa .*', orgs[, 'species']), ]
crossref <- keggLink('pau', 'enzyme')
path14 <- keggLink('path', 'pau')
path01 <- keggLink('path', 'pae')
