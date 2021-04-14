RiboTracker
===========
Dmytro Strunin
March, 2020

# Objectives

## Genome and Genomic Features annotation

- [x] Build `BSGenome` for *Pseudomonas aeruginosa P01* strain;
    - [ ] Resolve issues with *common organism name*;
- [x] Build `BSGenome` for *Pseudomonas aeruginosa P14* strain;
    - [ ] Resolve issues with *common organism name*;
- [x] Make both genomes standalone packages, which can be installed by standard R routine;
- [x] Build `TxDb` for *Pseudomonas aeruginosa P01* strain;
- [x] Build `TxDb` for *Pseudomonas aeruginosa P14* strain;

Both genomes can be installed (from source) and loaded. For loading into R session use:

```{r}
library(BSgenome.Paeruginosa.NCBI.ASM676v1)
library(BSgenome.Paeruginosa.NCBI.ASM1462v1)

gnm_pa01 <- BSgenome.Paeruginosa.NCBI.ASM676v1
gnm_pa14 <- BSgenome.Paeruginosa.NCBI.ASM1462v1
```

`TxDb` annotations require `GenomicFeatures` to be installed an loaded. Load databases as:

```{r}
txdb_pa01 <- loadDb('./dat/gnm/sql/TxDb.Paeruginosa.ASM676v1.sqlite')
txdb_pa14 <- loadDb('./dat/gnm/sql/TxDb.Paeruginosa.ASM1462v1.sqlite')
```
When package is build, load using `data()` routine.

**Issue**: genomes and `TxDb` objects have different sequence names (need to be fixed). Current solution:

```{r}
seqlevels(txdb_object) <- seqnames(genome_object)
```

## CDS extraction

- [x] Extract DNA sequences of *Pseudomonas aeruginosa P01* CDSs;
- [x] Extract DNA sequences of *Pseudomonas aeruginosa P14* CDSs;
- [ ] Check transcripts/CDS on both strands, just to be sure;
- [ ] Chek if sequence extraction follows the CDS orientation (+/- strand);

Since we are interested in relative translation speed, we can extract genomic coordinates only for `CDS`s from `TxDb` objects for PA01 and PA14 respectively, excluding irrelevant sequences of non-coding genes:

```{r}
cds_pa01 <- cds(txdb_pa01, use.names = T)
cds_pa14 <- cds(txdb_pa14, use.names = T)
```
Code listed above, returns `GRanges` objects, with *genomic coordinates* of all CDSs, strand orientation (**+**/**-**) and CDS's ID in metadata.

Using genomes and genomic features databases for each particular strain, we can extract CDSs as follows:

```{r}
dna_pa01 <- getSeq(gnm_pa01, cds_pa01)
dna_pa14 <- getSeq(gnm_pa14, cds_pa14)
```
`getSeq()` extracts sequence taking DNA strand orientation (+/-) into account. Just checking that it is correct:

```{r}
strand(cds_pa01[1])
strand(cds_pa01[5036])
(cds_01_coords <- c(start(cds_pa01[1]), end(cds_pa01[1])))
(cds_5036_coords <- c(start(cds_pa01[5036]), end(cds_pa01[5036])))
gnm_pa01$chr[cds_01_coords[1]:cds_01_coords[2]]
dna_pa01[1]
gnm_pa01$chr[cds_5036_coords[1]:cds_5036_coords[2]]
dna_pa01 <- gnm_pa01$chr
dna_pa01[5036]
```
```
// CDS 1, + strand
factor-Rle of length 1 with 1 run
  Lengths: 1
  Values : +
Levels(3): + - *
// CDS 5036 - strand
factor-Rle of length 1 with 1 run
  Lengths: 1
  Values : -
Levels(3): + - *
// Genomic coordinates of CDS 1
[1]  483 2027
//Genomic coordinates of CDS 5036
[1] 5162472 5163476
// Genomic DNA extracted by the coordinates, + strand
1545-letter "DNAString" instance
seq: GTGTCCGTGGAACTTTGGCAGCAGTGCGTGGATCTTCTCCGCGATGAGCTG...CGGATATCCGCGAGGACTACAAGAACCTGCTGCGTACCCTGACAACCTGA
// DNA extracted by GRanges + strand
A DNAStringSet instance of length 1
  width seq                                                                            names               
[1]  1545 GTGTCCGTGGAACTTTGGCAGCAGTGCGTGGATCTTCT...GGACTACAAGAACCTGCTGCGTACCCTGACAACCTGA NP_064721.1
// DNA extracted by the coordinates, + strand
1005-letter "DNAString" instance
seq: TCAGCGGCCTTCGACGATGCGTTCGAAGCCCTCGCGGATACTGTCCTCCGG...AGCACGGTGACGGGGATCGGGGCGTGGTTGTCGTTGGCTACTTCGGACAT
// DNA extracted by GRanges - strand
A DNAStringSet instance of length 1
  width seq                                                                            names               
[1]  1005 ATGTCCGAAGTAGCCAACGACAACCACGCCCCGATCCC...CCGCGAGGGCTTCGAACGCATCGTCGAAGGCCGCTGA NP_253294.1
// Reverse complement
```
### Split each CDS for both strains into triplets (codons)

- [ ] Store the triplets, or triplet names in `GRanges` object relevant to particular strain;
- [ ] Find a way to efficiently store/retrieve this information (currently: too much memory, relatively slow);

Since `codons()` function from `Biostrings` package is very slow, I applied considerably more efficient version:

```{r}
triplets <- function(x) {
  unlist(strsplit(toString(x), '(?<=.{3})', perl = TRUE), use.names = F)
}
```
Function above converts `DNAString` object into `character` string, and split it triplet-by-triplet, returning character vector. With `parLapply` from `parallel` it returns named list of codons for each CDS in `dna_pa01` and `dna_pa14` CDS sets:

```{r}
ncores <- detectCores()
cl     <- makeCluster(ncores)

cdn_pa01 <- parLapply(cl, dna_pa01, triplets)
cdn_pa14 <- parLapply(cl,  dna_pa14, triplets)

stopCluster(cl)
```
