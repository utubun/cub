## code to prepare `BSgenome.Ecoli.NCBI.K12.MG1655` dataset goes here
forgeBSgenomeDataPkg('Bsgenome.Ecoli.NCBI.K12.MG1655-seed')
system('R CMD build BSgenome.Ecoli.NCBI.K12.MG1655')
system('R CMD INSTALL BSgenome.Ecoli.NCBI.K12.MG1655_1.0.0.tar.gz')
usethis::use_data(BSgenome.Ecoli.NCBI.K12.MG1655, overwrite = TRUE)
