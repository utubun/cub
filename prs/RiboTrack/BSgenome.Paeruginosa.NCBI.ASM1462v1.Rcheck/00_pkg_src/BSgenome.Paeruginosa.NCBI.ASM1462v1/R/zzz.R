###
###

.pkgname <- "BSgenome.Paeruginosa.NCBI.ASM1462v1"

.seqnames <- 'CP000438.1'

.circ_seqs <- 'CP000438.1'

.mseqnames <- NULL

.onLoad <- function(libname, pkgname)
{
    if (pkgname != .pkgname)
        stop("package name (", pkgname, ") is not ",
             "the expected name (", .pkgname, ")")
    extdata_dirpath <- system.file("extdata", package=pkgname,
                                   lib.loc=libname, mustWork=TRUE)

    ## Make and export BSgenome object.
    bsgenome <- BSgenome(
        organism="Pseudomonas aeruginosa",
        common_name="P. aeruginosa PA14",
        genome="ASM1462v1",
        provider="NCBI",
        release_date="2013-10-06",
        source_url="https://www.ncbi.nlm.nih.gov/assembly/GCF_000014625.1",
        seqnames=.seqnames,
        circ_seqs=.circ_seqs,
        mseqnames=.mseqnames,
        seqs_pkgname=pkgname,
        seqs_dirpath=extdata_dirpath
    )

    ns <- asNamespace(pkgname)

    objname <- pkgname
    assign(objname, bsgenome, envir=ns)
    namespaceExport(ns, objname)

    old_objname <- "Paeruginosa"
    assign(old_objname, bsgenome, envir=ns)
    namespaceExport(ns, old_objname)
}

