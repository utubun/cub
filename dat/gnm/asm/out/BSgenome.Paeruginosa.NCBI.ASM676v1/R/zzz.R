###
###

.pkgname <- "BSgenome.Paeruginosa.NCBI.ASM676v1"

.seqnames <- 'AE004091.2'

.circ_seqs <- 'AE004091.2'

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
        organism="Pseudomonas aeruginosa (P. aeruginosa PAO1)",
        common_name="P. aeruginosa PA01",
        genome="ASM676v1",
        provider="NCBI",
        release_date="2006-07-07",
        source_url="https://www.ncbi.nlm.nih.gov/assembly/GCF_000006765.1",
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

