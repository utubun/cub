#' Creates a Counts instance from a given Cub object
#' @param object A Cub instance to summarize
#' @return Counts representing Relative Synonymous Codon Usage for a given Gene / Genome
#' @export
#' @examples
#' dat <- Cub(PAO14CDS[[1]])
#' rscu(dat)
setGeneric('rscu', function(object, ...) {
  standardGeneric('rscu')
})

#' Split DNAString / DNAString set into codons
#' @param character DNAString or DNAStringSet representing the Gene or Genome
#' @return character vector of codones, or list of characters for DNAStringSet
#' @export
#' @examples
#' dat <- PAO14CDS[[1]]
#' codons(dat)
setGeneric("codons", function(object, ...) {
  standardGeneric("codons")
})

#' Creates a RAC instance from a given Cub object
#' @param object A Cub instance to summarize
#' @return RAC representing Relative Adaptiveness of a Codon for a given Gene / Genome
#' @export
#' @examples
#' dat <- Cub(PAO14CDS[[1]])
#' rac(dat)
setGeneric('rac', function(object, ref, ...) {
  standardGeneric('rac')
})

#' Calculates CAI value for a given gene
#' @param object A Cub instance to summarize
#' @return RAC representing Relative Adaptiveness of a Codon for a given Gene / Genome
#' @export
#' @examples
#' dat <- Cub(PAO14CDS[[1]])
#' cai(dat)
setGeneric('cai', function(x, ref, ...) {
  standardGeneric('cai')
})

#' Counts codons' occurrence in a given Cub instance
#' @param object A Cub instance to summarize
#' @return A named integer vector representing count of a given ccodon in a DNAString
#' @export
#' @examples
#' dat <- Cub('A DNA', PA14CDS[[1]])
#' count(dat)
setGeneric('count', function(object, ...) {
  standardGeneric('count')
})

#' Calculates CAI value for a given gene
#' @param object A Cub instance to summarize
#' @return RAC representing Relative Adaptiveness of a Codon for a given Gene / Genome
#' @export
#' @examples
#' dat <- Cub(PAO14CDS[[1]])
#' cai(dat)
setGeneric('cplot', function(x, ...) {
  standardGeneric('cplot')
})
