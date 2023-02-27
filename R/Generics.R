#' Shows summary of a Cub object.
#' @param object A Cub instance to summarize
#' @return NULL
#' @export
#' @examples
#' dat <- Cub(PAO14CDS[[1]])
#' summary(dat)
setGeneric("summary", function(object, ...) {
  standardGeneric("summary")
})

#' Prints a Cub object.
#' @param object A Cub instance to summarize
#' @return NULL
#' @export
#' @examples
#' dat <- Cub(PAO14CDS[[1]])
#' print(dat)
setGeneric("print", function(object, ...) {
  standardGeneric("print")
})

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

#' Converts Counts instance into data.frame
#' @param Counts object to be converted
#' @return data.frame with columns aa: Amino Acid, codon: Codon code, rscu: Counts x(i, j)
#' for a given codon.
#' @export
#' @examples
#' dat <- rscu(PAO14CDS[[1]])
#' as.data.frame(dat)
setGeneric("as.data.frame", function(object, method, window) {
  standardGeneric("as.data.frame")
})

#' Plots Counts as a stucked bar chart, represented RSCU(i, j) for each codon, grouped by
#' amino acid.
#' @param x A Counts object to be plotted
#' @return A ggplot2::ggplot object
#' @export
#' @examples
#' dat <- rscu(PAO14)
#' plot(dat)
setGeneric("plot", function(object, value, scale, ...) {
  standardGeneric("plot")
})

#' Creates a RAC instance from a given Cub object
#' @param object A Cub instance to summarize
#' @return RAC representing Relative Adaptiveness of a Codon for a given Gene / Genome
#' @export
#' @examples
#' dat <- Cub(PAO14CDS[[1]])
#' rac(dat)
setGeneric('rac', function(object, ...) {
  standardGeneric('rac')
})

#' Calculates CAI value for a given gene
#' @param object A Cub instance to summarize
#' @return RAC representing Relative Adaptiveness of a Codon for a given Gene / Genome
#' @export
#' @examples
#' dat <- Cub(PAO14CDS[[1]])
#' cai(dat)
setGeneric('cai', function(object, ...) {
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
