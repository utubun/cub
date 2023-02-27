# class definition -------------------------------------------------------------

#' An S4 class to represent a Codon Usage Bias of a given DNA / Set of DNAs
#'
#' @slot id A length-one character vector, giving the identifier to your Cub
#' @slot data A character vector of codons on position 1, 4, 7, ... etc
#' @export
#' @examples
#' myCub <- Cub('A cub', PA14CDS[[1]])
#' myCub
setClass(
  'Cub',
  slots = c(
    id    = 'character',
    data  = 'character'
  ),
  prototype = list(
    id    = NA_character_,
    data  = NA_character_
  )
)

#' An S4 class to represent a Codon Usage Bias of a given Set of DNAs
#'
#' @slot data A list of character vectors of codons on position 1, 4, 7, ... etc
#' for the given DNAStringSet
#' @export
setClass(
  'CubSet',
  slots = c(
    data  = 'list'
  )
)

#' An S4 class to represent a Relative Synonymous Codon Usage counts for the given
#' DNAString
#'
#' @slot data A named list of numerics, representing RSCU for a given codon
#' @export
setClass(
  'Counts',
  slots = c(
    data  = 'list'
  )
)
