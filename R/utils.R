#' @export
#'
Cub <- function(id, x) {

  if(missing(id)) {
    stop('Argument id is rquired')
  }

  cdns <- tryCatch(
    {
      as.character(Biostrings::codons(x))
    },
    error = function(cond) {
      NA_character_
    }
  )

  return(
    new('Cub', id = id, data = cdns)
  )

}

#' Creates an instance of CubSet
#' @export
CubSet <- function(x) {
  return(new('CubSet', data = lapply(names(x), \(nm) { Cub(nm, x[[nm]]) })))

}

#' Creates a random DNA
#' @export
dnar <- function(n) {
  chr <- paste0('GTC', sample(names(Biostrings::GENETIC_CODE), (n - 1) %/% 3, T), collapse = '')
  dna <- Biostrings::DNAString(chr)
  return(dna)
}
