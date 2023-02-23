cub <- function(id, x) {
  if(missing(id)) {
    stop('Argument id is rquired')
  }

  if(!class(x) %in% c('character', 'DNAString')) {
    stop('Argument seq must be provided as character or DNAString')
  }

  if(class(x) == 'character') {
    x = tryCatch(
      {
        Biostrings::DNAString(x)
      },
      error = function(cond) {
        stop(sprintf('The string provided, contains non-DNA symbols. \\nError: %s', cond))
      }
    )

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
