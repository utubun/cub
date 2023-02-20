# class definition -------------------------------------------------------------
setClass(
  'Cub',
  representation(
    id    = 'character',
    data  = 'character'
  ),
  prototype(
    id    = NA_character_,
    data  = NA_character_
  )
)

Cub <- function(id, x) {
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
      c(NA_character_, attr = cond)
    }
  )

  return(
    new('Cub', id = id, data = cdns)
  )

}

# generic functions ------------------------------------------------------------

setGeneric("summary", function(object, ...) {
  standardGeneric("summary")
})

setGeneric("print", function(object, ...) {
  standardGeneric("print")
})

setGeneric("show", function(object, ...) {
  standardGeneric("show")
})

setGeneric('rscu', function(object, ...) {
  standardGeneric('rscu')
})

setGeneric("interpolate", function(object, method, window) {
  standardGeneric("interpolate")
})

setGeneric("codons", function(object, ...) {
  standardGeneric("codons")
})

setGeneric("cub", function(object, ...) {
  standardGeneric("cub")
})

# methods ----------------------------------------------------------------------

setMethod(summary, signature('Cub'), function(object, ...) {

  counts <- table(object@data)
  counts <- sort(setNames(as.vector(counts), names(counts)), decreasing = TRUE)
  counts <- sprintf('%s (%d)', names(counts), as.vector(counts))

  total  <- sprintf(' Total:    %d\n', length(object@data))
  cuniq  <- sprintf('Unique:   %d\n', length(unique(object@data)))
  domin  <- sprintf('Frequent: %s\n', paste0(head(counts, 5), collapse = ',\t'))
  rare   <- sprintf('Rare:     %s\n', paste0(rev(tail(counts, 5)), collapse = ',\t'))

  cat(total, cuniq, domin, rare)
})

setMethod('print', signature('Cub'), function(object, ...) {
  cat(
    sprintf('ID:\t%s\n', object@id),
    '\n',
    summary(object)
  )
})

setMethod('show', signature('Cub'), function(object) {
  print(object)
})

setMethod('rscu', signature('Cub'), function(object) {

  count <- table(object@data)
  count <- setNames(as.vector(count), names(count))

  res <- setNames(
    vector('numeric', length(Biostrings::GENETIC_CODE)),
    names(Biostrings::GENETIC_CODE)
  )

  res[names(count)] <- count

  return(RSCU(res))
})

setMethod(codons, signature('DNAStringSet'), function(object, ...) {
  # Tests for 'DNAStringSet or DNAString

  res <- lapply(object, function(dna) {

    cdns <- tryCatch(
      {
        as.character(Biostrings::codons(dna))
      },
      error = function(cond) {
        c(NA_character_, attr = cond)
      }
    )

    if(length(cdns > 1)) {
      count <- table(cdns)
      return(setNames(as.numeric(count), names(count)))
    }

    return(cdns)
  })

  return(res)

})
