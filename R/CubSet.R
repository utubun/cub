# class definition -------------------------------------------------------------
setClass(
  'CubSet',
  representation(
    data  = 'list'
  )
)

CubSet <- function(x) {

  new('CubSet', data = lapply(names(x), function(nm) { cub(nm, x[[nm]]) }))

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

setMethod(summary, signature('CubSet'), function(object, ...) {
  lapply(object@data, summary)
})

setMethod('print', signature('CubSet'), function(object, ...) {
  lapply(object@data, print)
})

setMethod('show', signature('CubSet'), function(object) {
  lapply(object@data, show)
})

setMethod('rscu', signature('CubSet'), function(object, each = TRUE) {
  standardGeneric("rscu")
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

setMethod('rscu', signature('CubSet'), function(object, each = FALSE) {

  if (each) {
    return(lapply(object@data, rscu))
  } else {
    res <- unlist(lapply(object@data, function(datum) { datum@data }))
    res <- table(res)
    res <- setNames(as.vector(res), names(res))
    return(RSCU(res))
  }
})
