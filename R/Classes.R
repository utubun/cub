# class definition -------------------------------------------------------------
setClass(
  'Cub',
  representation(
    info  = 'character',
    count = 'data.frame'
  ),
  prototype(
    info  = NA_character_,
    count = data.frame()
  )
)

# generic functions ------------------------------------------------------------

setGeneric("info", function(object) {
  standardGeneric("info")
})

setGeneric("print", function(object) {
  standardGeneric("print")
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

setAs('DNAStringSet', 'Cub', function(from, to) {
  cdns <- get_codons(from)
  new(
    'Cub',
    info = NA_character_,
    count = cdnCount(cdns)

  )
})

setMethod('cub', signature('Cub'), function(object, method = 'total') {
  if(method == 'total') {
    res <- object@count %>%
      group_by(gene) %>%
      mutate(frequency = count / sum(count)) %>%
      dplyr::select(gene, codon, frequency)
  }
  return(res)
})

setMethod("[", signature("RiboTrack"), function(x, i) {
  new(
    "RiboTrack",
    info   = x@info,
    freq   = x@freq[i]
  )
})

setMethod("info", signature(object = "Cub"), function(object) {
  object@info
})

setMethod("names", signature("Cub"), function(x) {
  cat(
    x@info,
    x@gene,
    '\n'
  )
  invisible(x@gene)
})

setMethod('-', signature('RiboTrack', 'numeric'), function(e1, e2) {
  new(
    'RiboTrack',
    info = e1@info,
    freq = NumericList(sapply(e1@freq, '-', y = e2))
  )
})

setMethod('+', signature('RiboTrack', 'numeric'), function(e1, e2) {
  new(
    'RiboTrack',
    info = e1@info,
    freq = NumericList(sapply(e1@freq, '+', y = e2))
  )
})

setMethod('*', signature('RiboTrack', 'numeric'), function(e1, e2) {
  new(
    'RiboTrack',
    info = e1@info,
    freq = NumericList(sapply(e1@freq, '*', y = e2))
  )
})

setMethod('/', signature('numeric', 'RiboTrack'), function(e1, e2) {
  new(
    'RiboTrack',
    info = e2@info,
    freq = NumericList(sapply(e2@freq, function(x) e1 / x))
  )
})

setMethod('scale', signature('RiboTrack'), function(x, center = TRUE, scale = FALSE) {
  freq <- NumericList(
    sapply(x@freq, function(i) scale(i, center = center, scale = scale))
  )
  new(
    'RiboTrack',
    info = x@info,
    freq = freq
  )
})

setMethod('cumsum', signature('RiboTrack'), function(x) {
  new(
    'RiboTrack',
    info = x@info,
    freq = NumericList(lapply(x@freq, cumsum))
  )
})

setMethod('spline', signature('RiboTrack'), function(x, n = 100) {
  y = lapply(
    x@freq,
    function(freq) {
      spline(seq(0, 100, length.out = length(freq)), freq, n = n)$y
    }
  )
  new(
    'RiboTrack',
    info = x@info,
    freq = NumericList(y)
  )
})

setMethod('smooth.spline', signature('RiboTrack'), function(x) {
  y = lapply(
    x@freq,
    function(freq) {
      smooth.spline(seq(0, 100, length.out = length(freq)), freq, nknots = 100)
    }
  )
  as.matrix(lapply(y, function(i) i$fit$coef))
})

setMethod('plot', signature('RiboTrack'), function(x, i, ...) {
  y = unlist(x[i]@freq)
  x = seq_along(y)
  plot(x, y, xlab = 'AA/Codone position', ylab = 'tRNA relative concentration', type = 'l', ...)
})

average <- function(x){
  x[is.infinite(x)] <- median(x[!is.infinite(x)])
  return(mean(x, na.rm = T))
}

geometric.mean <- function(x) {
  return(prod(x, na.rm = T)^(1/length(na.omit(x))))
}

harmonic.mean <- function(x) {
  x[x == 0] <- median(x, na.rm = T)
  x[is.infinite(x)] <- median(x[!is.infinite(x)], na.rm = T)
  return(1 / mean(1 / x, na.rm = T))
}


