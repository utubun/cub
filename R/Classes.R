# class definition -------------------------------------------------------------

setClass(
  'RiboTrack',
  representation(
    info  = 'character',
    freq  = 'CompressedNumericList'
    ),
  prototype(
    info  = NA_character_,
    freq  = NumericList()
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

# methods ----------------------------------------------------------------------

setMethod("[", signature("RiboTrack"), function(x, i) {
  new(
    "RiboTrack",
    info   = x@info,
    freq   = x@freq[i]
  )
})

setMethod("info", signature(object = "RiboTrack"), function(object) {
  object@info
})

setMethod("names", signature("RiboTrack"), function(x) {
  cat(
    x@info,
    head(names(x@freq))
  )
  invisible(names(x@freq))
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


# # ------------------------------------------------------------------------------
#
# pa14 <- new('RiboTrack')
# slotNames(pa14)
# pa14@info <- 'P.aeruginosa, PA14'
# pa14@freq <- cdn_freq_list_pa14
# info(pa14)
# x <- pa14[1]
# plot(x, 1)
# plot(pa14, i = 333, col = 'blue')
# plot(pa14, i = 'ABJ14969.1')
# x <- x - 1
# plot(x, i = 1, col = 'salmon')
# x <- cumsum(x)
# plot(x, i = 1, col = 'salmon')
# ctr = mean(unlist(pa14@freq))
#
# pa14_scaled <- scale(pa14, ctr, FALSE)
# pa14_sm_spline <- smooth.spline(pa14_scaled)
# pa14_spline <- spline(pa14_scaled, 100)
# pa14_cumsum <- cumsum(pa14_spline)

