# class definition -------------------------------------------------------------

setClass(
  'RiboTrack',
  representation(organism = 'character', strain = 'character', granges = 'GRanges'),
  prototype(organism = NA_character_, strain = NA_character_, granges = GRanges())
  )

# generic functions ------------------------------------------------------------

setGeneric("get_organism", function(object) {
  standardGeneric("get_organism")
})

setGeneric("get_strain", function(object) {
  standardGeneric("get_strain")
})

setGeneric("get_granges", function(object) {
  standardGeneric("get_granges")
})

setGeneric("get_cds_names", function(object, txdb) {
  standardGeneric("get_cds_names")
})

setGeneric("get_gene_names", function(object, txdb) {
  standardGeneric("get_gene_names")
})

setGeneric("interpolate", function(object, method, window) {
  standardGeneric("interpolate")
})

# methods ----------------------------------------------------------------------

setMethod("[", signature("RiboTrack"), function(x, i) {
  new(
    "RiboTrack",
    organism = x@organism,
    strain   = x@strain,
    granges  = x@granges[i]
  )
})

setMethod("get_organism", signature(object = "RiboTrack"), function(object) {
  object@organism
})

setMethod("get_strain", signature(object = "RiboTrack"),   function(object) {
  object@strain
})

setMethod("get_granges", signature(object = "RiboTrack"),   function(object) {
  object@granges
})

setMethod("get_cds_names", signature(object = "RiboTrack"), function(object, txdb) {
  select(txdb, as.character(object@granges$cds_id), 'CDSNAME', 'CDSID')[['CDSNAME']]
})

setMethod("get_gene_names", signature(object = "RiboTrack"), function(object, txdb) {
  select(txdb, as.character(object@granges$cds_id), 'GENEID', 'CDSID')[['GENEID']]
})

setMethod(
  "interpolate",
  signature(object = "RiboTrack"),
  function(
    object,
    method = c('none', 'geometric', 'average', 'harmonic'),
    window = 25
    ) {
    method <- match.arg(method)
    f = switch(
      method,
      none      = I,
      geometric = geometric.mean,
      average  = average,
      harmonic = harmonic.mean
    )

    x = object@granges$rel_tRNA_conc

    res <- sapply(x, function(i) {
      tmp <- data.table::frollapply(i, window, f)
      replace(tmp, is.na(tmp), f(i[is.na(tmp)]))
      })
    sp  <- sapply(res, splinefun)

    h.mean = sapply(res, harmonic.mean)
    a.mean = sapply(res, average)
    g.mean = sapply(res, geometric.mean)

    return(list(interpolation = res, sp_fun = sp, h.mean = h.mean, a.mean = a.mean, g.mean = g.mean))
  }
)

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


# ------------------------------------------------------------------------------

pa14 <- new('RiboTrack')
slotNames(pa14)
pa14@organism <- 'P.aeruginosa'
pa14@strain <- 'PA14'
pa14@granges <- sort(cds_pa14)
get_organism(pa14)
get_strain(pao1)
get_granges(pao1)

pa01 <- new('RiboTrack')
slotNames(pa01)
pa01@organism <- 'P.aeruginosa'
pa01@strain <- 'PA01'
pa01@granges <- sort(cds_pa01)
get_organism(pa01)
get_strain(pa01)
get_granges(pa01)

pa01_names <- get_gene_names(pa01, txdb_pa01)
pa14_names <- get_gene_names(pa14, txdb_pa14)

idx01 <- which(pa01_names %in% pa14_names)
idx14 <- which(pa14_names %in% pa01_names[idx01])

idx01 <- idx01[order(pa01_names)]
idx14 <- idx14[order(pa14_names)]

pa01_names <- sort(pa01_names)
pa14_names <- sort(pa14_names)

pa01.1 <- pa01[idx_pa01]
pa14.1 <- pa14[idx_pa14]

res01 <- interpolate(pa01.1, 'average', 25)
res14 <- interpolate(pa14.1, 'average', 25)

plot_res <- function(res, which, name = '') {
  x = seq_along(res[[1]][[which]])
  y = res[[1]][[which]]
  f = res[[2]][[which]]
  x_poly = c(1, x, length(x), 1)
  y_poly = c(1, f(x), 1, 1)
  plot(x, y, pch = '.', col = 'gray10', ylim = c(0, 5), xlab = 'Sequence position', ylab = 'Relative ribosome speed', main = name)
  curve(f, min(x), max(x), col = 'blue', add = T)
  polygon(x_poly, y_poly, col = 'red')
  abline(h = res[[4]][[which]], lty = 3, col = 'gray25')
  abline(h = 1, lty = 3, col = 'red')
}

ttest <- sapply(seq_along(pa01_names), function(i) t.test(res14[[1]][[i]], res01[[1]][[i]], alternative = 'less')$p.value <= 0.05)
mean(ttest)
testidx <- which(ttest)

i = testidx[5]
par(mfrow = c(2, 1))
plot_res(res01, i, pa01_names[i])
plot_res(res14, i, pa14_names[i])

plot(seq_along(res[[1]][[15]]), res[[1]][[15]], pch = '*', col = 'gray25', ylim = c(0, 2))
f = res[[2]][[15]]

curve(f, 0, 350, col = 'blue', add = T)
abline(h = c(res[[3]][[15]], res[[4]][[15]], res[[5]][[15]]), lty = 3, col = c('red', 'green', 'blue'))
abline(h = 1, col = 'red', lty = 3)
polygon(
  c(1, seq_along(res[[1]][[15]]), length(res[[1]][[15]]), 1),
  c(1, res[[2]][[15]](seq_along(res[[1]][[15]])), 1, 1),
  col = 'red'
)
