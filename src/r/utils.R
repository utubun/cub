
fun <- function(dna, ...) {
  res <- tryCatch(
    {
      as.character(codons(dna))
    },
    error = function(cond) {
      return(NA)
    }
  )

  return(res)
}

cdnCount <- function(x) {
  ncores <- parallel::detectCores()
  cl     <- makeCluster(ncores - 1)

  res <- parLapply(cl, x, function(datum) {
    tb <- table(datum)
    if(length(tb)) {
      tb <- tibble::tibble(codon = names(tb), count = as.numeric(tb))
    } else {
      tb <- list(codon = NA_character_, count = NA_integer_)
    }
    return(tb)
  })

  stopCluster(cl)
  rm(cl)

  data.table::rbindlist(res, fill = TRUE, idcol = 'gene')
}
