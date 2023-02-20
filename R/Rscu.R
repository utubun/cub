# class definition -------------------------------------------------------------
setClass(
  'RSCU',
  representation(
    data  = 'list'
  )
)

RSCU <- function(x) {

  if (!all(names(x) %in% names(Biostrings::GENETIC_CODE))) {
    stop('Wrong genetic code!')
  }

  if (any(!is.numeric(x)) | any(x < 0)) {
    stop('Only positive numeric counts allowed!')
  }

  if (any(is.na(x) | any(is.null(x)))) {
    warning('Replacing NA/NULL values with 0!')
    x[is.na(x) | is.null(x)] <- 0
  }

  aac_lookup <- split(names(GENETIC_CODE), GENETIC_CODE)

  res <- lapply(aac_lookup, function(cdns) {
    x[cdns] * length(cdns) / sum(x[cdns])
  })

  return(
    new('RSCU', data = res)
  )
}

setGeneric("as.data.frame", function(object, method, window) {
  standardGeneric("as.data.frame")
})

setGeneric("plot", function(object, method, window) {
  standardGeneric("plot")
})

setMethod('as.data.frame', signature('RSCU'), function(object) {
  do.call(
    rbind,
    lapply(
      names(object@data),
      function(nm) {
        data.frame(
          aa    = nm,
          codon = names(object@data[[nm]]),
          rscu  = object@data[[nm]]
        )
      }
    )
  )
})

setMethod('plot', signature('RSCU'), function(object) {

  dat <- as.data.frame(object) |>
    dplyr::filter(rscu > 0, aa != '*') |>
    dplyr::arrange(aa, dplyr::desc(rscu)) |>
    dplyr::mutate(aa = forcats::fct_reorder(aa, rscu, sum))

  ggplot2::ggplot(
      data = dat,
      mapping = ggplot2::aes(x = aa, y = rscu, label = codon)
    ) +
    ggplot2::geom_bar(stat = 'identity', position = 'stack', fill = 'gray35', color = 'gray85') +
    ggplot2::geom_text(
      stat = 'identity',
      position = ggplot2::position_stack(vjust = .5),
      size = 2.5,
      color = 'gray75'
    ) +
    ggplot2::scale_fill_grey() +
    ggplot2::coord_flip() +
    labs(x = 'Amino acid', y = 'RSCU')
})

