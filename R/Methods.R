#' Shows summary of a Cub object.
#' @param object A Cub instance to summarize
#' @return NULL
#' @export
#' @examples
#' dat <- Cub('A DNA', PAO1CDS[[1]])
#' summary(dat)
setMethod('summary', signature('Cub'), function(object, ...) {

  counts <- table(object@data)
  counts <- sort(setNames(as.vector(counts), names(counts)), decreasing = TRUE)
  counts <- sprintf('%s (%d)', names(counts), as.vector(counts))

  total  <- sprintf(' Total:    %d\n', length(object@data))
  cuniq  <- sprintf('Unique:   %d\n', length(unique(object@data)))
  domin  <- sprintf('Frequent: %s\n', paste0(head(counts, 5), collapse = ',\t'))
  rare   <- sprintf('Rare:     %s\n', paste0(rev(tail(counts, 5)), collapse = ',\t'))

  cat(total, cuniq, domin, rare)
})

#' Shows a Cub object.
#' @param object A Cub instance to summarize
#' @return NULL
#' @export
#' @examples
#' dat <- Cub('A DNA', PAO1CDS[[1]])
#' dat
setMethod('show', signature('Cub'), function(object) {
  object
})

#' Split DNAString / DNAString set into codons
#' @param character DNAString or DNAStringSet representing the Gene or Genome
#' @return character vector of codones, or list of characters for DNAStringSet
#' @export
#' @examples
#' dat <- PA14CDS[[1]]
#' codons(dat)
setMethod('codons', signature('DNAStringSet'), function(object, ...) {
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

#' Counts codons' occurrence in a given Cub instance
#' @param object A Cub instance to summarize
#' @return A named integer vector representing count of a given ccodon in a DNAString
#' @export
#' @examples
#' dat <- Cub('A DNA', PA14CDS[[1]])
#' count(dat)
setMethod(count, signature('Cub'), function(object) {

  cnts <- table(object@data)
  cnts <- setNames(as.vector(cnts), names(cnts))

  res <- setNames(
    vector('numeric', length(Biostrings::GENETIC_CODE)),
    names(Biostrings::GENETIC_CODE)
  )

  res[names(cnts)] <- cnts

  res <- split(res, Biostrings::GENETIC_CODE)

  return(new('Counts', data = res))
})


#' Creates a Count instance from a given Cub object
#' @param object A Cub instance to summarize
#' @return RSCU representing Relative Synonymous Codon Usage for a given Gene / Genome
#' @export
#' @examples
#' dat <- Cub('A DNA', PA14CDS[[1]])
#' rscu(dat)
setMethod('rscu', signature('Cub'), function(object) {

  cnts <- count(object)

  res <- lapply(cnts@data, function(datum) {
    datum * length(datum) / sum(datum)
  })

  return(
    new('Counts', data = res)
  )
})

#' Creates a Count instance from a given Cub object
#' @param object A Cub instance to summarize
#' @return RSCU representing Relative Synonymous Codon Usage for a given Gene / Genome
#' @export
#' @examples
#' dat <- Cub('A DNA', PA14CDS[[1]])
#' rscu(dat)
setMethod('rac', signature('Cub'), function(object, ref) {

  cnts_obs <- rscu(object)
  cnts_max <- sapply(ref@data, max, na.rm = TRUE)

  res <- sapply(names(cnts_obs@data), function(nm) {
    cnts_obs@data[[nm]] / cnts_max[[nm]]
  })

  names(res) <- names(cnts_obs@data)

  return(
    new('Counts', data = res)
  )
})

#' Calculates CAI value for a given gene
#' @param object A Cub instance to summarize
#' @return A numeric value for CAI
#' @export
#' @examples
#' dat <- Cub('A DNA', PA14CDS[[1]])
#' cai(dat)
setMethod('cai', signature('Cub'), function(x, ref) {

  w <- rac(x, ref)
  w <- unlist(setNames(w@data, NULL))
  w <- replace(w, w == 0, .5)

  w <- w[x@data]

  res <- (prod(w))^(1 / length(w))
  #res <- exp(1 / !is.na(w) * sum(log(w[!is.na(w)])))
  return(res)
})

#' Converts Counts instance into data.frame
#' @param Counts object to be converted
#' @return data.frame with columns aa: Amino Acid, codon: Codon code, rscu: RSCU x(i, j)
#' for a given codon.
#' @export
#' @examples
#' Cub('A DNA', PA14CDS[[1]]) |>
#'   rscu() |>
#'   as.data.frame() |>
#'   head()
setMethod('as.data.frame', signature('Counts'), function(x, ...) {
  res <- do.call(
    rbind,
    lapply(
      names(x@data),
      function(nm) {
        data.frame(
          aa    = nm,
          codon = names(x@data[[nm]]),
          val   = x@data[[nm]]
        )
      }
    )
  )
  rownames(res) <- seq(nrow(res))
  return(res)
})

#' Plots Counts as a stucked bar chart, represented RSCU(i, j) for each codon, grouped by
#' amino acid.
#' @param x A Counts object to be plotted
#' @return A ggplot2::ggplot object
#' @export
#' @examples
#' Cub('A DNA', PA14CDS[[1]]) |>
#'   rscu() |>
#'   plot()
setMethod('cplot', signature('Counts'), function(x, ...) {

  args <- list(...)

  scale <- args[['scale']]

  dat <- as.data.frame(x) |>
    dplyr::filter(val > 0, aa != '*')

  # if (!is.null(scale) & (scale == 'all' | scale == TRUE)) {
  #   dat <- dat |>
  #     dplyr::mutate(val = val / sum(val)) |>
  #     dplyr::arrange(aa, dplyr::desc(val)) |>
  #     dplyr::mutate(aa = forcats::fct_reorder(aa, val, sum))
  # } else if (!is.null(scale) & scale == 'aa') {
  #   dat <- dat |>
  #     dplyr::group_by(aa) |>
  #     dplyr::mutate(val = val / sum(val), maxval = max(val)) |>
  #     dplyr::ungroup() |>
  #     dplyr::mutate(aa = forcats::fct_reorder(aa, maxval), maxval = NULL) |>
  #     dplyr::arrange(aa, desc(val))
  # } else {
    dat <- dat |>
      dplyr::mutate(aa = purrr::map_chr(aa, ~AMINOACIDS[[.x]][[2]])) |>
      dplyr::arrange(aa, dplyr::desc(val)) |>
      dplyr::mutate(aa = forcats::fct_reorder(aa, val, sum))
  #}

  plt <- ggplot2::ggplot(
    data = dat,
    mapping = ggplot2::aes(x = aa, y = val, label = codon)
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
    ggplot2::labs(x = 'Amino acid', y = 'Count') +
    ggthemes::theme_few()

  # if (scale == TRUE | scale %in% c('all', 'aa')) {
  #   plt <- plt +
  #     ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 100))
  # }

  return(plt)
})

#' Shows summary of a Cub object.
#' @param object A CubSet instance to summarize
#' @return NULL
#' @export
#' @examples
#' dat <- Cub('A DNA', PA14CDS[[1]])
#' summary(dat)
setMethod('summary', signature('CubSet'), function(object, ...) {
  lapply(object@data, summary)
})

setMethod('show', signature('CubSet'), function(object) {
  lapply(object@data, show)
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

  return(new('Counts', data = res))

})

setMethod('count', signature('CubSet'), function(object, each = FALSE) {

  if (each) {
    pb <- progress::progress_bar$new(
      format = " calculating CUBs genewise [:bar] :percent eta: :eta",
      total = length(object@data),
      clear = FALSE,
      width = 80
    )
    result <- lapply(object@data, \(gene) {
      pb$tick()
      count(gene)
    })
    return(result)
  } else {
    res <- unlist(lapply(object@data, function(datum) { datum@data }))
    res <- new('Cub', data = res)
    res <- count(res)

    return(res)
  }
})

setMethod('rscu', signature('CubSet'), function(object, each = FALSE) {

  if (each) {
    return(lapply(object@data, rscu))
  } else {
    res <- unlist(lapply(object@data, function(datum) { datum@data }))
    res <- new('Cub', data = res)
    res <- rscu(res)

    return(res)
  }
})

setMethod('rac', signature('CubSet'), function(object, ref, each = FALSE) {

  if (each) {
    return(lapply(object@data, rac, ref = ref))
  } else {
    res <- unlist(lapply(object@data, function(datum) { datum@data }))
    res <- new('Cub', data = res)
    res <- rac(res, ref = ref)

    return(res)
  }
})

setMethod('cai', signature('CubSet'), function(x, ref) {

  return(sapply(x@data, \(datum) { cai(x = datum, ref = ref) } ))
})

#' Plots RSCU as a stucked bar chart, represented RSCU(i, j) for each codon, grouped by
#' amino acid.
#' @param x A RSCU object to be plotted
#' @return A ggplot2::ggplot object
#' @export
#' @examples
#' Cub('A DNA', PA14CDS[[1]]) |>
#'   rscu() |>
#'   plot()
setMethod('plot', signature('Cub'), function(x, value = 'count', scale = FALSE, ...) {

  dat <- switch(
    value,
    count = count(x),
    rscu  = rscu(x),
    rac   = rac(x),
    #cai   = cai(object),
  )

  dat <- as.data.frame(dat) |>
    dplyr::filter(val > 0, aa != '*')

  if (scale == 'all' | scale == TRUE) {
    dat <- dat |>
      dplyr::mutate(val = val / sum(val)) |>
      dplyr::arrange(aa, dplyr::desc(val)) |>
      dplyr::mutate(aa = forcats::fct_reorder(aa, val, sum))
  } else if (scale == 'aa') {
    dat <- dat |>
      dplyr::group_by(aa) |>
      dplyr::mutate(val = val / sum(val), maxval = max(val)) |>
      dplyr::ungroup() |>
      dplyr::mutate(aa = forcats::fct_reorder(aa, maxval), maxval = NULL) |>
      dplyr::arrange(aa, desc(val))
  } else {
    dat <- dat |>
      dplyr::arrange(aa, dplyr::desc(val)) |>
      dplyr::mutate(aa = forcats::fct_reorder(aa, val, sum))
  }

  plt <- ggplot2::ggplot(
    data = dat,
    mapping = ggplot2::aes(x = aa, y = val, label = codon)
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
    ggplot2::labs(x = 'Amino acid', y = toupper(value)) +
    ggthemes::theme_few()

  if (scale == TRUE | scale %in% c('all', 'aa')) {
    plt <- plt +
      ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 100))
  }

  return(plt)
})

#' Plots RSCU as a stucked bar chart, represented RSCU(i, j) for each codon, grouped by
#' amino acid.
#' @param x A RSCU object to be plotted
#' @return A ggplot2::ggplot object
#' @export
#' @examples
#' Cub('A DNA', PA14CDS[[1]]) |>
#'   rscu() |>
#'   plot()
setMethod('plot', signature('CubSet'), function(x, value = 'count', scale = FALSE, ...) {

  dat <- unlist(lapply(x@data, \(datum) { datum@data }))
  dat <- new('Cub', data = dat)

  plot(dat, value, scale)
})
