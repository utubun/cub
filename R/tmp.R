library(tidyverse)

cfreq <- read_csv('./dat/ftb/codons_rel_freq.csv')

# First glance

cfreq %<>%
  mutate_at(vars(anticodon, codon), ~ gsub('^(.+)-\\d+\\w*\\d*$', '\\1', .)) %>% # 1O1 with O instead of 0
  tidyr::separate(
    col  = anticodon,
    into = c('aa', 'anticodon'),
    sep  = '-'
  ) %>%
  dplyr::mutate(codon = gsub('^\\w+-(.+)$', '\\1', codon)) %>% # aa abbreviation
  dplyr::mutate(
    codon   = purrr::map(
      codon,
      ~ paste0(
        substr(., 1, 2),
        unlist(strsplit(substr(., 3, nchar(.)), split = '/'))
      )
    )
  ) %>%
  tidyr::unnest(cols = codon) %>%
  dplyr::group_by(codon) %>%
  dplyr::summarise(freq = mean(freq)) %>%
  dplyr::mutate(aa = as.character(translate(DNAStringSet(codon)))) %>%
  dplyr::select(aa, codon, freq) %>%
  bind_rows(
    data.frame(
      aa    = GENETIC_CODE[!(names(GENETIC_CODE) %in% .$codon)],
      codon = names(GENETIC_CODE[!(names(GENETIC_CODE) %in% .$codon)]),
      freq  = 1,
      stringsAsFactors = FALSE
    )
  ) %>%
  dplyr::filter(aa != '*')

cdn_freq_lookup <- as.list(cfreq$freq)

names(cdn_freq_lookup) <- cfreq$codon

readr::write_rds(cdn_freq_lookup, './dat/ftb/cdn_freq_lookup.Rds')

rm(list = c('cfreq', 'cdn_freq_lookup'))
