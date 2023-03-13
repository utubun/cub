## code to prepare `Sharp166` dataset goes here

Sharp166 <- data.frame(
  gene = c(
    'rpsU', 'rpsJ', 'rpsT', 'rpsL', 'rpsA', 'rpsB', 'rpsO', 'rpsG', 'rpmB',
    'rpmG', 'rpmH', 'rplK', 'rplJ', 'rplA', 'rplL', 'rplQ', 'rplC', 'lpp',
    'ompA', 'ompC', 'ompF', 'tufA', 'tufB', 'tsf',  'fusA', 'recA', 'dnaK',
    'rpoA', 'rpoB', 'rpoC', 'rpoD', 'uncA', 'uncD', 'uncE', 'alaS', 'metG',
    'glnS', 'glyS', 'gly2', 'thrS', 'trpS', 'tyrS', 'dnaG', 'araC', 'lacI',
    'trpR', 'cytR', 'deoR', 'galR', 'lexA', 'aceE', 'aceF', 'asnA', 'carB',
    'crp',  'cya',  'deoC', 'dnaA', 'dnaB', 'dnaN', 'dxi',  'fimA', 'frdA',
    'frdB', 'ftsA', 'fumA', 'gdhA', 'glgC', 'gltA', 'glyA', 'gnd',  'gpt',
    'gshw', 'himA', 'hsdS', 'ilvG', 'infB', 'infC', 'lamB', 'lep',  'lpd',
    'malE', 'malK', 'metK', 'mtlA', 'ndh',  'nrdA', 'nrdB', 'musA', 'papA',
    'phoS', 'pldA', 'plsB', 'polA', 'proA', 'rbsP', 'rho',  'sdhA', 'ssb',
    'sucA', 'sucB', 'thyA', 'tolC', 'trxA', 'uncB', 'uncF', 'uncG', 'alkA',
    'argF', 'argI', 'aroA', 'aroF', 'aroG', 'asd',  'carA', 'dam',  'dnaQ',
    'fol',  'ftsQ', 'ilvM', 'kdpA', 'kdpB', 'kdpC', 'lacA', 'lacY', 'lacZ',
    'lspA', 'lysA', 'malF', 'malB', 'metB', 'metL', 'motA', 'motB', 'pabA',
    'pabB', 'pbpB', 'pfkB', 'phoE', 'phr',  'proB', 'proC', 'purF', 'pyrB',
    'recF', 'rnh',  'sdhB', 'sulA', 'tap',  'tar',  'thrA', 'thrB', 'tnaA',
    'tonB', 'trg',  'trpA', 'trpB', 'trpC', 'trpD', 'trpE', 'tar',  'uncC',
    'uncH', 'uvrD', 'xylA'
  ),
  expression = factor(
    c(rep('very high', 27), rep('high', 15), rep('regulatory', 8), rep('other', 115)),
    levels = c('very high', 'high', 'regulatory', 'other')
  )
)

#org <- KEGGREST::keggList('organism')

usethis::use_data(Sharp166, overwrite = TRUE, compress = 'xz')
