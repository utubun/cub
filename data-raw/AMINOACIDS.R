## code to prepare `AMINOACIDS` dataset goes here
AMINOACIDS <- list(
  A = c('Ala', 'Alanine', 'NC(C)C(O)=O'),
  R = c('Arg', 'Arginine', 'NC(CCCNC(N)=N)C(O)=O'),
  N = c('Asn', 'Asparagine', 'NC(CC(N)=O)C(O)=O'),
  D = c('Asp', 'Aspartic acid', 'NC(CC(O)=O)C(O)=O'),
  C = c('Cys', 'Cysteine', 'NC(CS)C(O)=O'),
  E = c('Glu', 'Glutamic acid', 'NC(CCC(O)=O)C(O)=O'),
  Q = c('Gln', 'Glutamine', 'NC(CCC(N)=O)C(O)=O'),
  G = c('Gly', 'Glycine', 'NC([H])C(O)=O'),
  H = c('His', 'Histidine', 'NC(CC1CCNCN1)C(O)=O'),
  I = c('Ile', 'Isoleucine', 'NC(C(CC)C)C(O)=O'),
  L = c('Leu', 'Leucine', 'NC(CC(C)C)C(O)=O'),
  K = c('Lys', 'Lysine', 'NC(CCCCN)C(O)=O'),
  M = c('Met', 'Methionine', 'NC(CCSC)C(O)=O'),
  F = c('Phe', 'Phenylalanine', 'NC(CC1CCCCC1)C(O)=O'), # toupper must be upplied
  P = c('Pro', 'Proline', 'OC(C1CCCN1)=O'),
  S = c('Ser', 'Serine', 'NC(CO)C(O)=O'),
  T = c('Thr', 'Threonine', 'NC(C(C)O)C(O)=O'),
  W = c('Trp', 'Tryptophan', 'NC(CC1CNC2C1CCCC2)C(O)=O'),
  Y = c('Tyr', 'Tyrosine', 'NC(CC1CCC(O)CC1)C(O)=O'),
  V = c('Val', 'Valine', 'NC(C(C)C)C(O)=O'),
  U = c('Sec', 'Selenocysteine', 'C(C(C(=O)O)N)[Se]'),
  O = c('Pyl', 'Pyrrolysine', 'CC1CC=NC1C(=O)NCCCCC(C(=O)O)N'),
  HCY = c('Hcy','Homocysteine', 'C(CS)C(C(=O)O)N')
)

usethis::use_data(AMINOACIDS, overwrite = TRUE)
