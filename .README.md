# Code book

## Project directory structure

- `src` directory containes scripts for the data analysis and preprocessing;
  - `src/r` R scripts, related to the preprocessing and final analyses;
  - `src/sh` bash scripts for the data preprocessing;
- `dat` containe data related to the project:
  - `raw` raw, read-only data;
  - `int` intermediate data:
    - `gnm` genome assembies for use in R scripts;
    - `seq` sequences;
    - `trs` output of **tr**na**s**can program;
- `README.md` this file. 

## tRNAscan-SE

```bash
tRNAscan-SE -B -L -H -D -o /home/dst20/dev/ribotracker/dat/trnascan/${FNAME}.out -f /home/dst20/dev/ribotracker/dat/trnascan/${FNAME}.struct  -a /home/dst20/dev/ribotracker/dat/tmp/${FNAME}.fa -m /home/dst20/dev/ribotracker/dat/tmp/${FNAME}.stats -l /home/dst20/dev/ribotracker/dat/tmp/${FNAME}.log -d -Q $f;
```

### Options

- `-B` search for bacterial tRNAs;
- `-L` search using the legacy method (tRNAscan, EufindtRNA, COVE), use with `-E`, `-B`, `-A` and `-G` only;
- `-H` show breakdown of primary and secondary structure components to covariance model bit scores;
- `-D` disable pseudogene checking;
- `-o` save final results in `<file>`;
- `-f` save tRNA secondary structures to `<file>`;
- `-a` save predicted tRNA sequences in FASTA file format of `<file>`;
- `-m` save statistics summary for run in `<file>`;
- `-b` save result in BED file format of `<file>`;
- `-l` save log of program progress in `<file>`;
- `-d` show the progress;
- `-Q` do not prompt before overwriting pre-existing result files;