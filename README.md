# egads 
Electronically Guided Digestion Selection - A ddRAD-seq enzyme selection aid


## Installation

``` sh
git clone --recursive https://github.com/IGBB/egads
make -C egads/src
```

## Running

`egads` requires three arguments: 1) `--rare`: a comma-delimited list of rare
cutting restriction enzymes, 2) `--freq`: a comma-delimited list of frequent
cutting restriction enzymes, and 3) `--genome`: a fasta file, optionally
compressed with gzip. The default output is an html report (see example.html) to
`stdout`. The arguments `--html`, html report, and `--bed`, fragments in bed
format, can be used to direct the output elsewhere.

``` sh
./egads --rare <rare[,...]> --freq <freq[,...]> --genome <fasta[.gz]> > <output>
```

