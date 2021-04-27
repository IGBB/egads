# egads 
Electronically Guided Digestion Selection - A ddRAD-seq enzyme selection aid


## Installation

``` sh
git clone --recursive https://github.com/IGBB/egads
make -C egads/src
```

## Running

`egads` takes three positional arguments: 1) a comma-delimited list of rare
cutting restriction enzymes, 2) a comma-delimited list of frequent cutting
restriction enzymes, and 3) a fasta file, optionally compressed with gzip. It
outputs an html report (see example.html) to `stdout`.

``` sh
./egads <rare[,...]> <freq[,...]> <fasta[.gz]> > <output>
```

