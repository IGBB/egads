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

### Restriction Enzymes
The restriction enzymes are downloaded from REBASE (
http://rebase.neb.com/rebase/rebase.html ) when compiling `egads`. A complete
list of enzyme names can be found at http://rebase.neb.com/cgi-bin/acrolist;
however, it may be more helpful to look through the various subsets available at
http://rebase.neb.com/rebase/rebase.charts.html .

## Testing

To recreate the report used in the "Innovations in ddRAD-seq" paper (
https://doi.org/10.1016/j.ab.2022.115001 ), download [the *G. longicalyx*
genome](https://www.cottongen.org/cottongen_downloads/Gossypium_longicalyx/F1_NSF/assembly/longicalyx.fasta.gz)
and run the following (will take a few minutes):

``` sh
./egads --rare NsiI,HindII,PstI \
        --freq Sau3AI,BfaI,MspI,HinP1I \
        --genome longicalyx.fasta.gz \
        --html longicalyx.egads.html
```

In the example, I downloaded the genome in the same directory as the `egads`
executable, and will output the report `longicalyx.egads.html` in that
directory; however, any relative or absolute path should work.
