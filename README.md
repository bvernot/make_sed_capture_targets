# make_sed_capture_targets

This Python script queries an ENSEMBL database, and infers mammalian diversity in the region around 

## Input data

The script takes three files:
1. An alleles file giving the chromosome, position, and both expected alleles.
1. A bed file with information on repeats in the region. This file is not actually used in this step but the last three columns are:
⋅⋅* Number of bases in a window (in this case 104 bases) that are covered by simple repeats.
..* Number of bases in a window (in this case 104 bases) that are _unmasked_ by the Heng Li 35mer mask.
..* Is the target base masked by _unmasked_ by the Heng Li 35mer mask? (0 or 1)
1. Sites from the "ape" file, as defined here: https://dx.doi.org/10.17617/3.5h

```
    $ head snps.alleles
    1       752675  T       C
    1       812425  G       A
    1       812751  T       C
    1       813034  A       G
    1       814609  A       T
    1       821477  A       G
    1       821947  G       A
    1       822775  G       A
    1       826240  C       A
    1       834198  T       C
    
    $ head snps.rpt_map
    1       752674  752675  0       54      1
    1       812424  812425  0       33      1
    1       812750  812751  0       12      1
    1       813033  813034  0       12      1
    1       814608  814609  0       2       1
    1       821476  821477  0       26      1
    1       821946  821947  0       35      1
    1       822774  822775  0       33      1
    1       834197  834198  0       104     1
    1       834359  834360  0       67      1

    $ head snps.ape
    1       812425  G       G       G       A       A       A       G       R       G
    1       812751  T       C       C       C       C       C       T       T       C
    1       813034  A       G       G       G       G       G       R       G       A
    1       814609  A       A       A       A       T       T       A       W       W
    1       821477  A       G       G       A       A       N       A       A       G
    1       821947  G       G       G       G       G       N       R       G       G
    1       822775  G       G       G       G       G       G       A       A       G
    1       826240  C       T       T       T       T       T       C       M       M
    1       834198  T       T       T       N       T       T       C       C       T
    1       835831  G       G       G       N       G       G       G       G       A


    python3.5 query_ensembl_for_sites.py 25 52 \
                          snps.alleles snps.rpt_map snps.ape \
                          > snps.txt
```                          

This file contains a lot of debugging information - but to get a table of capture sites, grep for lines that include REPORT.

