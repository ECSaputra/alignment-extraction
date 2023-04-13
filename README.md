# Alignment Extraction

This repository stores scripts for extracting alignments of coding regions and non-coding elements.  Installation of [RPHAST](https://github.com/CshlSiepelLab/RPHAST) and the `maf_parse` binary file from [PHAST](http://compgen.cshl.edu/phast/) are needed.


## Extraction of coding region alignments

Alignment extraction is done in a chromosome-specific manner. The starting multiple sequence alignment must be chromosome-speficic and in a MAF format. If alignment file is not chromosome-specific, split the MAF file to chromosome-specific MAFs (example: use [MafFilter](https://jydu.github.io/maffilter/)).


