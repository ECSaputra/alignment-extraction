# Alignment Extraction

This repository stores scripts for extracting alignments of coding regions and non-coding elements.  Installation of [RPHAST](https://github.com/CshlSiepelLab/RPHAST) and the `maf_parse` binary file from [PHAST](http://compgen.cshl.edu/phast/) are needed.


## Extraction of *coding region* alignments

Alignment extraction is done in a **chromosome-specific** manner. The starting multiple sequence alignment must be chromosome-speficic and in a MAF format. If alignment file is not chromosome-specific, split the MAF file to chromosome-specific MAFs (example: use [MafFilter](https://jydu.github.io/maffilter/)) *BEFORE* proceeding to the steps below.

The input datasets for this procedure are:

* GTF file for *each* chromosome of interest

* MAF of the chromosome of interest, fragmented into smaller chunks

Example input and output files can be found in folder `data` and `output`.

### Step 1: Split chromosome-specific MAF into smaller fragments

To make the operation computationally feasible, start by splitting the MAF of a chromosome into smaller fragments of about 100kb-300kb, depending on the size of your MAF. Use the `maf_parse` function from PHAST to do this (see documentation [here](http://compgen.cshl.edu/phast/help-pages/maf_parse.txt)).

Below is an example command for splitting the MAF of chrY (`mouse24way_chrY.maf`) into fragments of 300kb. The resulting outputs can be seen in the folder `data/alignment/chrY/`.

```
maf_parse mouse24way_chrY.maf --split 300000 --out-root data/alignment/chrY/chrYsplit
```

### Step 2: Obtain the coordinates of coding regions (CDS) from GTF file

This step produces a GTF file that only contains coding region coordinates of genes in the chromosome of interest. This step is performed by the function `getCDSCoordinates.R`, which takes the following arguments:

* `-i, --inputpath`: path to input GTF file

* `-o, --outputpath`: path to output GTF file (CDS only)

The following is an example of how to run the function from the command line:

```
Rscript getCDSCoordinates.R -i data/genepred/mm10.ncbiRefSeq.chrY.gtf -o output/CDS-coordinates/mm10.ncbiRefSeq.coding.chrY.gtf
```


### Step 3: Get gene boundaries and CDS coordinates per gene

This step produces data structures necessary for the subsequent extraction of coding region alignments. This step is performed by the function `getGeneBoundaries.R`, which takes the following arguments:

* `-c, --chromosome`: chromosome of interest

* `-i, --genepred_cds_path`: GTF file path containing the CDS coordinates in the chromosome of interest (output from Step 2)

* `-o, --output_folder`: output folder path (end path with slash)

The following is an example of how to run the function from the command line:

```
Rscript getGeneBoundaries.R -c chrY -i output/CDS-coordinates/mm10.ncbiRefSeq.coding.chrY.gtf -o output/CDS-information/
```

### Step 4: Extract coding region alignments

This step extracts the alignments of coding regions of interest in the FASTA format. This step is performed by the function `getCodingAlignments.R`, which takes the following arguments:

* `-c, --chromosome`: chromosome of interest

* `-r, --refseq`: reference sequence of alignment

* `-i, --cds_info_folder`: CDS information folder (output folder path from Step 3), end path with slash

* `-o, --output_folder`: output folder for coding region alignments

* `-a, --mafFolderPath`: path to folder containing chromosome MAF fragments (end path with slash)

* `-p, --prefix`: prefix of MAF fragments

The following is an example of how to run the function from the command line. Resulting alignments are reverse-complemented when necessary. Example outputs can be seen in the folder `output/coding-region-alignment/chrY/`.

```
Rscript getCodingAlignments.R -c chrY -r mm10 -i output/CDS-information/ -o output/coding-region-alignment/chrY/ -a data/alignment/chrY/ -p chrYsplit
```

## Extraction of *non-coding region* alignments

Similar to the coding region extraction, non-coding region alignment extraction is also done in a **chromosome-specific** manner.

The input datasets for this procedure are:

* BED file containing coordinates of non-coding regions for *each* chromosome of interest

* MAF of the chromosome of interest, fragmented into smaller chunks

### Step 1: Split chromosome-specific MAF into smaller fragments

The same as Step 1 for coding regions.

### Step 2: Extract non-coding region alignments

This step extracts the alignments of non-coding regions of interest in the FASTA format. This step is performed by the function `getNonCodingAlignments.R`, which takes the following arguments:

* `-c, --chromosome`: chromosome of interest

* `-r, --refseq`: reference sequence of alignment

* `-i, --elementBedPath`: path to BED file containing coordinates of regions of interest

* `-o, --output_folder`: output folder for non-coding region alignments

* `-a, --mafFolderPath`: path to folder containing chromosome MAF fragments (end path with slash)

* `-p, --prefix`: prefix of MAF fragments

The following is an example of how to run the function from the command line. Example outputs can be seen in the folder `output/CNE-alignment/`.

```
Rscript getNonCodingAlignments.R -c chrY -r mm10 -i data/CNE-coordinates/mouseCNEs.chrY.bed -o output/CNE-alignment/ -a data/alignment/chrY/ -p chrYsplit
```