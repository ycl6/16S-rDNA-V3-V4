# 16S rDNA V3-V4 amplicon sequencing analysis

This GitHub repository includes codes and scripts that demonstrate the use of `dada2` and `phyloseq` (and associated tools and R packages) to analyze 16S rDNA amplicon sequencing data. An working example is included in the `example` folder.

**Disclaimer:**
Do not use any of the provided codes and scripts in production without fully understanding of the contents. I am not responsible for errors or omissions or for any consequences from use of the contents and make no warranty with respect to the currency, completeness, or accuracy of the contents from this GitHub repository.

## Install required R packages

```
install.packages(c("ggplot2", "data.table", "plyr", "dplyr", "phangorn", "ggbeeswarm", "ggrepel", "vegan", "GUniFrac"))
```

For R version >= 3.5 (Bioconductor version >= 3.8):
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("dada2", "phyloseq", "DECIPHER", "DESeq2", "ShortRead", "Biostrings"))
```

For older versions of R:
```
source("https://bioconductor.org/biocLite.R")
biocLite(c("dada2", "phyloseq", "DECIPHER", "DESeq2", "ShortRead", "Biostrings"))
```

## Workflow

> **Update**: Check out the new script **cutadapt-and-dada2-per-run-processing.R**, which combines running of cutadapt and DADA2 per-run-processing in one script

### 1. PCR primer trimming using [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
* **run_trimming.pl**
* Cutadapt can run on multiple CPU cores in parallel but is only available on Python 3
* To enable it, use the option `-j N` (or the long form `--cores=N`), where N is the number of cores to use

### 2. [DADA2](https://benjjneb.github.io/dada2/) amplicon sequence variant (ASV) inference
* Part 1 - **dada2-per-run-processing.R**
* Part 2 - **dada2-create-phyloseq-obj.R**

### 3. Phylogenetic tree construction (RAxML)
[**RAxML**](https://github.com/stamatak/standard-RAxML):

Download the source file `standard-RAxML-x.x.x.tar.gz`

Use `tar zxvf standard-RAxML-x.x.x.tar.gz` to extract file

Use `make -f Makefile.SSE3.PTHREADS.gcc` to compile specific version

[**raxml-ng**](https://github.com/amkozlov/raxml-ng):

Download the pre-compiled binary `raxml-ng_vx.x.x_linux_x86_64.zip` (change x.x.x to downloaded version)

To extract files to a directory:
```
mkdir raxml-ng_vx.x.x
unzip raxml-ng_vx.x.x_linux_x86_64.zip -d raxml-ng_vx.x.x
```

### 4. Explore microbiome profiles using [phyloseq](https://joey711.github.io/phyloseq/) 
* **phyloseq-analysis.R**

### 5. Identify and plot differential taxa features using LEfSe and GraPhlAn
* **lefse-analysis.R**
* Requires ***Python 2.7***
* [**LEfSe**](https://bitbucket.org/nsegata/lefse/downloads/): Download and unzip `nsegata-lefse-9adc3a62460e.zip`
  * Additional packages required: rpy2, numpy, matplotlib, argparse
  * Additional R libraries required: survival, mvtnorm, modeltools, coin, MASS
* [**export2graphlan**](https://github.com/segatalab/export2graphlan): Install using `conda install -c bioconda export2graphlan`
  * Additional packages required: pandas, scipy
* [**GraPhlAn**](https://bitbucket.org/nsegata/graphlan/wiki/Home) Install using `conda install -c bioconda graphlan`
  * Additional packages required: biopython, matplotlib

## Download Silva and NCBI taxonomy DB

1. Silva taxonomic training data formatted for DADA2 (Silva version 138)
https://zenodo.org/record/3731176#.XvdOiihKhaQ
* silva_nr_v138_train_set.fa.gz
* silva_species_assignment_v138.fa.gz

2. Create from NCBI 16SMicrobial blast DB

Install BLAST 
```
conda install -c bioconda blast
```
or
```
sudo apt-get install ncbi-blast+
```

Download NCBI's 16S rRNA BLAST DB
```
wget ftp://ftp.ncbi.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz
tar zxf 16S_ribosomal_RNA.tar.gz
```

Convert 16SMicrobial BLAST DB into FASTA format
```
blastdbcmd -db 16S_ribosomal_RNA -entry all -out 16SMicrobial.fa
gzip 16SMicrobial.fa
```

## Download taxa_summary.R
http://evomics.org/phyloseq/taxa_summary-r/

Use `gunzip taxa_summary.R.gz` to extract file
