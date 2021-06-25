# 16S rDNA V3-V4 amplicon sequencing analysis

This GitHub repository includes codes and scripts that demonstrate the use of `dada2` and `phyloseq` (and associated tools and R packages) to analyze 16S rDNA amplicon sequencing data. An working example is included in the `example` folder.

**Disclaimer:**
Do not use any of the provided codes and scripts in production without fully understanding of the contents. I am not responsible for errors or omissions or for any consequences from use of the contents and make no warranty with respect to the currency, completeness, or accuracy of the contents from this GitHub repository.

## Install required R packages

```
install.packages(c("ggplot2", "data.table", "R.utils", "plyr", "dplyr", "phangorn", "ape", "reshape2", 
	"gridExtra", "ggbeeswarm", "ggrepel", "vegan", "GUniFrac", "tibble"))
```

For R version >= 3.5 (Bioconductor version >= 3.8):
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("dada2", "phyloseq", "DECIPHER", "DESeq2", "ShortRead", "Biostrings", 
	"biomformat", "ALDEx2"))
```

For older versions of R:
```
source("https://bioconductor.org/biocLite.R")
biocLite(c("dada2", "phyloseq", "DECIPHER", "DESeq2", "ShortRead", "Biostrings", "biomformat", "ALDEx2"))
```

## Workflow

> **Update**: Check out the new script **cutadapt-and-dada2-per-run-processing.R**, which combines running of cutadapt and DADA2 per-run-processing in one script

### 1. PCR primer trimming using [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
* **run_trimming.pl**
* Cutadapt can run on multiple CPU cores in parallel but is only available on ***Python 3***

```
Usage: perl run_trimming.pl project_folder fastq_folder forward_primer_sequence reverse_primer_sequence
```

* **Example:**
  - project_folder: PRJEB27564
  - fastq_folder: raw
  - forward_primer_sequence: 5'-CCTACGGGNGGCWGCAG-3'
  - reverse_primer_sequence: 5'-GACTACHVGGGTATCTAATCC-3'

```
Usage: perl run_trimming.pl PRJEB27564 raw CCTACGGGNGGCWGCAG GACTACHVGGGTATCTAATCC
```

### 2a. [DADA2](https://benjjneb.github.io/dada2/) Part 1
* **dada2-per-run-processing.R** - `cutadapt` performed before-hand
* **cutadapt-and-dada2-per-run-processing.R** - Execute `cutadapt` in R by using `system2` function
* dada2 filtering and trimming `filterAndTrim`
* Amplicon sequence variant (ASV) inference

### 2b. [DADA2](https://benjjneb.github.io/dada2/) Part 2
* **dada2-create-phyloseq-obj.R**
* Taxonomy and species assignment
* Sample data preparation
* Phylogenetic tree construction

  1. [**RAxML**](https://github.com/stamatak/standard-RAxML):

  Download the source file `standard-RAxML-x.x.x.tar.gz`

  Use `tar zxvf standard-RAxML-x.x.x.tar.gz` to extract file

  Use `make -f Makefile.SSE3.PTHREADS.gcc` to compile specific version

  2. [**raxml-ng**](https://github.com/amkozlov/raxml-ng):

  Download the pre-compiled binary `raxml-ng_vx.x.x_linux_x86_64.zip` (change x.x.x to downloaded version)

  To extract files to a directory:
  ```
  mkdir raxml-ng_vx.x.x
  unzip raxml-ng_vx.x.x_linux_x86_64.zip -d raxml-ng_vx.x.x
  ```

* Phyloseq object construction

### 3. Explore microbiome profiles using [phyloseq](https://joey711.github.io/phyloseq/) 
* **phyloseq-analysis.R**

### 4. Identify and plot differential taxa features using LEfSe and GraPhlAn
* **lefse-analysis.R**
* **LEfSe** and **GraPhlAn** requires ***Python 2.7***
* One can create a Python 2 conda environment to work with LEfSe and GraPhlAn

```
# Create environment
conda create --name p27 python=2.7

# Activate environment
source activate p27

# Deactivate environment
conda deactivate
```

* [**LEfSe**](https://bitbucket.org/nsegata/lefse/downloads/): Download and unzip `nsegata-lefse-9adc3a62460e.zip`
  * Additional packages required: rpy2, numpy, matplotlib, argparse
  * Additional R libraries required: survival, mvtnorm, modeltools, coin, MASS
* [**export2graphlan**](https://github.com/segatalab/export2graphlan): Install using `conda install -c bioconda export2graphlan`
  * Additional packages required: pandas, scipy
* [**GraPhlAn**](https://bitbucket.org/nsegata/graphlan/wiki/Home) Install using `conda install -c bioconda graphlan`
  * Additional packages required: biopython, matplotlib

### 5. Infer functional capacity using [picrust2](https://github.com/picrust/picrust2) and statistical analysis using [ALDEx](https://github.com/ggloor/ALDEx_bioc)
* **picrust2-aldex-analysis.R**
* **picrust2** requires ***Python 3.5/3.6***
* Use `conda` to create conda environment with Python 3.6 and picrust2 2.3.0_b installed

```
# Create environment
conda create -n picrust2 -c bioconda -c conda-forge picrust2=2.3.0_b python=3.6

# Activate environment
conda activate picrust2

# Deactivate environment
conda deactivate
```

## Download Silva and NCBI taxonomy DB

1. Silva taxonomic training data formatted for DADA2 (Silva version 138.1)
https://zenodo.org/record/4587955#.YNWax3VKhkY

> Note: These files are intended for use in classifying prokaryotic 16S sequencing data and are not appropriate for classifying eukaryotic ASVs.

* silva_nr99_v138.1_train_set.fa.gz
* silva_species_assignment_v138.1.fa.gz

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
