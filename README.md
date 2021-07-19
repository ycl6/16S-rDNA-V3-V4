# 16S rDNA V3-V4 amplicon sequencing analysis

This GitHub repository includes codes and scripts that demonstrate the use of `dada2` and `phyloseq` (and associated tools and R packages) to analyze 16S rDNA amplicon sequencing data. An working example is included in the `example` folder.

**Disclaimer:**
Do not use any of the provided codes and scripts in production without fully understanding of the contents. I am not responsible for errors or omissions or for any consequences from use of the contents and make no warranty with respect to the currency, completeness, or accuracy of the contents from this GitHub repository.

## Set up conda environment

I strongly recommend using [conda](https://conda.io/projects/conda/en/latest/index.html) to manage the software and R packages required for this analysis.

- We will create an environment with R >= 4.1 and Python 3 to perform majority of the analysis.
- We will create an environment with Python 2.7 to run `export2graphlan` and `graphlan` as they do not support Python 3.
- If you are unable to install `picrust2` >= 2.4 or higher in the *main* environment with Python 3 due to conflicts, you will need to create a separate environment to run the latest `picrust2`.
- Use `conda activate my_env` to activate an environment named `my_env`, and `conda deactivate` to deactivate an environment. Learn more about managing environments [here](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

```sh
# Create main "16S" env
# Remove `picrust2=2.4` from the list if it causes conflicts
conda create -n 16S -c conda-forge -c bioconda r-base=4.1 python=3 curl libcurl openssl boost-cpp \
r-doparallel r-devtools r-rcurl r-httr r-magick r-png r-ggplot2 r-data.table r-phangorn r-ape r-gridextra \
r-ggbeeswarm r-ggrepel r-vegan r-tidyverse r-gtools r-r.utils bioconductor-dada2 bioconductor-phyloseq \
bioconductor-decipher bioconductor-deseq2 bioconductor-shortread bioconductor-biostrings \
bioconductor-biomformat bioconductor-aldex2 cutadapt raxml raxml-ng lefse picrust2=2.4

# Use "conda activate 16S" to activate this environment
```

```sh
# Create "graphlan" env
conda create -n graphlan export2graphlan graphlan

# Use "conda activate graphlan" to activate this environment
```

```sh
# Create "picrust2" env if it causes conflicts when setting up the "16S" env
conda create -n picrust2 picrust2=2.4

# Use "conda activate picrust2" to activate this environment
```

## Workflow

### 1. PCR primer trimming using [cutadapt](https://cutadapt.readthedocs.io/en/stable/)

* Conda environment: **16S**
* Script: "**run_trimming.pl**"
* Cutadapt can run on multiple CPU cores in parallel under ***Python 3***
* See "2a - Option 2" to run `cutadapt` in R

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

### 2. Run [DADA2](https://benjjneb.github.io/dada2/) workflow

* Conda environment: **16S**

#### 2a.

* Script:
  * Option 1: "**dada2-per-run-processing.R**" - Trimming already performed using `cutadapt` in terminal/console
  * Option 2: "**cutadapt-and-dada2-per-run-processing.R**" - Execute `cutadapt` in R by using `system2` function
* dada2 filtering and trimming `filterAndTrim`
* Amplicon sequence variant (ASV) inference

#### 2b.

* Script: "**dada2-create-phyloseq-obj.R**"
* Taxonomy and species assignment
* Sample data preparation
* Phylogenetic tree construction using [**RAxML**](https://github.com/stamatak/standard-RAxML) and [**raxml-ng**](https://github.com/amkozlov/raxml-ng)
* Phyloseq object construction

### 3. Explore microbiome profiles using [phyloseq](https://joey711.github.io/phyloseq/)

* Conda environment: **16S**
* Script: "**phyloseq-analysis.R**"

### 4. Identify and plot differential taxa features using LEfSe and GraPhlAn

* Conda environment: **16S** and **graphlan**
* Script: "**lefse-analysis.R**"
* Linear discriminant analysis (LDA) effect size ([**LEfSe**](https://github.com/SegataLab/lefse)) analysis
* Create cladogram using [**export2graphlan**](https://github.com/segatalab/export2graphlan) and [**GraPhlAn**](https://github.com/biobakery/graphlan). ***Require Python 2.7***

### 5. Infer functional capacity using [picrust2](https://github.com/picrust/picrust2) and statistical analysis using [ALDEx](https://github.com/ggloor/ALDEx_bioc)

* Conda environment: **16S** or **picrust2** (depending on where `picrust2` is installed)
* Script: "**picrust2-aldex-analysis.R**"

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
