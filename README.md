# 16S rDNA V3-V4 amplicon sequencing analysis

## Install required R packages

```
install.packages(c("ggplot2","data.table","dplyr","phangorn","ggbeeswarm","ggrepel","vegan"))
```

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("dada2")
BiocManager::install("phyloseq")
BiocManager::install("DECIPHER")
```

## Workflow

1. PCR primer trimming (***run_trimming.pl***)
* cutadapt: https://cutadapt.readthedocs.io/en/stable/

2. DADA2 part 1 (***dada2-per-run-processing.R***) and part 2 (***dada2-create-phyloseq-obj.R***)
* DADA2 GitHub: https://benjjneb.github.io/dada2/

3. phylogenetic tree construction (RAxML)
* https://github.com/stamatak/standard-RAxML (make -f Makefile.SSE3.PTHREADS.gcc)
* https://github.com/amkozlov/raxml-ng (下載已編譯版本 raxml-ng_vx.x.x_linux_x86_64.zip)

4. phyloseq (***phyloseq-analysis.R***)
* phyloseq  GitHub: https://joey711.github.io/phyloseq/

5. LEfSe and GraPhlAn (***lefse-analysis.R***)
* LEfSe Tutorial: https://bitbucket.org/biobakery/biobakery/wiki/lefse
* GraPhlAn Tutorial: https://bitbucket.org/nsegata/graphlan/wiki/Home

## Download Silva and NCBI taxonomy DB

1. Silva taxonomic training data formatted for DADA2 (Silva version 132)
https://zenodo.org/record/1172783#.XPFLeUdS8-U
* silva_nr_v132_train_set.fa.gz
* silva_species_assignment_v132.fa.gz

2. Create from NCBI 16SMicrobial blast DB

Install BLAST 
```
sudo apt-get install ncbi-blast+
```

Download NCBI's 16SMicrobial BLAST DB
```
wget ftp://ftp.ncbi.nih.gov/blast/db/16SMicrobial.tar.gz
tar zxf 16SMicrobial.tar.gz
```

Convert 16SMicrobial BLAST DB into FASTA format
```
blastdbcmd -db 16SMicrobial -entry all -out 16SMicrobial.fa
gzip 16SMicrobial.fa
```

## Download taxa_summary.R
http://evomics.org/phyloseq/taxa_summary-r/
