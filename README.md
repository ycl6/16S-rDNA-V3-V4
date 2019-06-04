# 16S rDNA V3-V4 amplicon sequencing analysis

1. PCR primer trimming (***run_trimming.pl***)
* cutadapt: https://cutadapt.readthedocs.io/en/stable/

2. DADA2 part 1 (***dada2-per-run-processing.R***) and part 2 (***dada2-create-phyloseq-obj.R***)
* DADA2 GitHub: https://benjjneb.github.io/dada2/

3. phylogenetic tree construction (RAxML)
* https://github.com/stamatak/standard-RAxML
* https://github.com/amkozlov/raxml-ng

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
```
wget ftp://ftp.ncbi.nih.gov/blast/db/16SMicrobial.tar.gz
tar zxf 16SMicrobial.tar.gz

blastdbcmd -db 16SMicrobial -entry all -out 16SMicrobial.fa
gzip 16SMicrobial.fa
```

## Download taxa_summary.R
http://evomics.org/phyloseq/taxa_summary-r/
