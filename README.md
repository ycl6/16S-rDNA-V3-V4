# 16S rDNA V3-V4 amplicon sequencing analysis

1. PCR primer trimming (run_trimming.pl)
2. dada2 part 1 (dada2-per-run-processing.R)
3. dada2 part 2 (dada2-create-phyloseq-obj.R)
4. phylogenetic tree construction (RAxML)
* https://github.com/stamatak/standard-RAxML
* https://github.com/amkozlov/raxml-ng
5. phyloseq
6. LEfSe

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
