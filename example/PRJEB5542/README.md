# Study: ERP004961 (PRJEB5542) 
### 16S metabarcoding of bacteria associated with cultured strains of the brown alga Ectocarpus sp.
Paper: https://www.nature.com/articles/ismej2015104
EBI ENA: https://www.ebi.ac.uk/ena/data/view/PRJEB5542

We will use this dataset as an example to perform 16S rDNA V3-V4 amplicon sequencing analysis. This dataset contain brown alga (*Ectocarpus sp.*) associated samples gather from different geographic locations. The V3 and V4 regions of the 16S rRNA gene fragment was amplified and sequenced on the Illumina MiSeq platform using V2 reagents (2x250 bp reads).

### Illumina Sequencing Primers

**341F**: CTATGGTAATTCGTCCTACGGRAGGCAGCAG

**806R**: AGTCAGTCAGCCGGACTACCRGGGTATCTAAT

## dada2 filterAndTrim() parameter for this study
```
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c(240,230), minLen=200, maxN=0, truncQ=2, maxEE=c(2,5),
                     rm.phix=TRUE, compress=TRUE, verbose=TRUE, multithread=TRUE)
```
