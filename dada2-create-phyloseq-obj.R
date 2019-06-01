#

library("dada2")
library("data.table")
library("DECIPHER")
library("phangorn")
library("phyloseq")
library("ggplot2")
options(width=190)

# Load sequence table from multiple runs
s01 <- readRDS("seqtab1.rds")
s02 <- readRDS("seqtab2.rds")

# Merge sequence table and remove chimeras
st.all <- mergeSequenceTables(s01, s02)
st.nochim <- removeBimeraDenovo(st.all, verbose=TRUE, multithread=TRUE)

SVformat = paste("%0",nchar(as.character(ncol(st.nochim))),"d", sep="")
svid <- paste0("SV", sprintf(SVformat, seq(ncol(st.nochim))))

# Assign taxonomy
ref1 <- "/path-to-db/silva_nr_v132_train_set.dada2.fa.gz"
ref2 <- "/path-to-db/silva_species_assignment_v132.dada2.fa.gz"
ref3 <- "/path-to-db/16SMicrobial.fa.gz"

taxtab <- assignTaxonomy(st.nochim, refFasta=ref1, minBoot=80, tryRC = TRUE, outputBootstraps=TRUE, verbose=TRUE, multithread=TRUE)
spec_silva <- assignSpecies(getSequences(st.nochim), ref2, allowMultiple = FALSE, tryRC = TRUE, verbose = TRUE)
spec_ncbi <- assignSpecies(getSequences(st.nochim), ref3, allowMultiple = FALSE, tryRC = TRUE, verbose = TRUE)

s_silva = data.frame(spec_silva)
rownames(s_silva) = svid

s_ncbi = data.frame(spec_ncbi)
rownames(s_ncbi) = svid
s_ncbi[grep("\\[",s_ncbi$Genus),]$Species = NA
s_ncbi[grep("\\[",s_ncbi$Genus),]$Genus = NA

s_merged = cbind(s_ncbi, s_silva)
colnames(s_merged) = c("nGenus","nSpecies","sGenus","sSpecies")
s_merged1 = s_merged[!is.na(s_merged$nSpecies),]	                        # NCBI assignment not empty
colnames(s_merged1)[1:2] = c("Genus","Species")
s_merged2 = s_merged[is.na(s_merged$nSpecies) & !is.na(s_merged$sSpecies),]	# no NCBI assignment but Silva not empty
colnames(s_merged2)[3:4] = c("Genus","Species")
s_merged3 = s_merged[is.na(s_merged$nSpecies) & is.na(s_merged$sSpecies),]	# no NCBI and Silva assignment
colnames(s_merged3)[3:4] = c("Genus","Species")

s_final = rbind(s_merged1[,c("Genus","Species")], s_merged2[,c("Genus","Species")], s_merged3[,c("Genus","Species")])
s_final = s_final[order(row.names(s_final)),]

s_final = as.matrix(s_final)

if("Genus" %in% colnames(taxtab$tax)) {
        gcol <- which(colnames(taxtab$tax) == "Genus")
} else { gcol <- ncol(taxtab$tax) }

matchGenera <- function(gen.tax, gen.binom, split.glyph="/") {
  if(is.na(gen.tax) || is.na(gen.binom)) { return(FALSE) }
  if((gen.tax==gen.binom) ||
     grepl(paste0("^", gen.binom, "[ _", split.glyph, "]"), gen.tax) ||
     grepl(paste0(split.glyph, gen.binom, "$"), gen.tax)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

gen.match <- mapply(matchGenera, taxtab$tax[,gcol], s_final[,1])
taxtab$tax <- cbind(taxtab$tax, s_final[,2])
colnames(taxtab$tax)[ncol(taxtab$tax)] <- "Species"
#print(paste(sum(!is.na(s_final[,2])), "out of", nrow(s_final), "were assigned to the species level."))

taxtab$tax[!gen.match,"Species"] <- NA
#print(paste("Of which", sum(!is.na(taxtab$tax[,"Species"])),"had genera consistent with the input table."))

# Prepare df
df <- data.frame(sequence=colnames(st.nochim), abundance=colSums(st.nochim))
SVformat = paste("%0",nchar(as.character(nrow(df))),"d", sep="")
df$id <- paste0("SV", sprintf(SVformat, seq(nrow(df))))
#uniquesToFasta(df, "dada2_out.fasta", id=df$id)

df = merge(df, as.data.frame(taxtab), by = "row.names")
rownames(df) = df$id
df = df[order(df$id),2:ncol(df)]

# Construct phylogenetic tree
seqs <- getSequences(st.nochim)
names(seqs) <- df$id # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

# Export alignment
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
write.phyDat(phang.align, file="alignment.fasta", format="fasta")
write.phyDat(phang.align, file="alignment.aln", format="phylip")

#######################################################
# Run below 2 commands outside R in Linux environment
#raxmlHPC-PTHREADS-SSE3 -T 2 -f E -p 1234 -x 5678 -m GTRCAT -N 1 -s alignment.aln -n raxml_tree_GTRCAT
#raxml-ng --evaluate --force --seed 1234 --log progress --threads 2 --msa alignment.fasta --model GTR+G --tree RAxML_fastTree.raxml_tree_GTRCAT --brlen scaled --prefix GTRCAT
#######################################################

# Import tree
raxml_tree <- read_tree("GTRCAT.raxml.mlTrees")

# Load sample info (sample.meta); Example content:
#  Sample_ID   File_ID    Batch     Group
#    Sample1  oldname1     run1         A
#    Sample2  oldname2     run1         B
#    Sample3  oldname3     run2         A
#    Sample4  oldname4     run2         B

samdf <- data.frame(fread("sample.meta", colClasses = "character"))

rownames(samdf) = samdf$Sample_ID
samdf$Sample_ID = as.factor(samdf$Sample_ID)
rownames(st.nochim) = as.character(samdf[match(rownames(st.nochim), samdf$File_ID),]$Sample_ID) # Update id if necessary
samdf$Sample_ID = factor(samdf$Sample_ID, levels=c(sort(levels(samdf$Sample_ID), decreasing=F)))
samdf$Batch = as.factor(samdf$Batch)
samdf$Group = as.factor(samdf$Group)

new_seqtab = st.nochim
new_taxtab = taxtab
colnames(new_seqtab) = df[match(colnames(new_seqtab), df$sequence),]$id
rownames(new_taxtab$tax) = df[match(rownames(new_taxtab$tax), df$sequence),]$id
tax = as.data.frame(new_taxtab$tax)
tax$Family = as.character(tax$Family)
tax$Genus = as.character(tax$Genus)
tax[grep("Family",tax$Family),]$Family = paste0(tax[grep("Family",tax$Family),]$Class,"_",tax[grep("Family",tax$Family),]$Family)
tax[is.na(tax$Genus) & !is.na(tax$Family),]$Genus = paste(tax[is.na(tax$Genus) & !is.na(tax$Family),]$Family,"_ge",sep="")

# Combine data into a phyloseq object
ps <- phyloseq(tax_table(as.matrix(tax)), sample_data(samdf), otu_table(new_seqtab, taxa_are_rows = FALSE), phy_tree(raxml_tree))

# Save current workspace
# save.image(file="image2.RData")
