#

library("dada2")
library("data.table")
library("phyloseq")
library("ggplot2")
options(width=190)

# Continued from DADA2 part 1 (dada2-per-run-processing.R or cutadapt-and-dada2-per-run-processing.R)

# Load saved workspace
# save.image(file="image.RData")

# Load sequence table from multiple runs
s01 = readRDS("seqtab1.rds")
s02 = readRDS("seqtab2.rds")

# Merge sequence table and remove chimeras
st.all = mergeSequenceTables(s01, s02)

# If this script is used independently and the "sample.names" is not carried over from previous step/script
sample.names = gsub(".1.fastq.gz", "", rownames(st.all))

# To sum values of same sample from multiple sequence table (i.e. when a sample was re-sequenced due to low depth)
# st.all = mergeSequenceTables(s01, s02, repeats = "sum")

# Remove bimeras (two-parent chimeras)
st.nochim = removeBimeraDenovo(st.all, verbose = TRUE, multithread = TRUE)
rownames(st.nochim) = sample.names

dim(seqtab)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# Track reads through the pipeline
# If you load saved workspace from previous step/script, you will have the necessary objects to build the track data.frame
# getN <- function(x) sum(getUniques(x))
# track = cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# colnames(track) = c("Trimmed", "Filtered", "denoisedF", "denoisedR", "merged", "nonchim")
# track = cbind(data.frame(SampleID = sample.names), track)
# write.table(track, file = "track.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# Assign taxonomy
# Note to change the PATH to reference training datasets accordingly
dbpath = "/path-to-db/"
ref1 = paste0(dbpath, "silva_nr_v138_train_set.fa.gz")
ref2 = paste0(dbpath, "silva_species_assignment_v138.fa.gz")
ref3 = paste0(dbpath, "16SMicrobial.fa.gz")

# Classifies sequences against SILVA reference training dataset
taxtab = assignTaxonomy(st.nochim, refFasta = ref1, minBoot = 80, tryRC = TRUE, outputBootstraps = TRUE, verbose = TRUE, multithread = TRUE)

# Taxonomic assignment to the species level by exact matching against SILVA and NCBI reference datasets
spec_silva = assignSpecies(getSequences(st.nochim), ref2, allowMultiple = FALSE, tryRC = TRUE, verbose = TRUE)
spec_ncbi = assignSpecies(getSequences(st.nochim), ref3, allowMultiple = FALSE, tryRC = TRUE, verbose = TRUE)

# Combine species-level taxonomic assignment from 2 reference sources
SVformat = paste("%0",nchar(as.character(ncol(st.nochim))),"d", sep = "")
svid = paste0("SV", sprintf(SVformat, seq(ncol(st.nochim))))

s_silva = as.data.frame(spec_silva, stringsAsFactors = FALSE)
rownames(s_silva) = svid

s_ncbi = as.data.frame(spec_ncbi, stringsAsFactors = FALSE)
rownames(s_ncbi) = svid
s_ncbi$Genus = gsub("\\[|\\]", "", s_ncbi$Genus)

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
        gcol = which(colnames(taxtab$tax) == "Genus")
} else { 
	gcol = ncol(taxtab$tax) 
}

matchGenera <- function(gen.tax, gen.binom, split.glyph = "/") {
	if(is.na(gen.tax) || is.na(gen.binom)) { return(FALSE) }
	if((gen.tax == gen.binom) || grepl(paste0("^", gen.binom, "[ _", split.glyph, "]"), gen.tax) || grepl(paste0(split.glyph, gen.binom, "$"), gen.tax)) {
		return(TRUE)
	} else {
		return(FALSE)
	}
}

gen.match = mapply(matchGenera, taxtab$tax[,gcol], s_final[,1])
taxtab$tax = cbind(taxtab$tax, s_final[,2])
colnames(taxtab$tax)[ncol(taxtab$tax)] = "Species"
print(paste(sum(!is.na(s_final[,2])), "out of", nrow(s_final), "were assigned to the species level."))

taxtab$tax[!gen.match,"Species"] = NA
print(paste("Of which", sum(!is.na(taxtab$tax[,"Species"])), "had genera consistent with the input table."))

# Prepare df
df = data.frame(sequence = colnames(st.nochim), abundance = colSums(st.nochim), stringsAsFactors = FALSE)
df$id = svid

df = merge(df, as.data.frame(taxtab), by = "row.names")
rownames(df) = df$id
df = df[order(df$id),2:ncol(df)]

# Construct phylogenetic tree
alignment = DECIPHER::AlignSeqs(Biostrings::DNAStringSet(setNames(df$sequence, df$id)), anchor = NA)

# Export alignment
phang.align = phangorn::phyDat(as(alignment, "matrix"), type = "DNA")
phangorn::write.phyDat(phang.align, file = "alignment.fasta", format = "fasta")
phangorn::write.phyDat(phang.align, file = "alignment.aln", format = "phylip")

# Phylogenetic tree construction
# Use system2() to execute commands in R
# Note to change the PATH to raxml and raxmlng accordingly
raxml = "/path-to-script/raxmlHPC-PTHREADS-SSE3"
raxmlng = "/path-to-script/raxml-ng"

system2(raxml, args = c("-T 2", "-f E", "-p 1234", "-x 5678", "-m GTRCAT", "-N 1", "-s alignment.aln", "-n raxml_tree_GTRCAT"))
system2(raxmlng, args = c("--evaluate", "--force", "--seed 1234", "--log progress", "--threads 2", "--msa alignment.fasta", "--model GTR+G", "--tree RAxML_fastTree.raxml_tree_GTRCAT", "--brlen scaled", "--prefix GTRCAT"))

# Import tree
# raxml-ng (v0.7.0): GTRCAT.raxml.mlTrees
# raxml-ng (v0.9.0): GTRCAT.raxml.bestTree
raxml_tree = read_tree("GTRCAT.raxml.bestTree")

# Load sample info (sample.meta)
samdf = data.frame(fread("sample.meta", colClasses = "character"))

# Example content:
#  Sample_ID     Run_ID    Batch    Group
#    Sample1   ERR00001     run1        A
#    Sample2   ERR00002     run1        B
#    Sample3   ERR00003     run2        A
#    Sample4   ERR00004     run2        B

rownames(samdf) = samdf$Sample_ID
samdf$Sample_ID = as.factor(samdf$Sample_ID)
samdf$Sample_ID = factor(samdf$Sample_ID, levels = c(sort(levels(samdf$Sample_ID), decreasing = F)))
samdf$Batch = as.factor(samdf$Batch)	# Encode a vector as a factor, update the factor level if necessary
samdf$Group = as.factor(samdf$Group)	# Encode a vector as a factor, update the factor level if necessary

#head(samdf, 10)

# Handoff to phyloseq
new_seqtab = st.nochim
colnames(new_seqtab) = df[match(colnames(new_seqtab), df$sequence),]$id

# Update rownames in new_seqtab from Run_ID to Sample_ID
rownames(new_seqtab) = as.character(samdf[match(rownames(new_seqtab), samdf$Run_ID),]$Sample_ID)

new_taxtab = taxtab
rownames(new_taxtab$tax) = df[match(rownames(new_taxtab$tax), df$sequence),]$id

tax = as.data.frame(new_taxtab$tax)
tax$Family = as.character(tax$Family)
tax$Genus = as.character(tax$Genus)

# Combine data into a phyloseq object
ps = phyloseq(tax_table(as.matrix(tax)), sample_data(samdf), otu_table(new_seqtab, taxa_are_rows = FALSE), phy_tree(raxml_tree))

ps

# Save current workspace
# save.image(file="image2.RData")
