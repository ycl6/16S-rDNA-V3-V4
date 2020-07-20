#

library("dada2")
library("ggplot2")
library("data.table")
library("DECIPHER")
library("phangorn")
library("phyloseq")
library("dplyr")  #"%>%"
library("ggbeeswarm")
library("ggrepel")
library("vegan")
options(width=190)

fastq = "trimmed"
filt = "filt"

fns <- sort(list.files(fastq, full.names = TRUE))
fnFs <- fns[grep(".1.fastq.gz", fns)]
fnRs <- fns[grep(".2.fastq.gz", fns)]

# Plot quality profile of fastq files
ii = 1:length(fnFs)
pdf("plotQualityProfile.pdf", width=10, height=10, pointsize=12)
for(i in ii) {
	print(plotQualityProfile(fnFs[i]) + ggtitle("Fwd"))
	print(plotQualityProfile(fnRs[i]) + ggtitle("Rev"))
}
dev.off()

if(!file_test("-d", filt)) dir.create(filt)

filtFs <- file.path(filt, basename(fnFs))
filtRs <- file.path(filt, basename(fnRs))

# Filtering and trimming
# Review "plotQualityProfile.pdf" to select the best paramters for truncLen
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
		truncLen=c(240,230), minLen=200, maxN=0, truncQ=2, maxEE=c(2,5), # Need to keep paramters consistent between runs
		rm.phix=TRUE, compress=TRUE, verbose=TRUE, multithread=TRUE)
#head(out)

# Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- gsub(".1.fastq.gz","",names(derepFs)) # Change sample name
names(derepRs) <- gsub(".2.fastq.gz","",names(derepRs)) # Change sample name

# Learn the Error Rates
dadaFs.lrn <- dada(derepFs, err=NULL, selfConsist = TRUE, multithread=TRUE)
dadaRs.lrn <- dada(derepRs, err=NULL, selfConsist = TRUE, multithread=TRUE)

pdf("plotErrors.pdf", width=10, height=10, pointsize=12)
plotErrors(dadaFs.lrn[[1]], nominalQ=TRUE)
plotErrors(dadaRs.lrn[[1]], nominalQ=TRUE)
dev.off()

# Pooled sample Inference
errF <- dadaFs.lrn[[1]]$err_out
dadaFs <- dada(derepFs, err=errF, pool=TRUE, multithread=TRUE)
errR <- dadaRs.lrn[[1]]$err_out
dadaRs <- dada(derepRs, err=errR, pool=TRUE, multithread=TRUE)

# Merge paired reads (Default: minOverlap = 20; maxMismatch = 0)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Construct a sample-by-sequence observation matrix
seqtab <- makeSequenceTable(mergers)

#> dim(seqtab)
#[1]    51 29335

# Save sequence table
saveRDS(seqtab, "seqtab.rds")

st.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#Identified 27562 bimeras out of 29335 input sequences.

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- as.data.frame(cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(st.nochim)))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- gsub(".1.fastq.gz","", rownames(track))
track$perc = (track$nonchim/track$input)*100

# We kept the majority of our raw reads, and there is no over-large drop associated with any single step
#> head(track)
#           input filtered denoisedF denoisedR merged nonchim     perc
#ERR440107 118425    93133     92982     92876  90953   86030 72.64513
#ERR440108 112545    87322     86979     87013  82564   64009 56.87414
#ERR440109 126457   101110    100824    100753  98170   90809 71.81018
#ERR440110 126010   102765    102685    102546 101526  101120 80.24760
#ERR440111 121743    97583     97479     97386  94921   92232 75.75959
#ERR440112 103303    77799     77608     77535  74405   65638 63.53930

SVformat = paste("%0",nchar(as.character(ncol(st.nochim))),"d", sep="")
svid <- paste0("SV", sprintf(SVformat, seq(ncol(st.nochim))))

# Assign taxonomy
ref1 <- "/path-to-db/silva_nr_v132_train_set.fa.gz"
ref2 <- "/path-to-db/silva_species_assignment_v132.fa.gz"
ref3 <- "/path-to-db/16SMicrobial.fa.gz"

taxtab <- assignTaxonomy(st.nochim, refFasta=ref1, minBoot=80, tryRC = TRUE, outputBootstraps=TRUE, verbose=TRUE, multithread=TRUE)
spec_silva <- assignSpecies(getSequences(st.nochim), ref2, allowMultiple = FALSE, tryRC = TRUE, verbose = TRUE)
#63 out of 1773 were assigned to the species level.
spec_ncbi <- assignSpecies(getSequences(st.nochim), ref3, allowMultiple = FALSE, tryRC = TRUE, verbose = TRUE)
#56 out of 1773 were assigned to the species level.

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
print(paste(sum(!is.na(s_final[,2])), "out of", nrow(s_final), "were assigned to the species level."))
#[1] "84 out of 1773 were assigned to the species level."

taxtab$tax[!gen.match,"Species"] <- NA
print(paste("Of which", sum(!is.na(taxtab$tax[,"Species"])),"had genera consistent with the input table."))
#[1] "Of which 70 had genera consistent with the input table."

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
#/path-to-script/raxmlHPC-PTHREADS-SSE3 -T 2 -f E -p 1234 -x 5678 -m GTRCAT -N 1 -s alignment.aln -n raxml_tree_GTRCAT
#/path-to-script/raxml-ng --evaluate --force --seed 1234 --log progress --threads 2 --msa alignment.fasta --model GTR+G --tree RAxML_fastTree.raxml_tree_GTRCAT --brlen scaled --prefix GTRCAT
#######################################################

# Import tree
raxml_tree <- read_tree("GTRCAT.raxml.bestTree")

# Load sample info
samdf <- data.frame(fread("sample.meta", colClasses = "character"))

#> head(samdf)
#  Sample.ID    Run.ID Experimental.factor Salinity Geographic.location Culture.origin
#1       R01 ERR440107                alga       32                 USA          Tampa
#2       R02 ERR440108                alga       32                 USA       Beaufort
#3       R03 ERR440109                alga       32      United Kingdom    Kingsbridge
#4       R04 ERR440110                alga       32        South Africa     Muizenberg
#5       R05 ERR440111                alga       32                 USA          Tampa
#6       R06 ERR440112                alga       32                 USA          Tampa

colnames(samdf)[1:2] = c("Sample_ID","File_ID")

rownames(samdf) = samdf$Sample_ID
samdf$Sample_ID = as.factor(samdf$Sample_ID)
rownames(st.nochim) = as.character(samdf[match(rownames(st.nochim), samdf$File_ID),]$Sample_ID) # Update sample name
samdf$Sample_ID = factor(samdf$Sample_ID, levels=c(sort(levels(samdf$Sample_ID), decreasing=F)))
samdf$Experimental.factor = as.factor(samdf$Experimental.factor)
samdf$Geographic.location = as.factor(samdf$Geographic.location)

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

#> ps
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1773 taxa and 51 samples ]
#sample_data() Sample Data:       [ 51 samples by 6 sample variables ]
#tax_table()   Taxonomy Table:    [ 1773 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1773 tips and 1771 internal nodes ]

ps0 <- subset_taxa(ps, Kingdom == "Bacteria" & !is.na(Phylum) & !is.na(Class) & Phylum!= "Cyanobacteria" & Order != "Chloroplast" & Family != "Mitochondria")

#> ps0
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1206 taxa and 51 samples ]
#sample_data() Sample Data:       [ 51 samples by 6 sample variables ]
#tax_table()   Taxonomy Table:    [ 1206 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1206 tips and 1205 internal nodes ]

source("/path-to-script/taxa_summary.R", local = TRUE)

tdt = data.table(tax_table(ps0), TotalCounts = taxa_sums(ps0), SV = taxa_names(ps0))
mdt = fast_melt(ps0)
prevdt = mdt[, list(Prevalence = sum(count > 0), TotalCounts = sum(count), MaxCounts = max(count)), by = TaxaID]

addPhylum = unique(copy(mdt[, list(TaxaID, Phylum)]))
# Join by TaxaID
setkey(prevdt, TaxaID)
setkey(addPhylum, TaxaID)
prevdt <- addPhylum[prevdt]
setkey(prevdt, Phylum)

pdf("Prevalence_TotalCounts.pdf", width=15, height=8, pointsize=12)
ggplot(prevdt, aes(Prevalence, TotalCounts, color = Phylum)) + geom_point(size = 4, alpha = 0.6) + scale_y_log10() +
facet_wrap(~Phylum, nrow = 2) + theme(legend.position="none") + xlab("Prevalence [No. Samples]") + ylab("Total Abundance")
dev.off()

#  Define prevalence threshold as 1 sample
prevalenceThreshold = 0
abundanceThreshold = 5
maxThreshold = 5

# Execute prevalence & abundance filter, using `prune_taxa()` function
keepTaxa = prevdt[(Prevalence > prevalenceThreshold & TotalCounts > abundanceThreshold & MaxCounts > maxThreshold), TaxaID]
ps1 = prune_taxa(keepTaxa, ps0)
phy_tree(ps1) = root(phy_tree(ps1), sample(taxa_names(ps1), 1), resolve.root = TRUE)

# Remove sample R63
ps1 = prune_samples(sample_names(ps1) != "R63", ps1)
ps1 = prune_taxa(taxa_sums(ps1) > 0, ps1)

#> ps1
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 606 taxa and 50 samples ]
#sample_data() Sample Data:       [ 50 samples by 6 sample variables ]
#tax_table()   Taxonomy Table:    [ 606 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 606 tips and 605 internal nodes ]

#> head(sample_data(ps1))
#Sample Data:        [6 samples by 6 sample variables]:
#    Sample_ID   File_ID Experimental.factor Salinity Geographic.location Culture.origin
#R01       R01 ERR440107                alga       32                 USA          Tampa
#R02       R02 ERR440108                alga       32                 USA       Beaufort
#R03       R03 ERR440109                alga       32      United Kingdom    Kingsbridge
#R04       R04 ERR440110                alga       32        South Africa     Muizenberg
#R05       R05 ERR440111                alga       32                 USA          Tampa
#R06       R06 ERR440112                alga       32                 USA          Tampa

# Create FASTA file
uniquesToFasta(df[rownames(df) %in% rownames(tax_table(ps1)),], "expr.otu.fasta", id=df[rownames(df) %in% rownames(tax_table(ps1)),]$id)

# Create OTU and Taxonomy tables
write.table(as.data.table(otu_table(ps1), keep.rownames=T), file="expr.otu_table.txt", 
  sep = "\t", quote = F, row.names = F, col.names = T)
write.table(as.data.table(tax_table(ps1), keep.rownames=T), file="expr.tax_table.txt", 
  sep = "\t", quote = F, row.names = F, col.names = T)

abtax = reshape2::dcast(ps1 %>% psmelt(), OTU+Kingdom+Phylum+Class+Order+Family+Genus+Species~Sample_ID, value.var = "Abundance")
write.table(abtax, file="expr.abundance.all.txt", sep = "\t", quote = F, row.names = F, col.names = T)
abtax = reshape2::dcast(ps1 %>% transform_sample_counts(function(x) {x/sum(x)} ) 
  %>% psmelt(), OTU+Kingdom+Phylum+Class+Order+Family+Genus+Species~Sample_ID, value.var = "Abundance")
write.table(abtax, file="expr.relative_abundance.all.txt", sep = "\t", quote = F, row.names = F, col.names = T)

write.table(reshape2::dcast(ps1 %>% tax_glom(taxrank = "Phylum") %>% psmelt(), Phylum~Sample_ID, value.var = "Abundance"), 
	file="expr.abundance.abphy.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(reshape2::dcast(ps1 %>% tax_glom(taxrank = "Class") %>% psmelt(), Class~Sample_ID, value.var = "Abundance"), 
	file="expr.abundance.abcls.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(reshape2::dcast(ps1 %>% tax_glom(taxrank = "Family") %>% psmelt(), Family~Sample_ID, value.var = "Abundance"), 
	file="expr.abundance.abfam.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(reshape2::dcast(ps1 %>% tax_glom(taxrank = "Genus") %>% psmelt(), Genus~Sample_ID, value.var = "Abundance"), 
	file="expr.abundance.abgen.txt", sep = "\t", quote = F, row.names = F, col.names = T)

write.table(reshape2::dcast(ps1 %>% tax_glom(taxrank = "Phylum") %>% transform_sample_counts(function(x) {x/sum(x)} ) 
  %>% psmelt(), Phylum~Sample_ID, value.var = "Abundance"), 
	file="expr.relative_abundance.abphy.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(reshape2::dcast(ps1 %>% tax_glom(taxrank = "Class") %>% transform_sample_counts(function(x) {x/sum(x)} ) 
  %>% psmelt(), Class~Sample_ID, value.var = "Abundance"), 
	file="expr.relative_abundance.abcls.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(reshape2::dcast(ps1 %>% tax_glom(taxrank = "Family") %>% transform_sample_counts(function(x) {x/sum(x)} ) 
  %>% psmelt(), Family~Sample_ID, value.var = "Abundance"), 
	file="expr.relative_abundance.abfam.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(reshape2::dcast(ps1 %>% tax_glom(taxrank = "Genus") %>% transform_sample_counts(function(x) {x/sum(x)} ) 
  %>% psmelt(), Genus~Sample_ID, value.var = "Abundance"), 
	file="expr.relative_abundance.abgen.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# Define colors for the common taxa in the taxonomic hierarchy
hexcolor.phylum = data.frame(fread("/path-to-script/hexcolor2.phylum"))
hexcolor.class = data.frame(fread("/path-to-script/hexcolor2.class"))
hexcolor.family = data.frame(fread("/path-to-script/hexcolor2.family"))

hexcolor.phylum = setNames(as.character(hexcolor.phylum$HEX), hexcolor.phylum$Phylum)
hexcolor.class = setNames(as.character(hexcolor.class$HEX), hexcolor.class$Class)
hexcolor.family = setNames(as.character(hexcolor.family$HEX), hexcolor.family$Family)

ps1.rp = transform_sample_counts(ps1, function(OTU) OTU/sum(OTU))	# Transform to proportions/relative abundances

N = 100 # Top N taxa
topN <- names(sort(taxa_sums(ps1), decreasing=TRUE))[1:N]
ps1.topN <- transform_sample_counts(ps1, function(OTU) OTU/sum(OTU))
ps1.topN <- prune_taxa(topN, ps1.topN)

pdf("plot_bar.pdf", width=10, height=12, pointsize=12)
plot_bar(ps1.rp, x="Sample_ID", fill="Phylum", title=paste(nrow(tax_table(ps1.rp)), "Taxa colored by Phylum")) + 
geom_bar(stat = "identity", size = 0.1, color = "black") + facet_wrap(~Experimental.factor, scales="free_x", nrow=1) + 
guides(fill = guide_legend(ncol = 1)) + scale_fill_manual(values=hexcolor.phylum) + scale_y_continuous(breaks = seq(0,1,0.1)) + 
theme(legend.title = element_text(size=12), legend.text = element_text(size=10), 
axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1), axis.text.y=element_text(size=10)) + 
xlab("Sample") + ylab("Relative abundance")

plot_bar(ps1.rp, x="Sample_ID", fill="Class", title=paste(nrow(tax_table(ps1.rp)), "Taxa colored by Class")) + 
geom_bar(stat = "identity", size = 0.1, color = "black") + facet_wrap(~Experimental.factor, scales="free_x", nrow=1) + 
guides(fill = guide_legend(ncol = 1)) + scale_fill_manual(values=hexcolor.class) + scale_y_continuous(breaks = seq(0,1,0.1)) + 
theme(legend.title = element_text(size=12), legend.text = element_text(size=10), 
axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1), axis.text.y=element_text(size=10)) + 
xlab("Sample") + ylab("Relative abundance")

plot_bar(ps1.topN, x="Sample_ID", fill="Phylum", title=paste("Top",N, "Taxa colored by Phylum")) + 
geom_bar(stat = "identity", size = 0.1, color = "black") + facet_wrap(~Experimental.factor, scales="free_x", nrow=1) + 
guides(fill = guide_legend(ncol = 1)) + scale_fill_manual(values=hexcolor.phylum) + scale_y_continuous(breaks = seq(0,1,0.1)) + 
theme(legend.title = element_text(size=12), legend.text = element_text(size=10), 
axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1), axis.text.y=element_text(size=10)) + 
xlab("Sample") + ylab("Relative Abundance")

plot_bar(ps1.topN, x="Sample_ID", fill="Class", title=paste("Top",N, "Taxa colored by Class")) + 
geom_bar(stat = "identity", size = 0.1, color = "black") + facet_wrap(~Experimental.factor, scales="free_x", nrow=1) + 
guides(fill = guide_legend(ncol = 1)) + scale_fill_manual(values=hexcolor.class) + scale_y_continuous(breaks = seq(0,1,0.1)) + 
theme(legend.title = element_text(size=12), legend.text = element_text(size=10), 
axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1), axis.text.y=element_text(size=10)) + 
xlab("Sample") + ylab("Relative Abundance")

plot_bar(ps1.topN, x="Sample_ID", fill="Family", title=paste("Top",N, "Taxa colored by Family")) + 
geom_bar(stat = "identity", size = 0.1, color = "black") + facet_wrap(~Experimental.factor, scales="free_x", nrow=1) + 
guides(fill = guide_legend(ncol = 1)) + scale_fill_manual(values=hexcolor.family) + scale_y_continuous(breaks = seq(0,1,0.1)) + 
theme(legend.title = element_text(size=12), legend.text = element_text(size=10), 
axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1), axis.text.y=element_text(size=10)) + 
xlab("Sample") + ylab("Relative Abundance")
dev.off()

# alpha diversity
png("plot_richness.png", width = 10, height = 6, units = "in", res = 300)
plot_richness(ps1, x="Sample_ID", measures=c("Observed","Chao1","Shannon", "Simpson"), color = "Experimental.factor", nrow=2) + 
geom_point(size=2) + theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 6, angle=90, vjust=0.5, hjust=1), 
axis.text.y = element_text(size = 10)) + labs(x = "Sample ID", y = "Alpha Diversity Measure")
dev.off()

ad = estimate_richness(ps1, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
rownames(ad) = sample_data(ps1)$Sample_ID
ad = merge(data.frame(sample_data(ps1)), ad, by = "row.names")

#> head(ad)
#  Row.names Sample_ID   File_ID Experimental.factor Salinity Geographic.location Culture.origin Observed    Chao1 se.chao1  Shannon   Simpson
#1       R01       R01 ERR440107                alga       32                 USA          Tampa      163 208.3158 18.65533 1.526441 0.4886394
#2       R02       R02 ERR440108                alga       32                 USA       Beaufort      148 262.8333 46.21770 2.324726 0.8363093
#3       R03       R03 ERR440109                alga       32      United Kingdom    Kingsbridge      121 227.2500 43.25843 1.338918 0.4935119
#4       R04       R04 ERR440110                alga       32        South Africa     Muizenberg       76 205.0000 61.23312 2.065741 0.7929804
#5       R05       R05 ERR440111                alga       32                 USA          Tampa       99 141.0000 18.79592 1.461219 0.6646137
#6       R06       R06 ERR440112                alga       32                 USA          Tampa      150 256.9375 39.56555 1.965727 0.7466008

ad = melt(ad[,c("Sample_ID","Experimental.factor","Geographic.location","Observed","Chao1","Shannon","Simpson")], id = c("Sample_ID","Experimental.factor","Geographic.location"))

location.colors = c("Australia" = "brown", "France" = "darkorange", "Peru" = "gold", "South Africa" = "green1", "Sweden" = "cyan", "United Kingdom" = "violet", "USA" = "purple")

png("plot_richness.boxplot.png", width = 7, height = 5, units = "in", res = 300)
ggplot(ad, aes(Experimental.factor, value)) + geom_boxplot(outlier.shape = NA, size = 0.8, width = 0.8) + 
geom_quasirandom(aes(color = Geographic.location), size = 1.5) + theme_bw() + facet_wrap(~ variable, scales = "free_y", nrow = 1) + scale_color_manual(values=location.colors) +
theme(legend.position="top", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
labs(x = "Experimental factor", y = "Alpha Diversity Measure")
dev.off()

# beta diversity
unifracs = GUniFrac::GUniFrac(otu_table(ps1), phy_tree(ps1), alpha=c(0, 0.5, 1))$unifracs

# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
        exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

ps1.vsd <- ps1
dds <- phyloseq_to_deseq2(ps1, ~ Experimental.factor + Geographic.location)
geoMeans = apply(counts(dds), 1, gm_mean)
dds <- DESeq2::estimateSizeFactors(dds, geoMeans = geoMeans)
dds <- estimateDispersions(dds)
abund.vsd <- getVarianceStabilizedData(dds)
abund.vsd[abund.vsd < 0] <- 0   # set negative values to 0
otu_table(ps1.vsd) <- otu_table(abund.vsd, taxa_are_rows = TRUE)

ord_UW <- ordinate(ps1, "PCoA", as.dist(unifracs[, , "d_UW"]))		# Unweighted UniFrac
ord_1 <- ordinate(ps1, "PCoA", as.dist(unifracs[, , "d_1"]))		# Weighted UniFrac
ord_VAW <- ordinate(ps1, "PCoA", as.dist(unifracs[, , "d_VAW"]))	# Variance adjusted weighted UniFrac
ord_G5 <- ordinate(ps1, "PCoA", as.dist(unifracs[, , "d_0.5"]))		# GUniFrac with alpha 0.5
ord_bray <- ordinate(ps1, "NMDS", phyloseq::distance(ps1.vsd, "bray"))	# Bray

pdf("plot_ordination.pdf", width=6, height=5, pointsize=12)
plot_ordination(ps1, ord_UW, type = "samples", color = "Experimental.factor", title = "PCoA on unweighted UniFrac") + theme_bw() + 
coord_fixed(ratio = 1) + geom_text_repel(aes(label = Sample_ID), size = 1, nudge_y = -0.01, segment.size = 0.2, box.padding = 0.1) +
stat_ellipse(aes(color=Experimental.factor, group=Experimental.factor), type = "t", size = 0.5, linetype = "dashed")

plot_ordination(ps1, ord_1, type = "samples", color = "Experimental.factor", title = "PCoA on weighted UniFrac") + theme_bw() + 
coord_fixed(ratio = 1) + geom_text_repel(aes(label = Sample_ID), size = 1, nudge_y = -0.005, segment.size = 0.2, box.padding = 0.1) +
stat_ellipse(aes(color=Experimental.factor, group=Experimental.factor), type = "t", size = 0.5, linetype = "dashed")

plot_ordination(ps1, ord_VAW, type = "samples", color = "Experimental.factor", title = "PCoA on Variance adjusted weighted UniFrac") + theme_bw() + 
coord_fixed(ratio = 1) + geom_text_repel(aes(label = Sample_ID), size = 1, nudge_y = -0.01, segment.size = 0.2, box.padding = 0.2) +
stat_ellipse(aes(color=Experimental.factor, group=Experimental.factor), type = "t", size = 0.5, linetype = "dashed")

plot_ordination(ps1, ord_G5, type = "samples", color = "Experimental.factor", title = "PCoA on GUniFrac with alpha 0.5") + theme_bw() + 
coord_fixed(ratio = 1) + geom_text_repel(aes(label = Sample_ID), size = 1, nudge_y = -0.0125, segment.size = 0.2, box.padding = 0.1) +
stat_ellipse(aes(color=Experimental.factor, group=Experimental.factor), type = "t", size = 0.5, linetype = "dashed")

plot_ordination(ps1, ord_bray, type = "samples", color = "Experimental.factor", title = "NMDS on Bray-Curtis distance (VST)") + theme_bw() + 
coord_fixed(ratio = 1) + geom_text_repel(aes(label = Sample_ID), size = 1, nudge_y = -0.015, segment.size = 0.2, box.padding = 0.1) +
stat_ellipse(aes(color=Experimental.factor, group=Experimental.factor), type = "t", size = 0.5, linetype = "dashed")
dev.off()

# Permanova test using the adonis function
# adonis tests are significant, reject the null hypothesis that our conditions have the same centroid.
adonis(as.dist(unifracs[, , "d_UW"]) ~ Experimental.factor, as(sample_data(ps1), "data.frame"), permutations = 1000)
#                    Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)
#Experimental.factor  1    0.3457 0.34568  2.9311 0.05755 0.000999 ***
#Residuals           48    5.6609 0.11793         0.94245
#Total               49    6.0065                 1.00000

adonis(as.dist(unifracs[, , "d_1"]) ~ Experimental.factor, as(sample_data(ps1), "data.frame"), permutations = 1000)
#                    Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)
#Experimental.factor  1    0.5416 0.54158  4.9349 0.09323 0.001998 **
#Residuals           48    5.2678 0.10975         0.90677
#Total               49    5.8094                 1.00000

adonis(as.dist(unifracs[, , "d_VAW"]) ~ Experimental.factor, as(sample_data(ps1), "data.frame"), permutations = 1000)
#                    Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
#Experimental.factor  1    0.1654 0.165415  2.0749 0.04144 0.0969 .
#Residuals           48    3.8267 0.079723         0.95856
#Total               49    3.9921                  1.00000

adonis(as.dist(unifracs[, , "d_0.5"]) ~ Experimental.factor, as(sample_data(ps1), "data.frame"), permutations = 1000)
#                    Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)
#Experimental.factor  1    0.5579 0.55789  3.2534 0.06348 0.000999 ***
#Residuals           48    8.2309 0.17148         0.93652
#Total               49    8.7888                 1.00000

adonis(phyloseq::distance(ps1.vsd, "bray") ~ Experimental.factor, as(sample_data(ps1), "data.frame"), permutations = 1000)
#                    Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)
#Experimental.factor  1    1.1118 1.11180  2.8532 0.05611 0.000999 ***
#Residuals           48   18.7042 0.38967         0.94389
#Total               49   19.8160                 1.00000

# Homogeneity of dispersion test
# betadisper results are not significant, meaning we cannot reject the null hypothesis that our groups have the same dispersions (i.e. having same dispersions)
# We can be more confident that our adonis result are real results, and not due to differences in group dispersions
permutest(betadisper(as.dist(unifracs[, , "d_UW"]), as(sample_data(ps1), "data.frame")$Experimental.factor), permutations = 1000)
#          Df   Sum Sq   Mean Sq      F N.Perm  Pr(>F)
#Groups     1 0.003907 0.0039068 4.3993   1000 0.03097 *
#Residuals 48 0.042626 0.0008880

permutest(betadisper(as.dist(unifracs[, , "d_1"]), as(sample_data(ps1), "data.frame")$Experimental.factor), permutations = 1000)
#          Df  Sum Sq  Mean Sq      F N.Perm  Pr(>F)
#Groups     1 0.05794 0.057936 4.3123   1000 0.04595 *
#Residuals 48 0.64489 0.013435

permutest(betadisper(as.dist(unifracs[, , "d_VAW"]), as(sample_data(ps1), "data.frame")$Experimental.factor), permutations = 1000)
#          Df  Sum Sq  Mean Sq      F N.Perm   Pr(>F)
#Groups     1 0.30256 0.302557 19.066   1000 0.000999 ***
#Residuals 48 0.76171 0.015869

permutest(betadisper(as.dist(unifracs[, , "d_0.5"]), as(sample_data(ps1), "data.frame")$Experimental.factor), permutations = 1000)
#          Df  Sum Sq  Mean Sq      F N.Perm  Pr(>F)
#Groups     1 0.03814 0.038140 5.3178   1000 0.03596 *
#Residuals 48 0.34427 0.007172

permutest(betadisper(phyloseq::distance(ps1.vsd, "bray"), as(sample_data(ps1), "data.frame")$Experimental.factor), permutations = 1000)
#          Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.01236 0.0123639 1.7611   1000 0.2188
#Residuals 48 0.33698 0.0070205

source("/path-to-script/lefse.R", local = TRUE)

# Create directory if not exist
if(!dir.exists("lefse")) { dir.create("lefse") }

tax1 = lefse_1name_obj(ps1, sample_data(ps1)$Experimental.factor)
lefse1 = lefse_obj(ps1)
lefse1 = rbind(tax1, lefse1)
lefse1$name = gsub("-","_",lefse1$name)
lefse1$name = gsub("Escherichia/Shigella","Escherichia_Shigella",lefse1$name)
lefse1$name = gsub("Bacteria\\|Actinobacteria","Bacteria|Actinobacteria_P",lefse1$name)
lefse1$name = gsub("Bacteria\\|Actinobacteria_P\\|Actinobacteria","Bacteria|Actinobacteria_P|Actinobacteria_C",lefse1$name)
write.table(lefse1, file="lefse/expr1.lefse_table.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#######################################################
# Run below commands outside R in Linux environment
# /path-to-script/nsegata-lefse-82605a2ae7b7/format_input.py lefse/expr1.lefse_table.txt lefse/expr1.lefse_table.in -c 1 -u 2 -o 1000000
# /path-to-script/nsegata-lefse-82605a2ae7b7/run_lefse.py lefse/expr1.lefse_table.in lefse/expr1.lefse_table.res -b 100 -a 0.05 -w 0.05 -l 2
#Number of significantly discriminative features: 51 ( 51 ) before internal wilcoxon
#Number of discriminative features with abs LDA score > 2.0 : 51

# create a "lefse.expr.colors" file to define group color (see example file)

# /path-to-script/graphlan_commit_6ca8735/export2graphlan/export2graphlan.py -i lefse/expr1.lefse_table.txt -o lefse/expr1.lefse_table.res -t lefse/expr1.graphlan_tree.txt -a lefse/expr1.graphlan_annot.txt --external_annotations 2,3,4,5,6 --fname_row 0 --skip_rows 1 --biomarkers2colors lefse/lefse.expr.colors
# /path-to-script/graphlan_commit_6ca8735/graphlan_annotate.py --annot lefse/expr1.graphlan_annot.txt lefse/expr1.graphlan_tree.txt lefse/expr1.graphlan_outtree.txt
# /path-to-script/graphlan_commit_6ca8735/graphlan.py --dpi 150 lefse/expr1.graphlan_outtree.txt lefse/expr1.graphlan.png --external_legends --size 8 --pad 0.2
# perl /path-to-script/lefse.pl lefse/expr1.lefse_table.res lefse/expr1.graphlan_outtree.txt lefse/expr1.out
#######################################################

# plot res in R
library(grid)
Sys.setlocale("LC_COLLATE","C")

res1 = data.frame(fread("lefse/expr1.out"))
res1 = res1[order(res1$order),]
names(res1)[5] = "Group"
res1$taxon = paste0(res1$taxon, "(",res1$rank, ";", " ",res1$order, ")")
res1$taxon = as.factor(res1$taxon)
res1$taxon = factor(res1$taxon, levels = unique(res1$taxon[order(res1$order, decreasing=T)]))
res1$Group = as.factor(res1$Group)
levels(res1$taxon) = gsub("Escherichia_Shigella","Escherichia/Shigella",levels(res1$taxon))
levels(res1$taxon) = gsub("_UCG_","_UCG-",levels(res1$taxon))

#> head(res1)
#                                                                            fulltaxon                    taxon order rank  Group      lda     pvalue logmaxpct
#31            Bacteria.Actinobacteria_P.Actinobacteria_C.Micrococcales.Micrococcaceae     Micrococcaceae(F; A)     A    F   alga 3.231928 0.02321769  1.727602
#36     Bacteria.Actinobacteria_P.Actinobacteria_C.Micrococcales.Micrococcaceae.Rothia             Rothia(G; B)     B    G   alga 3.632929 0.02229896  1.440856
#43                       Bacteria.Actinobacteria_P.Actinobacteria_C.Pseudonocardiales  Pseudonocardiales(O; C)     C    O   alga 3.200686 0.04604858  1.512401
#15    Bacteria.Actinobacteria_P.Actinobacteria_C.Pseudonocardiales.Pseudonocardiaceae Pseudonocardiaceae(F; D)     D    F   alga 3.286615 0.04604858  1.512401
#37       Bacteria.Bacteroidetes.Bacteroidia.Cytophagales.Cyclobacteriaceae.Fabibacter         Fabibacter(G; E)     E    G medium 3.648937 0.01492188  3.937590
#13 Bacteria.Bacteroidetes.Bacteroidia.Flavobacteriales.Crocinitomicaceae.Salinirepens       Salinirepens(G; F)     F    G medium 3.617960 0.02000840  1.835670

plot_lefse = ggplot(res1, aes(taxon, lda, fill = Group)) + geom_bar(stat = "identity", width = 0.7, size = 0.5) + coord_flip() + theme_bw() +
facet_wrap(~ Group, ncol = 1, scales = "free_y") + theme(legend.position="right", axis.text.x=element_text(size = 12), axis.text.y=element_text(face = "bold", size = 8),
strip.text.x = element_text(face = "bold", size = 12)) + labs(title = "Linear discriminant analysis (LDA)\nEffect Size Analysis", x = "Taxon", y = "LDA")
groupN = res1 %>% group_by(Group) %>% summarise(count = length(unique(taxon)))
gt = ggplotGrob(plot_lefse)
panelI.1 <- gt$layout$t[grepl("panel", gt$layout$name)]
gt$heights[panelI.1] <- unit(groupN$count,"null")
dev.off()

png("lefse/expr1.lefse_table.png", width = 8, height = 8, units = "in", res = 300)
grid.draw(gt)
dev.off()

##### Signifiant genera based on LEfSe results #####
abgen = ps1 %>% tax_glom(taxrank = "Genus") %>% transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt()
abgen = abgen[abgen$Genus %in% gsub("\\(.+","",as.character(res1[res1$rank == "G",]$taxon)),]

#> length(unique(as.character(abgen$Genus)))
#[1] 24

library(plyr)

anno1 = ddply(abgen, .(Genus), summarise, y = max(Abundance)*100*1.1)
anno1 = merge(anno1, cbind(Genus = gsub("\\(.+","",as.character(res1[res1$rank == "G",]$taxon)), res1[res1$rank == "G",6:7]), by = "Genus")
anno1$label = paste0("P=",formatC(anno1$pvalue, format = "f"), " (LDA=",formatC(anno1$lda, format = "f"),")")

png("expr1.lefse.boxplot.png", width=12, height=10, units = "in", res = 300)
ggplot(abgen, aes(Experimental.factor, Abundance*100, color = Experimental.factor)) + geom_boxplot(outlier.shape = NA, size = 0.8, width = 0.8) +
geom_quasirandom(size = 1, color = "black") + theme_bw() + facet_wrap(~ Genus, scales = "free_y", nrow = 4) +
theme(legend.position="none", axis.title.x = element_blank(), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
strip.text.x = element_text(face = "bold", size = 8)) + labs(y = "Relative abundance (%)") +
geom_text(data=anno1, aes(x = 1.5, y = y, label = label), size = 3, color = "black")
dev.off()

##### Signifiant species #####
abspc = ps1 %>% tax_glom(taxrank = "Species") %>% transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt()
abspc$Species = paste(abspc$Genus, abspc$Species)
abspc$Species = gsub("_[a-zA-Z0-9_]*","",abspc$Species)

#> length(unique(as.character(abspc$Species)))
#[1] 59

result <- data.frame()
for (Species in unique(as.character(abspc$Species))) {
        test.result = kruskal.test(Abundance ~ Experimental.factor, abspc[abspc$Species == Species,])
        result[Species, 'p'] <- test.result$p.value
}
result$p.adjusted <- p.adjust(result$p, method="fdr")

#> table(result$p < 0.05)
#FALSE  TRUE
#   46    13

#> table(result$p.adjusted < 0.1)
#FALSE  TRUE
#   57     2

abspc = abspc[abspc$Species %in% rownames(result[result$p < 0.05,]),]

anno2 = ddply(abspc, .(Species), summarise, y = max(Abundance)*100*1.1)
anno2 = merge(anno2, result[result$p < 0.05,], by.x = "Species", by.y = "row.names")
anno2$label = paste0("P=",formatC(anno2$p, format = "f"))

png("species_targets.boxplot.png", width=12, height=8, units = "in", res = 300)
ggplot(abspc, aes(Experimental.factor, Abundance*100, color = Experimental.factor)) + geom_boxplot(outlier.shape = NA, size = 0.8, width = 0.8) +
geom_quasirandom(size = 1, color = "black") + theme_bw() + facet_wrap(~ Species, scales = "free_y", nrow = 3) +
theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_text(size = 10), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
strip.text.x = element_text(face = "bold", size = 8)) + labs(title = "Candidate bacterial species", y = "Relative abundance (% species-level assigned)") +
geom_text(data=anno2, aes(x = 1.5, y = y, label = label), size = 3, color = "black")
dev.off()

# Save current workspace
# save.image(file="image.RData")
