#

library("data.table")
library("phyloseq")
library("ggplot2")
library("dplyr")  #"%>%"
library("vegan")
options(width=190)

# Continued from DADA2 part 2 (dada2-create-phyloseq-obj.R)

ps0 <- subset_taxa(ps, Kingdom == "Bacteria" & !is.na(Phylum) & !is.na(Class) & Phylum!= "Cyanobacteria" & Order != "Chloroplast" & Family != "Mitochondria")

source("taxa_summary.R", local = TRUE)
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

#  Define prevalence threshold as 1 samples
prevalenceThreshold = 0
abundanceThreshold = 5
maxThreshold = 5

# Execute prevalence & abundance filter, using `prune_taxa()` function
keepTaxa = prevdt[(Prevalence > prevalenceThreshold & TotalCounts > abundanceThreshold & MaxCounts > maxThreshold), TaxaID]
ps1 = prune_taxa(keepTaxa, ps0)
phy_tree(ps1) = root(phy_tree(ps1), sample(taxa_names(ps1), 1), resolve.root = TRUE)

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
hexcolor.phylum = data.frame(fread("hexcolor.phylum"))
hexcolor.class = data.frame(fread("hexcolor.class"))
hexcolor.family = data.frame(fread("hexcolor.family"))
hexcolor.genus = data.frame(fread("hexcolor.genus"))

hexcolor.phylum = setNames(as.character(hexcolor.phylum$HEX), hexcolor.phylum$Phylum)
hexcolor.class = setNames(as.character(hexcolor.class$HEX), hexcolor.class$Class)
hexcolor.family = setNames(as.character(hexcolor.family$HEX), hexcolor.family$Family)
hexcolor.genus = setNames(as.character(hexcolor.genus$HEX), hexcolor.genus$Genus)

ps1.rp = transform_sample_counts(ps1, function(OTU) OTU/sum(OTU))	# Transform to proportions/relative abundances

N = 200
topN <- names(sort(taxa_sums(ps1), decreasing=TRUE))[1:N]
ps1.topN <- transform_sample_counts(ps1, function(OTU) OTU/sum(OTU))
ps1.topN <- prune_taxa(topN, ps1.topN)

pdf("plot_bar.pdf", width=10, height=12, pointsize=12)
plot_bar(ps1.rp, x="Sample_ID", fill="Phylum", title=paste(nrow(tax_table(ps1.rp)), "Taxa colored by Phylum")) + 
geom_bar(stat = "identity", size = 0.1, color = "black") + facet_wrap(~Group, scales="free_x", nrow=1) + 
guides(fill = guide_legend(ncol = 1)) + scale_fill_manual(values=hexcolor.phylum) + scale_y_continuous(breaks = seq(0,1,0.1)) + 
theme(legend.title = element_text(size=12), legend.text = element_text(size=10), 
axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1), axis.text.y=element_text(size=10)) + 
xlab("Sample") + ylab("Relative abundance")

plot_bar(ps1.rp, x="Sample_ID", fill="Class", title=paste(nrow(tax_table(ps1.rp)), "Taxa colored by Class")) + 
geom_bar(stat = "identity", size = 0.1, color = "black") + facet_wrap(~Group, scales="free_x", nrow=1) + 
guides(fill = guide_legend(ncol = 1)) + scale_fill_manual(values=hexcolor.class) + scale_y_continuous(breaks = seq(0,1,0.1)) + 
theme(legend.title = element_text(size=12), legend.text = element_text(size=10), 
axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1), axis.text.y=element_text(size=10)) + 
xlab("Sample") + ylab("Relative abundance")

plot_bar(ps1.topN, x="Sample_ID", fill="Phylum", title=paste("Top",N, "Taxa colored by Phylum")) + 
geom_bar(stat = "identity", size = 0.1, color = "black") + facet_wrap(~Group, scales="free_x", nrow=1) + 
guides(fill = guide_legend(ncol = 1)) + scale_fill_manual(values=hexcolor.phylum) + scale_y_continuous(breaks = seq(0,1,0.1)) + 
theme(legend.title = element_text(size=12), legend.text = element_text(size=10), 
axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1), axis.text.y=element_text(size=10)) + 
xlab("Sample") + ylab("Relative Abundance")

plot_bar(ps1.topN, x="Sample_ID", fill="Class", title=paste("Top",N, "Taxa colored by Class")) + 
geom_bar(stat = "identity", size = 0.1, color = "black") + facet_wrap(~Group, scales="free_x", nrow=1) + 
guides(fill = guide_legend(ncol = 1)) + scale_fill_manual(values=hexcolor.class) + scale_y_continuous(breaks = seq(0,1,0.1)) + 
theme(legend.title = element_text(size=12), legend.text = element_text(size=10), 
axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1), axis.text.y=element_text(size=10)) + 
xlab("Sample") + ylab("Relative Abundance")

plot_bar(ps1.topN, x="Sample_ID", fill="Family", title=paste("Top",N, "Taxa colored by Family")) + 
geom_bar(stat = "identity", size = 0.1, color = "black") + facet_wrap(~Group, scales="free_x", nrow=1) + 
guides(fill = guide_legend(ncol = 1)) + scale_fill_manual(values=hexcolor.family) + scale_y_continuous(breaks = seq(0,1,0.1)) + 
theme(legend.title = element_text(size=12), legend.text = element_text(size=10), 
axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1), axis.text.y=element_text(size=10)) + 
xlab("Sample") + ylab("Relative Abundance")

plot_bar(ps1.topN, x="Sample_ID", fill="Genus", title=paste("Top",N, "Taxa colored by Genus")) + 
geom_bar(stat = "identity", size = 0.1, color = "black") + facet_wrap(~Group, scales="free_x", nrow=1) + 
guides(fill = guide_legend(ncol = 1)) + scale_fill_manual(values=hexcolor.genus) + scale_y_continuous(breaks = seq(0,1,0.1)) + 
theme(legend.title = element_text(size=12), legend.text = element_text(size=10), 
axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1), axis.text.y=element_text(size=10)) + 
xlab("Sample") + ylab("Relative Abundance")
dev.off()

# alpha diversity
png("plot_richness.png", width = 10, height = 5, units = "in", res = 300)
plot_richness(ps1, x="Sample_ID", measures=c("Observed","Chao1","Shannon", "Simpson"), color = "Group", nrow=1) + 
geom_point(size=2) + theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 10), 
axis.text.y = element_text(size = 10)) + labs(x = "Sample ID", y = "Alpha Diversity Measure")
dev.off()

ad = estimate_richness(ps1, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
rownames(ad) = sample_data(ps1)$Sample_ID
ad = merge(data.frame(sample_data(ps1)), ad, by = "row.names")
ad = melt(ad[,c("Sample_ID","Group","Observed","Chao1","Shannon","Simpson")], id = c("Sample_ID","Group"))

png("plot_richness.boxplot.png", width = 10, height = 5, units = "in", res = 300)
ggplot(ad, aes(Group, value, color = Group)) + geom_boxplot(outlier.shape = NA, size = 0.8, width = 0.8) + 
geom_quasirandom(size = 1, color = "black") + theme_bw() + facet_wrap(~ variable, scales = "free_y", nrow = 1) + 
theme(legend.position="top", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
labs(x = "Group", y = "Alpha Diversity Measure")
dev.off()

# beta diversity
unifracs = GUniFrac::GUniFrac(otu_table(ps1), phy_tree(ps1), alpha=c(0, 0.5, 1))$unifracs

# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
        exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

ps1.vsd <- ps1
dds <- phyloseq_to_deseq2(ps1, ~ Group + Batch)
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
plot_ordination(ps1, ord_UW, type = "samples", color = "Group", title = "PCoA on unweighted UniFrac") + theme_bw() + 
coord_fixed(ratio = 1) + geom_text_repel(aes(label = Sample_ID), size = 2, nudge_y = -0.01, segment.size = 0.2, box.padding = 0.1) +
stat_ellipse(aes(color=Group, group=Group), type = "t", size = 0.5, linetype = "dashed")

plot_ordination(ps1, ord_1, type = "samples", color = "Group", title = "PCoA on weighted UniFrac") + theme_bw() + 
coord_fixed(ratio = 1) + geom_text_repel(aes(label = Sample_ID), size = 2, nudge_y = -0.005, segment.size = 0.2, box.padding = 0.1) +
stat_ellipse(aes(color=Group, group=Group), type = "t", size = 0.5, linetype = "dashed")

plot_ordination(ps1, ord_VAW, type = "samples", color = "Group", title = "PCoA on Variance adjusted weighted UniFrac") + theme_bw() + 
coord_fixed(ratio = 1) + geom_text_repel(aes(label = Sample_ID), size = 2, nudge_y = -0.01, segment.size = 0.2, box.padding = 0.2) +
stat_ellipse(aes(color=Group, group=Group), type = "t", size = 0.5, linetype = "dashed")

plot_ordination(ps1, ord_G5, type = "samples", color = "Group", title = "PCoA on GUniFrac with alpha 0.5") + theme_bw() + 
coord_fixed(ratio = 1) + geom_text_repel(aes(label = Sample_ID), size = 2, nudge_y = -0.0125, segment.size = 0.2, box.padding = 0.1) +
stat_ellipse(aes(color=Group, group=Group), type = "t", size = 0.5, linetype = "dashed")

plot_ordination(ps1, ord_bray, type = "samples", color = "Group", title = "NMDS on Bray-Curtis distance (VST)") + theme_bw() + 
coord_fixed(ratio = 1) + geom_text_repel(aes(label = Sample_ID), size = 2, nudge_y = -0.015, segment.size = 0.2, box.padding = 0.1) +
stat_ellipse(aes(color=Group, group=Group), type = "t", size = 0.5, linetype = "dashed")
dev.off()

# Permanova test using the adonis function
# adonis tests are significant, reject the null hypothesis that our conditions have the same centroid.
adonis(as.dist(unifracs[, , "d_UW"]) ~ Group, as(sample_data(ps1), "data.frame"), permutations = 1000)
adonis(as.dist(unifracs[, , "d_1"]) ~ Group, as(sample_data(ps1), "data.frame"), permutations = 1000)
adonis(as.dist(unifracs[, , "d_VAW"]) ~ Group, as(sample_data(ps1), "data.frame"), permutations = 1000)
adonis(as.dist(unifracs[, , "d_0.5"]) ~ Group, as(sample_data(ps1), "data.frame"), permutations = 1000)
adonis(phyloseq::distance(ps1.vsd, "bray") ~ Group, as(sample_data(ps1), "data.frame"), permutations = 1000)

# Homogeneity of dispersion test
# betadisper results are not significant, meaning we cannot reject the null hypothesis that our groups have the same dispersions (i.e. having same dispersions)
# We can be more confident that our adonis result are real results, and not due to differences in group dispersions
permutest(betadisper(as.dist(unifracs[, , "d_UW"]), as(sample_data(ps1), "data.frame")$Group), permutations = 1000)
permutest(betadisper(as.dist(unifracs[, , "d_1"]), as(sample_data(ps1), "data.frame")$Group), permutations = 1000)
permutest(betadisper(as.dist(unifracs[, , "d_VAW"]), as(sample_data(ps1), "data.frame")$Group), permutations = 1000)
permutest(betadisper(as.dist(unifracs[, , "d_0.5"]), as(sample_data(ps1), "data.frame")$Group), permutations = 1000)
permutest(betadisper(phyloseq::distance(ps1.vsd, "bray"), as(sample_data(ps1), "data.frame")$Group), permutations = 1000)

# Save current workspace
# save.image(file="image3.RData")
