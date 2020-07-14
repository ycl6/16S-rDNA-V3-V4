#

library("data.table")
library("phyloseq")
library("DESeq2")
library("ggplot2")
library("ggbeeswarm")
library("ggrepel")
library("vegan")
library("dplyr")
options(width=190)

# Continued from DADA2 part 2 (dada2-create-phyloseq-obj.R)

# Load saved workspace
# save.image(file="image2.RData")

# Note to change the PATH to the "taxa_summary.R" script accordingly
source("/path-to-script/taxa_summary.R", local = TRUE)

ps0 = subset_taxa(ps, Kingdom == "Bacteria" & !is.na(Phylum) & !is.na(Class) & Family != "Mitochondria")

# Create a total counts data.table
tdt = data.table(tax_table(ps0), TotalCounts = taxa_sums(ps0), SV = taxa_names(ps0))
tdt

ggplot(tdt, aes(TotalCounts)) + geom_histogram(bins = 50) + theme_bw() + 
	ggtitle("Histogram of Total Counts")

# taxa cumulative sum
taxcumsum = tdt[, .N, by = TotalCounts]
setkey(taxcumsum, TotalCounts)
taxcumsum[, CumSum := cumsum(N)]

# Plot the cumulative sum of ASVs against the total counts
pCumSum = ggplot(taxcumsum, aes(TotalCounts, CumSum)) + geom_point() + theme_bw() + 
	xlab("Filtering Threshold") + ylab("ASV Filtered")

gridExtra::grid.arrange(pCumSum, pCumSum + xlim(0, 500), 
	pCumSum + xlim(0, 100), pCumSum + xlim(0, 50), nrow = 2, 
	top = "ASVs that would be filtered vs. minimum taxa counts threshold")

# Create a prevalence data.table
mdt = fast_melt(ps0)
mdt

prevdt = mdt[, list(Prevalence = sum(count > 0), TotalCounts = sum(count), MaxCounts = max(count)), by = TaxaID]
prevdt

ggplot(prevdt, aes(Prevalence)) + geom_histogram(bins = 50) + theme_bw() + 
	ggtitle("Histogram of Taxa Prevalence")

ggplot(prevdt, aes(MaxCounts)) + geom_histogram(bins = 100) + xlim(0, 500) + theme_bw() + 
	ggtitle("Histogram of Maximum TotalCounts")

# taxa cumulative sum
prevcumsum = prevdt[, .N, by = Prevalence]
setkey(prevcumsum, Prevalence)
prevcumsum[, CumSum := cumsum(N)]

# Plot the cumulative sum of ASVs against the prevalence
pPrevCumSum = ggplot(prevcumsum, aes(Prevalence, CumSum)) + geom_point(size = 2, alpha = 0.5) + 
	theme_bw() + xlab("Filtering Threshold") + ylab("ASVs Filtered") + 
	ggtitle("ASVs that would be filtered vs. minimum sample count threshold")

pPrevCumSum

# Plot Prevalence vs. Total Count
pt1 = ggplot(prevdt, aes(Prevalence, TotalCounts)) + geom_point(size = 2, alpha = 0.5) + 
	scale_y_log10() + theme_bw() + xlab("Prevalence [No. Samples]") + ylab("TotalCounts [Taxa]")

# Plot Prevalence vs. Total Count (colored by phylum)
addPhylum = unique(copy(mdt[, list(TaxaID, Phylum)]))

# Join by TaxaID
setkey(prevdt, TaxaID)
setkey(addPhylum, TaxaID)
prevdt = addPhylum[prevdt]
setkey(prevdt, Phylum)

pt2 = ggplot(prevdt, aes(Prevalence, TotalCounts, color = Phylum)) + geom_point(size = 1, alpha = 0.5) + 
	scale_y_log10() + theme_bw() + facet_wrap(~Phylum, nrow = 3) + theme(legend.position="none") + 
	xlab("Prevalence [No. Samples]") + ylab("Total Abundance")

pdf("Prevalence_TotalCounts.pdf", width = 15, height = 8, pointsize = 12)
print(pt1)
print(pt2)
dev.off()

# Define filter threshold
# Review the Prevalence, Total Count and MaxCounts plots to select the best paramters
prevalenceThreshold = 5
abundanceThreshold = 10
maxThreshold = 5

# Execute prevalence & abundance filter, using `prune_taxa()` function

keepTaxa = prevdt[(Prevalence > prevalenceThreshold & TotalCounts > abundanceThreshold & MaxCounts > maxThreshold), TaxaID]
ps1 = prune_taxa(keepTaxa, ps0)
phy_tree(ps1) = ape::root(phy_tree(ps1), sample(taxa_names(ps1), 1), resolve.root = TRUE)

ps1

# Create FASTA file
# The df data.frame is created in previous step
dada2::uniquesToFasta(df[rownames(df) %in% taxa_names(ps1),], "expr.asv.fasta", ids = df[rownames(df) %in% taxa_names(ps1),]$id)

# Create BIOM file
biomformat::write_biom(biomformat::make_biom(data = t(as.matrix(otu_table(ps1)))), "expr.biom")

# Create ASV and Taxonomy tables
write.table(as.data.table(otu_table(ps1), keep.rownames=T), file="expr.otu_table.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(as.data.table(tax_table(ps1), keep.rownames=T), file="expr.tax_table.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# Create raw abundance tables
write.table(reshape2::dcast(ps1 %>% psmelt() %>% arrange(OTU) %>% rename(ASV = OTU), 
	ASV+Kingdom+Phylum+Class+Order+Family+Genus+Species~Sample_ID, value.var = "Abundance"), 
	file = "expr.abundance.all.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(reshape2::dcast(ps1 %>% tax_glom(taxrank = "Phylum") %>% psmelt(), Phylum~Sample_ID, value.var = "Abundance"), 
	file = "expr.abundance.abphy.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(reshape2::dcast(ps1 %>% tax_glom(taxrank = "Class") %>% psmelt(), Class~Sample_ID, value.var = "Abundance"), 
	file = "expr.abundance.abcls.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(reshape2::dcast(ps1 %>% tax_glom(taxrank = "Family") %>% psmelt(), Family~Sample_ID, value.var = "Abundance"), 
	file = "expr.abundance.abfam.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(reshape2::dcast(ps1 %>% tax_glom(taxrank = "Genus") %>% psmelt(), Genus~Sample_ID, value.var = "Abundance"), 
	file = "expr.abundance.abgen.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# Create relative abundance tables
write.table(reshape2::dcast(ps1 %>% transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt() %>% arrange(OTU) %>% rename(ASV = OTU), 
	ASV+Kingdom+Phylum+Class+Order+Family+Genus+Species~Sample_ID, value.var = "Abundance"), 
	file = "expr.relative_abundance.all.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(reshape2::dcast(ps1 %>% tax_glom(taxrank = "Phylum") %>% transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt(), 
	Phylum~Sample_ID, value.var = "Abundance"), file = "expr.relative_abundance.abphy.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(reshape2::dcast(ps1 %>% tax_glom(taxrank = "Class") %>% transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt(), 
	Class~Sample_ID, value.var = "Abundance"), file = "expr.relative_abundance.abcls.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(reshape2::dcast(ps1 %>% tax_glom(taxrank = "Family") %>% transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt(), 
	Family~Sample_ID, value.var = "Abundance"), file = "expr.relative_abundance.abfam.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(reshape2::dcast(ps1 %>% tax_glom(taxrank = "Genus") %>% transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt(), 
	Genus~Sample_ID, value.var = "Abundance"), file = "expr.relative_abundance.abgen.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# Plot abundance

# Define colors for the common taxa in the taxonomic hierarchy
hexcolor.phylum = data.frame(fread("/path-to-script/hexcolor.phylum"))
hexcolor.class = data.frame(fread("/path-to-script/hexcolor.class"))
hexcolor.family = data.frame(fread("/path-to-script/hexcolor.family"))
hexcolor.genus = data.frame(fread("/path-to-script/hexcolor.genus"))

hexcolor.phylum = setNames(as.character(hexcolor.phylum$HEX), hexcolor.phylum$Phylum)
hexcolor.class = setNames(as.character(hexcolor.class$HEX), hexcolor.class$Class)
hexcolor.family = setNames(as.character(hexcolor.family$HEX), hexcolor.family$Family)
hexcolor.genus = setNames(as.character(hexcolor.genus$HEX), hexcolor.genus$Genus)

# Transform to proportions (relative abundances)
ps1.rp = transform_sample_counts(ps1, function(OTU) OTU/sum(OTU))

# Top N taxa
N = 200
topN = names(sort(taxa_sums(ps1), decreasing=TRUE))[1:N]
ps1.topN = transform_sample_counts(ps1, function(OTU) OTU/sum(OTU))
ps1.topN = prune_taxa(topN, ps1.topN)

ps1.topN

pdf("plot_bar.pdf", width = 10, height = 12, pointsize = 12)
plot_bar(ps1.rp, x = "Sample_ID", fill = "Phylum", title = paste(ntaxa(ps1.rp), "Taxa colored by Phylum")) + 
	geom_bar(stat = "identity", size = 0.1, color = "black") + 
	facet_wrap(~ Group, scales = "free_x", nrow = 1) + guides(fill = guide_legend(ncol = 1)) + 
	scale_fill_manual(values = hexcolor.phylum) + scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
	theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10), 
		axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1), 
		axis.text.y = element_text(size = 10)) + xlab("Sample") + ylab("Relative abundance")

plot_bar(ps1.rp, x = "Sample_ID", fill = "Class", title = paste(ntaxa(ps1.rp), "Taxa colored by Class")) + 
	geom_bar(stat = "identity", size = 0.1, color = "black") + 
	facet_wrap(~ Group, scales = "free_x", nrow = 1) + guides(fill = guide_legend(ncol = 1)) + 
	scale_fill_manual(values = hexcolor.class) + scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
	theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10), 
		axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1), 
		axis.text.y = element_text(size = 10)) + xlab("Sample") + ylab("Relative abundance")

plot_bar(ps1.topN, x = "Sample_ID", fill = "Phylum", title = paste("Top",N, "Taxa colored by Phylum")) + 
	geom_bar(stat = "identity", size = 0.1, color = "black") + 
	facet_wrap(~ Group, scales = "free_x", nrow = 1) + guides(fill = guide_legend(ncol = 1)) + 
	scale_fill_manual(values = hexcolor.phylum) + scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
	theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10), 
		axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1), 
		axis.text.y = element_text(size = 10)) + xlab("Sample") + ylab("Relative abundance")

plot_bar(ps1.topN, x = "Sample_ID", fill = "Class", title = paste("Top",N, "Taxa colored by Class")) + 
	geom_bar(stat = "identity", size = 0.1, color = "black") + 
	facet_wrap(~ Group, scales = "free_x", nrow = 1) + guides(fill = guide_legend(ncol = 1)) + 
	scale_fill_manual(values = hexcolor.class) + scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
	theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10), 
		axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1), 
		axis.text.y = element_text(size = 10)) + xlab("Sample") + ylab("Relative abundance")

plot_bar(ps1.topN, x = "Sample_ID", fill = "Family", title = paste("Top",N, "Taxa colored by Family")) + 
	geom_bar(stat = "identity", size = 0.1, color = "black") + 
	facet_wrap(~ Group, scales = "free_x", nrow = 1) + guides(fill = guide_legend(ncol = 1)) + 
	scale_fill_manual(values = hexcolor.family) + scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
	theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10), 
		axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1), 
		axis.text.y = element_text(size = 10)) + xlab("Sample") + ylab("Relative abundance")

plot_bar(ps1.topN, x = "Sample_ID", fill = "Genus", title = paste("Top",N, "Taxa colored by Genus")) + 
	geom_bar(stat = "identity", size = 0.1, color = "black") + 
	facet_wrap(~ Group, scales = "free_x", nrow = 1) + guides(fill = guide_legend(ncol = 1)) + 
	scale_fill_manual(values = hexcolor.genus) + scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
	theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10), 
		axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1), 
		axis.text.y = element_text(size = 10)) + xlab("Sample") + ylab("Relative abundance")
dev.off()

# alpha diversity

# Select alpha-diversity measures
betaM = c("Observed", "Chao1", "Shannon", "Simpson")

png("plot_richness-sample.png", width = 10, height = 5, units = "in", res = 300)
plot_richness(ps1, x = "Sample_ID", measures = betaM, color = "Group", nrow = 1) + 
	geom_point(size = 2) + theme_bw() + theme(legend.position = "top", axis.title.x = element_blank(), 
		axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1), 
		axis.text.y = element_text(size = 10)) + labs(x = "Sample ID", y = "Alpha Diversity Measure")
dev.off()

png("plot_richness-group.png", width = 6, height = 5, units = "in", res = 300)
plot_richness(ps1, x = "Group", measures = betaM, color = "Group", nrow = 1) + 
	geom_point(size = 2) + theme_bw() + theme(legend.position = "none", axis.title.x = element_blank(), 
		axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1), 
		axis.text.y = element_text(size = 10)) + labs(x = "Group", y = "Alpha Diversity Measure")
dev.off()

# Or, plot with ggplot function
ad = estimate_richness(ps1, measures = betaM)
ad = merge(data.frame(sample_data(ps1)), ad, by = "row.names")
ad = reshape2::melt(ad[,c("Sample_ID", "Group", betaM)], id = c("Sample_ID", "Group"))

png("plot_richness.boxplot.png", width = 6, height = 5, units = "in", res = 300)
ggplot(ad, aes(Group, value, color = Group)) + geom_boxplot(outlier.shape = NA, size = 0.8, width = 0.8) + 
	geom_quasirandom(size = 1, color = "black") + facet_wrap(~ variable, scales = "free_y", nrow = 1) + 
	theme_bw() + theme(legend.position = "none", 
		axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1), 
		axis.text.y = element_text(size = 10)) + labs(x = "Group", y = "Alpha Diversity Measure")
dev.off()

# beta diversity
unifracs = GUniFrac::GUniFrac(otu_table(ps1), phy_tree(ps1), alpha = c(0, 0.5, 1))$unifracs

# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm = TRUE){
	exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

# Perform variance stabilizing transformation
ps1.vsd = ps1
dds = phyloseq_to_deseq2(ps1, ~ Group + Batch)
geoMeans = apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds, geoMeans = geoMeans)
dds = estimateDispersions(dds)
abund.vsd = getVarianceStabilizedData(dds)
abund.vsd[abund.vsd < 0] = 0	# set negative values to 0
otu_table(ps1.vsd) = otu_table(abund.vsd, taxa_are_rows = TRUE)

# Create distance objects
dist_un = as.dist(unifracs[, , "d_UW"])       # Unweighted UniFrac
attr(dist_un, "method") = "Unweighted UniFrac"

dist_wu = as.dist(unifracs[, , "d_1"])        # Weighted UniFrac
attr(dist_wu, "method") = "Weighted UniFrac"

dist_vu = as.dist(unifracs[, , "d_VAW"])      # Variance-adjusted-weighted UniFrac
attr(dist_vu, "method") = "Variance-adjusted-weighted UniFrac"

dist_gu = as.dist(unifracs[, , "d_0.5"])      # GUniFrac with alpha 0.5
attr(dist_gu, "method") = "GUniFrac with alpha 0.5"

dist_bc = phyloseq::distance(ps1.vsd, "bray") # Bray-Curtis

# Perform ordination
ord_un = ordinate(ps1, method = "PCoA", distance = dist_un)
ord_wu = ordinate(ps1, method = "PCoA", distance = dist_wu)
ord_vu = ordinate(ps1, method = "PCoA", distance = dist_vu)
ord_gu = ordinate(ps1, method = "PCoA", distance = dist_gu)
set.seed(12345)
ord_bc = ordinate(ps1, method = "NMDS", distance = dist_bc)

pdf("plot_ordination.pdf", width = 6, height = 5, pointsize = 12)
plot_ordination(ps1, ord_un, type = "samples", color = "Group", title = "PCoA on unweighted UniFrac") + 
	theme_bw() + coord_fixed(ratio = 1) + geom_text_repel(aes(label = Sample_ID), size = 2, segment.size = 0.2) +
	stat_ellipse(aes(group = Group), type = "t", size = 0.5, linetype = "dashed")

plot_ordination(ps1, ord_wu, type = "samples", color = "Group", title = "PCoA on weighted UniFrac") + 
	theme_bw() + coord_fixed(ratio = 1) + geom_text_repel(aes(label = Sample_ID), size = 2, segment.size = 0.2) +
	stat_ellipse(aes(group = Group), type = "t", size = 0.5, linetype = "dashed")

plot_ordination(ps1, ord_vu, type = "samples", color = "Group", title = "PCoA on Variance adjusted weighted UniFrac") + 
	theme_bw() + coord_fixed(ratio = 1) + geom_text_repel(aes(label = Sample_ID), size = 2, segment.size = 0.2) +
	stat_ellipse(aes(group = Group), type = "t", size = 0.5, linetype = "dashed")

plot_ordination(ps1, ord_gu, type = "samples", color = "Group", title = "PCoA on GUniFrac with alpha 0.5") + 
	theme_bw() + coord_fixed(ratio = 1) + geom_text_repel(aes(label = Sample_ID), size = 2, segment.size = 0.2) +
	stat_ellipse(aes(group = Group), type = "t", size = 0.5, linetype = "dashed")

plot_ordination(ps1, ord_bc, type = "samples", color = "Group", title = "NMDS on Bray-Curtis distance (VST)") + 
	theme_bw() + coord_fixed(ratio = 1) + geom_text_repel(aes(label = Sample_ID), size = 2, segment.size = 0.2) +
	stat_ellipse(aes(group = Group), type = "t", size = 0.5, linetype = "dashed")
dev.off()

# Permanova test using the adonis function
# adonis tests are significant, reject the null hypothesis that our conditions have the same centroid.
adonis(dist_un ~ Group, data.frame(sample_data(ps1)), permutations = 1000)

adonis(dist_wu ~ Group, data.frame(sample_data(ps1)), permutations = 1000)

adonis(dist_vu ~ Group, data.frame(sample_data(ps1)), permutations = 1000)

adonis(dist_gu ~ Group, data.frame(sample_data(ps1)), permutations = 1000)

adonis(dist_bc ~ Group, data.frame(sample_data(ps1)), permutations = 1000)

# Homogeneity of dispersion test
# betadisper results are not significant, meaning we cannot reject the null hypothesis that our groups have the same dispersions (i.e. having same dispersions)
# We can be more confident that our adonis result are real results, and not due to differences in group dispersions
disp1 = betadisper(dist_un, group = data.frame(sample_data(ps1))$Group)
permutest(disp1, permutations = 1000)

disp2 = betadisper(dist_wu, group = data.frame(sample_data(ps1))$Group)
permutest(disp2, permutations = 1000)

disp3 = betadisper(dist_vu, group = data.frame(sample_data(ps1))$Group)
permutest(disp3, permutations = 1000)

disp4 = betadisper(dist_gu, group = data.frame(sample_data(ps1))$Group)
permutest(disp4, permutations = 1000)

disp5 = betadisper(dist_bc, group = data.frame(sample_data(ps1))$Group)
permutest(disp5, permutations = 1000)

# Save current workspace
# save.image(file="image3.RData")
