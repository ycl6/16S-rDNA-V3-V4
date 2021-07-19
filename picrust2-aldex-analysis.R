#

# Run picrust2

#######################################################	
# Run below commands outside R in Linux environment
# If you have picrust2 installed on a conda environment, activate the required environment, e.g.
# $ conda activate picrust2
#
# Requires fasta file "expr.asv.fasta" and biom file "expr.biom" generated during phyloseq analysis
# The `--output` argument to specify the output folder for final files
# The `--processes N` argument to specify the number of CPUs to run picrust2 in parallel
#
#----- START: Run in terminal/console -----#
# Run picrust2_pipeline.py
/path-to-picrust2_pipeline.py --study_fasta expr.asv.fasta --input expr.biom \
	--output picrust2_out_stratified --processes 2 --stratified --remove_intermediate --verbose

# Locate and take note the directory that keeps the picrust2 mapfiles
locate description_mapfiles
#----- END: Run in terminal/console -----#
#
# Example location:
# /home/user/miniconda3/envs/picrust2/lib/python3.6/site-packages/picrust2/default_files/description_mapfiles
# /home/user/miniconda3/envs/picrust2/lib/python3.6/site-packages/picrust2/default_files/description_mapfiles/KEGG_modules_info.tsv.gz
# /home/user/miniconda3/envs/picrust2/lib/python3.6/site-packages/picrust2/default_files/description_mapfiles/KEGG_pathways_info.tsv.gz
# ...
# Deactivate the environment if required
# $ conda deactivate
#######################################################	

library("data.table")	# Also requires R.utils to read gz and bz2 files
library("phyloseq")
library("ALDEx2")
library("tidyverse")
options(width=190)

# Continued from phyloseq workflow (phyloseq-analysis.R)

# Load saved workspace
# load("image3.RData")

# Note to change the PATH to the "picrust2" output folder accordingly
# Full path is not necessary if `R` is executed in the directory one-level above "picrust2_out_stratified"
picrust2 = "picrust2_out_stratified"

# Note to change the PATH to the picrust2 mapfiles accordingly
# Use the "description_mapfiles" PATH you located with the `locate` command above
mapfile = "/home/user/miniconda3/envs/picrust2/lib/python3.6/site-packages/picrust2/default_files/description_mapfiles"

# Set picrust2 output file paths
p2_EC = paste0(picrust2, "/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz")
p2_KO = paste0(picrust2, "/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz")
p2_PW = paste0(picrust2, "/pathways_out/path_abun_unstrat.tsv.gz")

# Set picrust2 map file paths
mapfile_EC = paste0(mapfile, "/ec_level4_info.tsv.gz")
mapfile_KO = paste0(mapfile, "/ko_info.tsv.gz")
mapfile_PW = paste0(mapfile, "/metacyc_pathways_info.txt.gz")

# Load map files
mapEC = as.data.frame(fread(mapfile_EC, header = FALSE))
colnames(mapEC) = c("function","description")
mapKO = as.data.frame(fread(mapfile_KO, header = FALSE, sep = "\t"))
colnames(mapKO) = c("function","description")
mapPW = as.data.frame(fread(mapfile_PW, header = FALSE))
colnames(mapPW) = c("pathway","description")

# Load picrust2 output files
p2EC = as.data.frame(fread(p2_EC))
rownames(p2EC) = p2EC$"function"
p2EC = as.matrix(p2EC[,-1])
p2EC = round(p2EC)

p2KO = as.data.frame(fread(p2_KO))
rownames(p2KO) = p2KO$"function"
p2KO = as.matrix(p2KO[,-1])
p2KO = round(p2KO)

p2PW = as.data.frame(fread(p2_PW))
rownames(p2PW) = p2PW$"pathway"
p2PW = as.matrix(p2PW[,-1])
p2PW = round(p2PW)

# Perform compositional data analysis (CoDA)
# Allow statistical analysis between 2 groups/conditions
# The default "mc.samples" option is 128
set.seed(12345)
aldex2_EC = aldex(p2EC, sample_data(ps1)$Group, mc.samples = 500, test = "t", effect = TRUE, denom = "iqlr", verbose = TRUE)
head(aldex2_EC)

set.seed(12345)
aldex2_KO = aldex(p2KO, sample_data(ps1)$Group, mc.samples = 500, test = "t", effect = TRUE, denom = "iqlr", verbose = TRUE)
head(aldex2_KO)

set.seed(12345)
aldex2_PW = aldex(p2PW, sample_data(ps1)$Group, mc.samples = 500, test = "t", effect = TRUE, denom = "iqlr", verbose = TRUE)
head(aldex2_PW)

# Check estimated effect size range
# ALDEx2 authors suggest that an effect size of 1 or greater can be used as significance cutoff
png("images/ALDEx2_picrust2_effect.png", width = 6, height = 6, units = "in", res = 300)
par(mfrow = c(2,2))
hist(aldex2_EC$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "EC")
hist(aldex2_KO$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "KO")
hist(aldex2_PW$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "Pathway")
dev.off()

# Create MW (fold-change to variance/effect) and MA (Bland-Altman) plots
# The "test" option can be "welch" or "wilcox"
png("ALDEx2_picrust2_MW_MA.png", width = 6, height = 8, units = "in", res = 300)
par(mfrow = c(3,2))
aldex.plot(aldex2_EC, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
title(main = "(EC) MW Plot")

aldex.plot(aldex2_EC, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, called.cex = 0.6, cutoff = 0.05, xlab = "Log-ratio abundance", ylab = "Difference")
title(main = "(EC) MA Plot")

aldex.plot(aldex2_KO, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
title(main = "(KO) MW Plot")

aldex.plot(aldex2_KO, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, called.cex = 0.6, cutoff = 0.05, xlab = "Relative abundance", ylab = "Difference")
title(main = "(KO) MA Plot")

aldex.plot(aldex2_PW, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
title(main = "(PW) MW Plot")

aldex.plot(aldex2_PW, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, called.cex = 0.6, cutoff = 0.05, xlab = "Relative abundance", ylab = "Difference")
title(main = "(PW) MA Plot")
dev.off()

# Observe relationship between effect, difference, and P values
png("ALDEx2_picrust2_P_adjP.png", width = 6, height = 8, units = "in", res = 300)
par(mfrow = c(3,2))
plot(aldex2_EC$effect, aldex2_EC$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, xlab = "Effect size", ylab = "P value", main = "(EC) Effect size plot")
points(aldex2_EC$effect, aldex2_EC$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_EC$diff.btw, aldex2_EC$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, xlab = "Difference", ylab = "P value", main = "(EC) Volcano plot")
points(aldex2_EC$diff.btw, aldex2_EC$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")

plot(aldex2_KO$effect, aldex2_KO$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, xlab = "Effect size", ylab = "P value", main = "(KO) Effect size plot")
points(aldex2_KO$effect, aldex2_KO$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_KO$diff.btw, aldex2_KO$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, xlab = "Difference", ylab = "P value", main = "(KO) Volcano plot")
points(aldex2_KO$diff.btw, aldex2_KO$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")

plot(aldex2_PW$effect, aldex2_PW$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, xlab = "Effect size", ylab = "P value", main = "(PW) Effect size plot")
points(aldex2_PW$effect, aldex2_PW$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_PW$diff.btw, aldex2_PW$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, xlab = "Difference", ylab = "P value", main = "(PW) Volcano plot")
points(aldex2_PW$diff.btw, aldex2_PW$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
dev.off()

# Merge with map file data
df_EC = aldex2_EC %>% tibble::rownames_to_column(var = "EC") %>% inner_join(mapEC, by = c("EC" = "function")) %>% arrange(EC)
df_KO = aldex2_KO %>% tibble::rownames_to_column(var = "KO") %>% inner_join(mapKO, by = c("KO" = "function")) %>% arrange(KO)
df_PW = aldex2_PW %>% tibble::rownames_to_column(var = "Pathway") %>% inner_join(mapPW, by = c("Pathway" = "pathway")) %>% arrange(Pathway)

# Output to file
write.table(df_EC, file = "ALDEx2_picrust2_EC_results.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(df_KO, file = "ALDEx2_picrust2_KO_results.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(df_PW, file = "ALDEx2_picrust2_Pathway_results.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# Save current workspace
# save.image(file = "image4.RData")
