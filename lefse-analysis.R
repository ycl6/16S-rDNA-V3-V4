#

library("data.table")
library("phyloseq")
library("ggplot2")
library("grid")
library("tidyverse")
Sys.setlocale("LC_COLLATE","C")
options(width=190)

# Continued from phyloseq workflow (phyloseq-analysis.R)

# Load saved workspace
# load("image3.RData")

# Note to change the PATH to the "lefse.R" script accordingly
source("/path-to-script/lefse.R", local = TRUE)

# Create directory if not exist
if(!dir.exists("lefse")) { dir.create("lefse") }

# Prepare input data
# Add 'subject = TRUE' to include subject in the data.frame
tax1 = lefse_1name_obj(ps1, sample_data(ps1)$Group)
lefse1 = lefse_obj(ps1)
lefse1 = rbind(tax1, lefse1)

# Replace unsupported chars with underscore
lefse1$name = gsub(" ","_",lefse1$name)
lefse1$name = gsub("-","_",lefse1$name)
lefse1$name = gsub("/","_",lefse1$name)

# For **Silva v132** users, apply below to fix identifical phylum & class names (Actinobacteria and Deferribacteres)
lefse1 = fix_taxa_silva132(lefse1)

# Output the prepared LEfSe input to file
write.table(lefse1, file="lefse/expr1.lefse_table.txt", sep = "\t", quote = F, row.names = F, col.names = F)

# Set the full-paths to LEfSe’s python scripts if they are not in PATH
lefse-format_input = "/path-to-lefse-format_input.py"	# e.g. /home/ngs/lefse/lefse-format_input.py
run_lefse = "/path-to-run_lefse.py"			# e.g. /home/ngs/lefse/run_lefse.py

# Run LEfSe

system2(lefse_format_input,
	args = c("lefse/expr1.lefse_table.txt", "lefse/expr1.lefse_table.in", "-c 1", "-o 1000000"))


# set Kruskal-Wallis alpha (-a) to 1 to allow returning of all P-value to perform adjustment later
system2(run_lefse,
	args = c("lefse/expr1.lefse_table.in", "lefse/expr1.lefse_table.res", "-b 100", "-a 1", "-l 1"), stdout = TRUE)

# Perform multiple testing correction on KW P-values
q = 0.1	# fdr threshold at 0.1

expr1 = data.frame(fread("lefse/expr1.lefse_table.res"))
expr1 = pcorrection(expr1, q)

# You can use table() to check number of taxa pass the threshold
table(expr1$V3 != "")

# Write new result file with the P-value column replaced by the FDR
write.table(expr1[,c(1:4,6)], file = "lefse/expr1.lefse_table.res.padj", sep = "\t", quote = F, row.names = F, col.names = F)

# Plot cladogram

#######################################################
# Run below commands outside R in Linux environment
# Both export2graphlan & graphlan require Python 2.7
# If you have export2graphlan & graphlan installed on a conda environment, activate the required environment, e.g.
# $ conda activate graphlan
#
# Requires "lefse/expr1.lefse_table.txt" and "lefse/expr1.lefse_table.res.padj" generated above
# Requries "lefse/expr1.colors" with color setting
# Run export2graphlan & graphlan with full-paths if they are not in PATH
#
#----- START: Run in terminal/console -----#
# Run export2graphlan
/path-to-export2graphlan.py -i lefse/expr1.lefse_table.txt -o lefse/expr1.lefse_table.res.padj \
        -t lefse/expr1.graphlan_tree.txt -a lefse/expr1.graphlan_annot.txt --external_annotations 2,3,4,5,6 \
        --fname_row 0 --biomarkers2colors lefse/expr1.colors

# Run graphlan
/path-to-graphlan_annotate.py --annot lefse/expr1.graphlan_annot.txt lefse/expr1.graphlan_tree.txt \
	lefse/expr1.graphlan_outtree.txt

# Convert 'lsqb' and 'rsqb' back to square bracket symbols
sed 's/lsqb/[/' lefse/expr1.graphlan_outtree.txt | sed 's/rsqb/]/' > lefse/expr1.graphlan_outtree_sqb.txt
sed 's/lsqb/[/' lefse/expr1.lefse_table.res.padj | sed 's/rsqb/]/' > lefse/expr1.lefse_table.res_sqb.padj

# Create cladogram
/path-to-graphlan.py --dpi 150 lefse/expr1.graphlan_outtree_sqb.txt lefse/expr1.graphlan.png \
	--external_legends --size 8 --pad 0.2

# Convert to readble output
perl /path-to-lefse.pl lefse/expr1.lefse_table.res_sqb.padj lefse/expr1.graphlan_outtree_sqb.txt lefse/expr1.out
#----- END: Run in terminal/console -----#
#
# Deactivate the environment if required
# $ conda deactivate
#######################################################

# Plot results

res1 = data.frame(fread("lefse/expr1.out"))
res1 = res1[order(res1$order),]
names(res1)[c(5:7)] = c("Group","LDA","FDR")
res1$taxon = paste0(res1$taxon, "(",res1$rank, ";", " ",res1$order, ")")
res1$taxon = as.factor(res1$taxon)
res1$taxon = factor(res1$taxon, levels = unique(res1$taxon[order(res1$order, decreasing = T)]))
res1$Group = as.factor(res1$Group)
levels(res1$taxon) = gsub("Escherichia_Shigella","Escherichia/Shigella",levels(res1$taxon))
levels(res1$taxon) = gsub("_UCG_","_UCG-",levels(res1$taxon))

png("lefse/expr1.lefse_table.png", width = 8, height = 8, units = "in", res = 300)
ggplot(res1, aes(taxon, LDA, fill = Group)) + geom_bar(stat = "identity", width = 0.7, size = 0.5) +
	coord_flip() + theme_bw() + facet_grid(rows = vars(Group), scales = "free_y", space = "free_y") +
	scale_fill_manual(values = c("control" = "blue", "patient" = "red")) +
	theme(legend.position = "none",
	      axis.text.y = element_text(face = "bold", size = 8),
	      strip.placement = "outside",
	      strip.text.y = element_text(face = "bold", angle = 0)) +
	labs(title = "Linear discriminant analysis (LDA)\nEffect Size Analysis", x = "Taxon", y = "LDA")
dev.off()
