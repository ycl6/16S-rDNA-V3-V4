#

library("data.table")
library("phyloseq")
library("ggplot2")
library("dplyr")
library("grid")
Sys.setlocale("LC_COLLATE","C")
options(width=190)

# Continued from phyloseq workflow (phyloseq-analysis.R)

# Load saved workspace
# save.image(file="image3.RData")

# Note to change the PATH to the "lefse.R" script accordingly
source("/path-to-script/lefse.R", local = TRUE)

# Create directory if not exist
if(!dir.exists("lefse")) { dir.create("lefse") }

# Prepare input data
tax1 = lefse_1name_obj(ps1, sample_data(ps1)$Group)
lefse1 = lefse_obj(ps1)
lefse1 = rbind(tax1, lefse1)
lefse1$name = gsub("-","_",lefse1$name)
lefse1$name = gsub("/","_",lefse1$name)

# For **Silva v132** users, apply below to fix identifical phylum & class names (Actinobacteria and Deferribacteres)
lefse1 = fix_taxa_silva132(lefse1)

write.table(lefse1, file="lefse/expr1.lefse_table.txt", sep = "\t", quote = F, row.names = F, col.names = F)

# Set the full-path to the "nsegata-lefse" folder
nsegata_lefse = "/path-to-script/nsegata-lefse-9adc3a62460e/"

# Run LEfSe
system2(paste0(nsegata_lefse, "format_input.py"), 
	args = c("lefse/expr1.lefse_table.txt", "lefse/expr1.lefse_table.in", "-c 1", "-u 2", "-o 1000000"))
# set KW alpha to 1 to allow returning of all P-value to perform adjustment later
system2(paste0(nsegata_lefse, "run_lefse.py"), 
	args = c("lefse/expr1.lefse_table.in", "lefse/expr1.lefse_table.res", "-b 100", "-a 1", "-l 2"))

# Perform multiple testing correction on KW P-values
q = 0.1	# fdr threshold at 0.1

expr1 = data.frame(fread("lefse/expr1.lefse_table.res"))
expr1 = pcorrection(expr1, q)

# You can use table() to check number of taxa pass the threshold
table(expr1$V3 != "")

# Write new result file with the P-value column replaced by the FDR
write.table(expr1[,c(1:4,6)], file = "lefse/expr1.lefse_table.res.padj", sep = "\t", quote = F, row.names = F, col.names = F)

# Set the full-paths to `export2graphlan`, `GraPhlAn` and `lefse.pl` if they are not in PATH
export2graphlan = "/path-to-script/export2graphlan.py"
graphlan_annotat = "/path-to-script/graphlan_annotate.py"
graphlan = "/path-to-script/graphlan.py"
graphlan_parser = "/path-to-script/lefse.pl"

# Plot cladogram
# Update the coloring file "lefse/expr1.colors" if necessary. The colors should be defined in HSV (hue, saturation, value) scale

system2(export2graphlan,
        args = c("-i lefse/expr1.lefse_table.txt", "-o lefse/expr1.lefse_table.res.padj",
                 "-t lefse/expr1.graphlan_tree.txt", "-a lefse/expr1.graphlan_annot.txt",
                 "--external_annotations 2,3,4,5,6", "--fname_row 0", "--skip_rows 1",
                 "--biomarkers2colors lefse/expr1.colors"))

system2(graphlan_annotat,
        args = c("--annot lefse/expr1.graphlan_annot.txt",
                 "lefse/expr1.graphlan_tree.txt",
                 "lefse/expr1.graphlan_outtree.txt"))

system2(graphlan,
        args = c("--dpi 150", "lefse/expr1.graphlan_outtree.txt", "lefse/expr1.graphlan.png",
                 "--external_legends", "--size 8", "--pad 0.2"))

# Convert to readble output
system2("perl", args = c(graphlan_parser, "lefse/expr1.lefse_table.res.padj",
                         "lefse/expr1.graphlan_outtree.txt", "lefse/expr1.out"))

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

plot_lefse = ggplot(res1, aes(taxon, LDA, fill = Group)) + geom_bar(stat = "identity", width = 0.7, size = 0.5) + coord_flip() + theme_bw() +
facet_wrap(~ Group, ncol = 1, scales = "free_y") + theme(legend.position = "right", axis.text.x = element_text(size = 12), axis.text.y = element_text(face = "bold", size = 8),
strip.text.x = element_text(face = "bold", size = 12)) + labs(title = "Linear discriminant analysis (LDA)\nEffect Size Analysis", x = "Taxon", y = "LDA")

groupN = res1 %>% group_by(Group) %>% summarise(count = length(unique(taxon)))
gt = ggplotGrob(plot_lefse)
panelI.1 = gt$layout$t[grepl("panel", gt$layout$name)]
gt$heights[panelI.1] = unit(groupN$count, "null")
dev.off()

png("lefse/expr1.lefse_table.png", width = 8, height = 8, units = "in", res = 300)
grid.draw(gt)
dev.off()
