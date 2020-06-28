#

library("data.table")
library("phyloseq")
library("ggplot2")
library("dplyr")  #"%>%"
library("ggbeeswarm")
library("ggrepel")
options(width=190)

# Continued from phyloseq workflow (phyloseq-analysis.R)

# lefse
# Update path to "lefse.R" if necessary
source("/path-to-script/lefse.R", local = TRUE)

# Create directory if not exist
if(!dir.exists("lefse")) { dir.create("lefse") }

tax1 = lefse_1name_obj(ps1, sample_data(ps1)$Group)
lefse1 = lefse_obj(ps1)
lefse1 = rbind(tax1, lefse1)
lefse1$name = gsub("-","_",lefse1$name)
lefse1$name = gsub("/","_",lefse1$name)
write.table(lefse1, file="lefse/expr1.lefse_table.txt", sep = "\t", quote = F, row.names = F, col.names = F)

# Run LEfSe
# Update path to LEfSe python scripts if necessary
nsegata_lefse = "/path-to-script/nsegata-lefse-9adc3a62460e/"
system2(paste0(nsegata_lefse, "format_input.py"), args = c("lefse/expr1.lefse_table.txt", "lefse/expr1.lefse_table.in", "-c 1", "-u 2", "-o 1000000"))
system2(paste0(nsegata_lefse, "run_lefse.py"), args = c("lefse/expr1.lefse_table.in", "lefse/expr1.lefse_table.res", "-b 100", "-a 1", "-l 2"))	# set KW alpha to 1 to allow returning of all P-value to perform adjustment

# Perform multiple testing correction on KW P-values
q = 0.1	# fdr threshold at 0.1

res0 = data.frame(fread("lefse/expr1.lefse_table.res"))
res0$q = p.adjust(res0$V5, method = "fdr")
res0[is.na(res0$q),]$q = "-"

# You may encounter error message if there are no matching rows to replace, which is alright
res0[res0$q > q,]$V3 = ""
res0[res0$q > q,]$V4 = ""
res0[is.na(res0$V4),]$V4 = ""
res0[res0$q > q,]$V5 = "-"
res0[res0$q > q,]$q = "-"

# You can use table() to check number of taxa pass the threshold
# table(res0$V3 != "")

# Write new result file with the P-value column replaced by the FDR
write.table(res0[,c(1:4,6)], file = "lefse/expr1.lefse_table.res.padj", sep = "\t", quote = F, row.names = F, col.names = F)

# Plot cladogram
# If export2graphlan and GraPhlAn are not in PATH, you need to specify the full path to the scripts
# Update path to "lefse.pl" if necessary
system2("export2graphlan.py", args = c("-i lefse/expr1.lefse_table.txt", "-o lefse/expr1.lefse_table.res.padj", "-t lefse/expr1.graphlan_tree.txt", "-a lefse/expr1.graphlan_annot.txt", "--external_annotations 2,3,4,5,6", "--fname_row 0", "--skip_rows 1", "--biomarkers2colors lefse/expr1.colors"))
system2("graphlan_annotate.py", args = c("--annot lefse/expr1.graphlan_annot.txt", "lefse/expr1.graphlan_tree.txt", "lefse/expr1.graphlan_outtree.txt"))
system2("graphlan.py", args = c("--dpi 150", "lefse/expr1.graphlan_outtree.txt", "lefse/expr1.graphlan.png", "--external_legends", "--size 8", "--pad 0.2"))
system2("perl", args = c("/path-to-script/lefse.pl", "lefse/expr1.lefse_table.res.padj", "lefse/expr1.graphlan_outtree.txt", "lefse/expr1.out"))

# plot res in R
library(grid)
Sys.setlocale("LC_COLLATE","C")

res1 = data.frame(fread("lefse/expr1.out"))
res1 = res1[order(res1$order),]
names(res1)[c(5:7)] = c("Group","LDA","FDR")
res1$taxon = paste0(res1$taxon, "(",res1$rank, ";", " ",res1$order, ")")
res1$taxon = as.factor(res1$taxon)
res1$taxon = factor(res1$taxon, levels = unique(res1$taxon[order(res1$order, decreasing=T)]))
res1$Group = as.factor(res1$Group)
levels(res1$taxon) = gsub("Escherichia_Shigella","Escherichia/Shigella",levels(res1$taxon))
levels(res1$taxon) = gsub("_UCG_","_UCG-",levels(res1$taxon))

plot_lefse = ggplot(res1, aes(taxon, LDA, fill = Group)) + geom_bar(stat = "identity", width = 0.7, size = 0.5) + coord_flip() + theme_bw() +
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
