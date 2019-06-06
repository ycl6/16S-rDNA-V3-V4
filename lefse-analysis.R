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
source("/path-to-script/lefse.R", local = TRUE)

# Create directory if not exist
if(!dir.exists("lefse")) { dir.create("lefse") }

tax1 = lefse_1name_obj(ps1, sample_data(ps1)$Group)
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
#
# create a "lefse.expr.colors" file to define group color (see example file)
#
# /path-to-script/graphlan_commit_6ca8735/export2graphlan/export2graphlan.py -i lefse/expr1.lefse_table.txt -o lefse/expr1.lefse_table.res -t lefse/expr1.graphlan_tree.txt -a lefse/expr1.graphlan_annot.txt --external_annotations 2,3,4,5,6 --fname_row 0 --skip_rows 1 --biomarkers2colors lefse/lefse.expr.colors
# /path-to-script/graphlan_commit_6ca8735/graphlan_annotate.py --annot lefse/expr1.graphlan_annot.txt lefse/expr1.graphlan_tree.txt lefse/expr1.graphlan_outtree.txt
# /path-to-script/graphlan_commit_6ca8735/graphlan.py --dpi 150 lefse/expr1.graphlan_outtree.txt lefse/expr1.graphlan.png --external_legends --size 8 --pad 0.2
# perl /path-to-script/lefse.pl lefse/expr1.lefse_table.res lefse/expr1.graphlan_outtree.txt lefse/expr1.out
#######################################################

# plot_res in R
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
