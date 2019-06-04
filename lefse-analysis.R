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
source("lefse.R", local = TRUE)

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
# /tools/nsegata-lefse-82605a2ae7b7/format_input.py lefse/expr1.lefse_table.txt lefse/expr1.lefse_table.in -c 1 -u 2 -o 1000000
# /tools/nsegata-lefse-82605a2ae7b7/run_lefse.py lefse/expr1.lefse_table.in lefse/expr1.lefse_table.res -b 100 -a 0.05 -w 0.05 -l 2
#
# create a "lefse.expr.colors" file to define group color (see example file)
#
# /tools/graphlan_commit_6ca8735/export2graphlan/export2graphlan.py -i lefse/expr1.lefse_table.txt -o lefse/lefse.expr1.lefse_table.res -t lefse/expr1.graphlan_tree.txt -a lefse/expr1.graphlan_annot.txt --external_annotations 2,3,4,5,6 --fname_row 0 --skip_rows 1 --biomarkers2colors lefse/lefse.expr.colors
# /tools/graphlan_commit_6ca8735/graphlan_annotate.py --annot lefse/expr1.graphlan_annot.txt lefse/expr1.graphlan_tree.txt lefse/expr1.graphlan_outtree.txt
# /tools/graphlan_commit_6ca8735/graphlan.py --dpi 150 lefse/expr1.graphlan_outtree.txt lefse/expr1.graphlan.png --external_legends --size 8 --pad 0.2
#######################################################
