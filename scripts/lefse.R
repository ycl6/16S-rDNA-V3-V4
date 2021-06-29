# By I-Hsuan Lin
# Return lefse objects

# Use subject = TRUE to include subject in output
lefse_1name_obj <- function(ps, group, subject = FALSE) {
	tax = as.data.frame(rbind(c("class",as.vector(group))), stringsAsFactors = FALSE)
	colnames(tax) = c("name", sample_names(ps))

	if(subject) {
		tax = rbind(tax, c("subject", sample_names(ps)))
	}

	return(tax)
}

# Use subject = TRUE to include subject in output
lefse_2name_obj <- function(ps, group, subgroup, subject = FALSE) {
	tax = as.data.frame(rbind(c("class",as.vector(group)), c("subclass", as.vector(subgroup))), stringsAsFactors = FALSE)
	colnames(tax) = c("name", sample_names(ps))

	if(subject) {
		tax = rbind(tax, c("subject", sample_names(ps)))
	}

	return(tax)
}

lefse_obj <- function(ps) {
	tax1 = ps %>% tax_glom(taxrank = "Kingdom") %>% transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt()
	tax1 = tax1[order(tax1$Sample),]
	tax1 = reshape2::dcast(tax1, Kingdom~Sample, value.var = "Abundance")
	names(tax1)[1] = "name"

	tax2 = ps %>% tax_glom(taxrank = "Phylum") %>% transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt()
	tax2$name = paste(tax2$Kingdom, tax2$Phylum, sep = "|")
	tax2 = tax2[order(tax2$Sample, tax2$name),]
	tax2 = reshape2::dcast(tax2, name~Sample, value.var = "Abundance")

	tax3 = ps %>% tax_glom(taxrank = "Class") %>% transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt()
	tax3$Class = as.character(tax3$Class)
	tax3[grep("^Subgroup_",tax3$Class),]$Class = paste(tax3[grep("^Subgroup_",tax3$Class),]$Phylum,tax3[grep("^Subgroup_",tax3$Class),]$Class, sep="_")
	tax3$Class = as.factor(tax3$Class)
	tax3$name = paste(tax3$Kingdom, tax3$Phylum, tax3$Class, sep = "|")
	tax3 = tax3[order(tax3$Sample, tax3$name),]
	tax3 = reshape2::dcast(tax3, name~Sample, value.var = "Abundance")

	tax4 = ps %>% tax_glom(taxrank = "Order") %>% transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt()
	tax4$Class = as.character(tax4$Class)
	tax4[grep("^Subgroup_",tax4$Class),]$Class = paste(tax4[grep("^Subgroup_",tax4$Class),]$Phylum,tax4[grep("^Subgroup_",tax4$Class),]$Class, sep="_")
	tax4$Class = as.factor(tax4$Class)
	tax4$Order = as.character(tax4$Order)
	tax4[grep("^Subgroup_",tax4$Order),]$Order = paste(tax4[grep("^Subgroup_",tax4$Order),]$Phylum,tax4[grep("^Subgroup_",tax4$Order),]$Order, sep="_")
	tax4$Order = as.factor(tax4$Order)
	tax4$name = paste(tax4$Kingdom, tax4$Phylum, tax4$Class, tax4$Order, sep = "|")
	tax4 = tax4[order(tax4$Sample, tax4$name),]
	tax4 = reshape2::dcast(tax4, name~Sample, value.var = "Abundance")

	tax5 = ps %>% tax_glom(taxrank = "Family") %>% transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt()
	tax5$Class = as.character(tax5$Class)
	tax5[grep("^Subgroup_",tax5$Class),]$Class = paste(tax5[grep("^Subgroup_",tax5$Class),]$Phylum,tax5[grep("^Subgroup_",tax5$Class),]$Class, sep="_")
	tax5$Class = as.factor(tax5$Class)
	tax5$Order = as.character(tax5$Order)
	tax5[grep("^Subgroup_",tax5$Order),]$Order = paste(tax5[grep("^Subgroup_",tax5$Order),]$Phylum,tax5[grep("^Subgroup_",tax5$Order),]$Order, sep="_")
	tax5$Order = as.factor(tax5$Order)
	tax5$Family = as.character(tax5$Family)
	tax5[grep("^Subgroup_",tax5$Family),]$Family = paste(tax5[grep("^Subgroup_",tax5$Family),]$Phylum,tax5[grep("^Subgroup_",tax5$Family),]$Family, sep="_")
	tax5[grep("^UCG-",tax5$Family),]$Family = paste(tax5[grep("^UCG-",tax5$Family),]$Order,tax5[grep("^UCG-",tax5$Family),]$Family, sep="_")
	tax5[grep("^uncultured",tax5$Family),]$Family = paste(tax5[grep("^uncultured",tax5$Family),]$Order,tax5[grep("^uncultured",tax5$Family),]$Family, sep="_")
	tax5$Family = as.factor(tax5$Family)
	tax5$name = paste(tax5$Kingdom, tax5$Phylum, tax5$Class, tax5$Order, tax5$Family, sep = "|")
	tax5 = tax5[order(tax5$Sample, tax5$name),]
	tax5 = reshape2::dcast(tax5, name~Sample, value.var = "Abundance")

	tax6 = ps %>% tax_glom(taxrank = "Genus") %>% transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt()
	tax6$Class = as.character(tax6$Class)
	tax6[grep("^Subgroup_",tax6$Class),]$Class = paste(tax6[grep("^Subgroup_",tax6$Class),]$Phylum,tax6[grep("^Subgroup_",tax6$Class),]$Class, sep="_")
	tax6$Class = as.factor(tax6$Class)
	tax6$Order = as.character(tax6$Order)
	tax6[grep("^Subgroup_",tax6$Order),]$Order = paste(tax6[grep("^Subgroup_",tax6$Order),]$Phylum,tax6[grep("^Subgroup_",tax6$Order),]$Order, sep="_")
	tax6$Order = as.factor(tax6$Order)
	tax6$Family = as.character(tax6$Family)
	tax6[grep("^Subgroup_",tax6$Family),]$Family = paste(tax6[grep("^Subgroup_",tax6$Family),]$Phylum,tax6[grep("^Subgroup_",tax6$Family),]$Family, sep="_")
	tax6[grep("^UCG-",tax6$Family),]$Family = paste(tax6[grep("^UCG-",tax6$Family),]$Phylum,tax6[grep("^UCG-",tax6$Family),]$Family, sep="_")
	tax6[grep("^uncultured",tax6$Family),]$Family = paste(tax6[grep("^uncultured",tax6$Family),]$Order,tax6[grep("^uncultured",tax6$Family),]$Family, sep="_")
	tax6$Family = as.factor(tax6$Family)
	tax6$Genus = as.character(tax6$Genus)
	tax6[grep("^UCG-",tax6$Genus),]$Genus = paste(tax6[grep("^UCG-",tax6$Genus),]$Family,tax6[grep("^UCG-",tax6$Genus),]$Genus, sep="_")
	tax6[grep("^uncultured",tax6$Genus),]$Genus = paste(tax6[grep("^uncultured",tax6$Genus),]$Family,tax6[grep("^uncultured",tax6$Genus),]$Genus, sep="_")
	tax6$Genus = gsub("_uncultured_uncultured","_uncultured",tax6$Genus)
	tax6$Genus = as.factor(tax6$Genus)
	tax6$name = paste(tax6$Kingdom, tax6$Phylum, tax6$Class, tax6$Order, tax6$Family, tax6$Genus, sep = "|")
	tax6 = tax6[order(tax6$Sample, tax6$name),]
	tax6 = reshape2::dcast(tax6, name~Sample, value.var = "Abundance")

	lefse = rbind(tax1, tax2, tax3, tax4, tax5, tax6)
	lefse = lefse[order(lefse$name),]
	lefse$name = as.character(lefse$name)

	return(lefse)
}

# To fix identifical phylum & class names (Actinobacteria and Deferribacteres) for Silva v132 users
fix_taxa_silva132 <- function(DF) {
	# fix Actinobacteria
	DF$name = gsub("Bacteria\\|Actinobacteria","Bacteria|Actinobacteria_P", DF$name)
	DF$name = gsub("Bacteria\\|Actinobacteria_P\\|Actinobacteria","Bacteria|Actinobacteria_P|Actinobacteria_C", DF$name)
	# fix Deferribacteres
	DF$name = gsub("Bacteria\\|Deferribacteres","Bacteria|Deferribacteres_P", DF$name)
	DF$name = gsub("Bacteria\\|Deferribacteres_P\\|Deferribacteres","Bacteria|Deferribacteres_P|Deferribacteres_C", DF$name)

	return(DF)
}

pcorrection <- function(res, q) {
	res$q = p.adjust(res$V5, method = "fdr")
	res[is.na(res$q),]$q = "-"

	res[res$q > q, "V3"] = ""
	res[res$q > q, "V4"] = ""
	res[is.na(res$V4), "V4"] = ""
	res[res$q > q, "V5"] = "-"
	res[res$q > q,  "q"] = "-"
	
	return(res)
}
