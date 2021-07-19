# By I-Hsuan Lin
# Return lefse objects

# Use subject = TRUE to include subject
lefse_1name_obj <- function(ps, group, subject = FALSE) {
	tax = as.data.frame(rbind(c("class",as.vector(group))), stringsAsFactors = FALSE)
	colnames(tax) = c("name", sample_names(ps))

	if(subject) {
		tax = rbind(tax, c("subject", sample_names(ps)))
	}
	return(tax)
}

# Use subject = TRUE to include subject
lefse_2name_obj <- function(ps, group, subgroup, subject = FALSE) {
	tax = as.data.frame(rbind(c("class",as.vector(group)), c("subclass", as.vector(subgroup))), stringsAsFactors = FALSE)
	colnames(tax) = c("name", sample_names(ps))

	if(subject) {
		tax = rbind(tax, c("subject", sample_names(ps)))
	}
	return(tax)
}

lefse_obj <- function(ps) {
	if ( !isNamespaceLoaded("tidyverse") ) {
		print("Load package: tidyverse")
		library("tidyverse")
	}

	tax1 = ps %>% tax_glom(taxrank = "Kingdom") %>% transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt() %>% 
		arrange(Sample) %>% select(Kingdom, Sample, Abundance) %>% spread(Sample, Abundance) %>% rename(name = Kingdom)

	tax2 = ps %>% tax_glom(taxrank = "Phylum") %>% transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt() %>% 
		unite("name", Kingdom:Phylum, sep = "|", remove = FALSE) %>% arrange(Sample, name) %>% 
		select(name, Sample, Abundance) %>% spread(Sample, Abundance)

	tax3 = ps %>% tax_glom(taxrank = "Class") %>% transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt()
	tax3[grep("^Subgroup", tax3$Class),]$Class = paste(tax3[grep("^Subgroup", tax3$Class),]$Phylum, 
							   tax3[grep("^Subgroup", tax3$Class),]$Class, sep = "_")
	tax3[tax3$Class == "Incertae Sedis",]$Class = paste(tax3[tax3$Class == "Incertae Sedis",]$Phylum,
							    tax3[tax3$Class == "Incertae Sedis",]$Class, sep = "_")
	tax3 = tax3 %>% unite("name", Kingdom:Class, sep = "|", remove = FALSE) %>% arrange(Sample, name) %>% 
		select(name, Sample, Abundance) %>% spread(Sample, Abundance)

	tax4 = ps %>% tax_glom(taxrank = "Order") %>% transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt()
	tax4[grep("^Subgroup", tax4$Order),]$Order = paste(tax4[grep("^Subgroup", tax4$Order),]$Class, 
							   tax4[grep("^Subgroup", tax4$Order),]$Order, sep = "_")
	tax4[tax4$Class == "Incertae Sedis",]$Class = paste(tax4[tax4$Class == "Incertae Sedis",]$Class,
							    tax4[tax4$Class == "Incertae Sedis",]$Order, sep = "_")
	tax4 = tax4 %>% unite("name", Kingdom:Order, sep = "|", remove = FALSE) %>% arrange(Sample, name) %>% 
		select(name, Sample, Abundance) %>% spread(Sample, Abundance)

	tax5 = ps %>% tax_glom(taxrank = "Family") %>% transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt()
	tax5[grep("^Subgroup", tax5$Family),]$Family = paste(tax5[grep("^Subgroup ", tax5$Family),]$Order, 
							     tax5[grep("^Subgroup", tax5$Family),]$Family, sep = "_")
	tax5[tax5$Family == "Incertae Sedis",]$Family = paste(tax5[tax5$Family == "Incertae Sedis",]$Order, 
							      tax5[tax5$Class == "Incertae Sedis",]$Family, sep = "_")
	tax5[tax5$Family == "Unknown Family",]$Family = paste(tax5[tax5$Family == "Unknown Family",]$Order, 
							      tax5[tax5$Class == "Unknown Family",]$Family, sep = "_")
	tax5[grep("^UCG-", tax5$Family),]$Family = paste(tax5[grep("^UCG-", tax5$Family),]$Order, 
							 tax5[grep("^UCG-", tax5$Family),]$Family, sep = "_")
	tax5[grep("^type", tax5$Family),]$Family = paste(tax5[grep("^type", tax5$Family),]$Order, 
							 tax5[grep("^type", tax5$Family),]$Family, sep = "_")
	tax5[grep("^Clade", tax5$Family),]$Family = paste(tax5[grep("^Clade", tax5$Family),]$Order, 
							  tax5[grep("^Clade", tax5$Family),]$Family, sep = "_")
	tax5$Family = gsub("\\[(.+)\\]","lsqb\\1rsqb", tax5$Family) # replace square brackets (lsqb & rsqb)
	tax5 = tax5 %>% unite("name", Kingdom:Family, sep = "|", remove = FALSE) %>% arrange(Sample, name) %>% 
		select(name, Sample, Abundance) %>% spread(Sample, Abundance)

	tax6 = ps %>% tax_glom(taxrank = "Genus") %>% transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt()
	tax6[grep("^Subgroup", tax6$Genus),]$Genus = paste(tax6[grep("^Subgroup ", tax6$Genus),]$Family, 
							   tax6[grep("^Subgroup", tax6$Genus),]$Genus, sep = "_")
	tax6[tax6$Genus == "Incertae Sedis",]$Genus = paste(tax6[tax6$Genus == "Incertae Sedis",]$Family, 
							    tax6[tax6$Class == "Incertae Sedis",]$Genus, sep = "_")
	tax6[grep("^UCG-", tax6$Genus),]$Genus = paste(tax6[grep("^UCG-", tax6$Genus),]$Family, 
						       tax6[grep("^UCG-", tax6$Genus),]$Genus, sep = "_")
	tax6[grep("^Clade", tax6$Genus),]$Genus = paste(tax6[grep("^Clade", tax6$Genus),]$Order, 
							tax6[grep("^Clade", tax6$Genus),]$Genus, sep = "_")
	tax6$Genus = gsub("\\[(.+)\\]","lsqb\\1rsqb", tax6$Genus) # replace square brackets (lsqb & rsqb)
	tax6 = tax6 %>% unite("name", Kingdom:Genus, sep = "|", remove = FALSE) %>% arrange(Sample, name) %>% 
		select(name, Sample, Abundance) %>% spread(Sample, Abundance)

	lefse = rbind(tax1, tax2, tax3, tax4, tax5, tax6) %>% arrange(name)

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

