#

library("dada2")
library("ggplot2")
options(width = 190)

# Set path
trimmed = "trimmed"	# cutadapt trimmed fastq files
filt = "filt"		# dada2 trimmed fastq files

fns = sort(list.files(trimmed, full.names = TRUE))
fnFs = fns[grep("1.fastq.gz", fns)]
fnRs = fns[grep("2.fastq.gz", fns)]
sample.names = gsub(".1.fastq.gz", "", basename(fnFs))

# Plot quality profile of fastq files
ii = 1:length(sample.names)
pdf("plotQualityProfile.pdf", width = 8, height = 8, pointsize = 12)
for(i in ii) {
	message(paste0("[", i ,"/", length(sample.names), "] ", sample.names[i]))
	print(plotQualityProfile(fnFs[i]) + ggtitle("Fwd"))
	print(plotQualityProfile(fnRs[i]) + ggtitle("Rev"))
}
dev.off()

# Set paths to the dada2-filterd files
filtFs = file.path(filt, basename(fnFs))
filtRs = file.path(filt, basename(fnRs))

# Perform filtering and trimming
# Review "plotQualityProfile.pdf" to select the best paramters for truncLen
out = filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
	# Need to keep paramters consistent between runs of the same study
	truncLen = c(260,200), minLen = 200, maxN = 0, truncQ = 2, maxEE = c(2,5),
	rm.phix = TRUE, compress = TRUE, verbose = TRUE, multithread = TRUE)

out = as.data.frame(out)
rownames(out) = sample.names
head(out, 10)

# Dereplication and learn the error rates (Default: nbases = 1e8)
# derepFastq() has been intergrated into learnErrors()
errF = learnErrors(filtFs, multithread = TRUE)
errR = learnErrors(filtRs, multithread = TRUE)

pdf("plotErrors.pdf", width = 10, height = 10, pointsize = 12)
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
dev.off()

# Sample Inference
# By default, the `dada` function processes each sample independently (pool = FALSE)
# Use `pool = TRUE` or `pool = pseudo` (recommended) if samples are from an extremely diverse community (e.g. soil)
dadaFs = dada(filtFs, err = errF, pool = FALSE, multithread = TRUE)
dadaRs = dada(filtRs, err = errR, pool = FALSE, multithread = TRUE)

# Merge paired reads (Default: minOverlap = 12; maxMismatch = 0)
mergers = mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

# Construct sequence table
seqtab = makeSequenceTable(mergers)

table(nchar(getSequences(seqtab)))

# Save sequence table
saveRDS(seqtab, "seqtab.rds") # or as an example, use seqtab[c(1:5),] to save data for a subset of the first 5 samples

# Save current workspace
# save.image(file="image.RData")
