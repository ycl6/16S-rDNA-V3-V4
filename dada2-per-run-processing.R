#

library("dada2")
library("ggplot2")
options(width = 190)

fastq = "trimmed"
filt = "filt"

fns = sort(list.files(fastq, full.names = TRUE))
fnFs = fns[grep(".1.fastq.gz", fns)]
fnRs = fns[grep(".2.fastq.gz", fns)]
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

if(!file_test("-d", filt)) dir.create(filt)

filtFs <- file.path(filt, basename(fnFs))
filtRs <- file.path(filt, basename(fnRs))

# Filtering and trimming
# Review "plotQualityProfile.pdf" to select the best paramters for truncLen
out = filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     # Need to keep paramters consistent between runs of the same study
                     truncLen = c(260,200), minLen = 200, maxN = 0, truncQ = 2, maxEE = c(2,5),
                     rm.phix = TRUE, compress = TRUE, verbose = TRUE, multithread = TRUE)

out = as.data.frame(out)
rownames(out) = sample.names
head(out, 10)

# Learn the Error Rates
# The derepFastq function used in past workflow has been intergrated into learnErrors function
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

# Merge paired reads (Default: minOverlap = 20; maxMismatch = 0)
mergers = mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

# Construct a sample-by-sequence observation matrix
seqtab = makeSequenceTable(mergers)

table(nchar(getSequences(seqtab)))

# Save sequence table
saveRDS(seqtab, "seqtab.rds") # or as an example, use seqtab[c(1:5),] to save data for a subset of 5 samples

# Save current workspace
# save.image(file="image.RData")
