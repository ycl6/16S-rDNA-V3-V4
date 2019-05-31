#

library("dada2")
library("ggplot2")
options(width=190)

fastq = "trimmed"
filt = "filt"

fns <- sort(list.files(fastq, full.names = TRUE))
fnFs <- fns[grep(".1.fastq.gz", fns)]
fnRs <- fns[grep(".2.fastq.gz", fns)]

# Plot quality profile of fastq files
ii = 1:length(fnFs)
pdf("plotQualityProfile.pdf", width=10, height=10, pointsize=12)
for(i in ii) {
        print(plotQualityProfile(fnFs[i]) + ggtitle("Fwd"))
        print(plotQualityProfile(fnRs[i]) + ggtitle("Rev"))
}
dev.off()

if(!file_test("-d", filt)) dir.create(filt)

filtFs <- file.path(filt, basename(fnFs))
filtRs <- file.path(filt, basename(fnRs))

# Filtering and trimming
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(260,200), minLen=200, maxN=0, truncQ=2, maxEE=c(2,5), rm.phix=TRUE, compress=TRUE, verbose=TRUE, multithread=TRUE)
#head(out)

# Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- gsub(".1.fastq.gz","",names(derepFs)) # Change sample name
names(derepRs) <- gsub(".2.fastq.gz","",names(derepRs)) # Change sample name

# Learn the Error Rates
dadaFs.lrn <- dada(derepFs, err=NULL, selfConsist = TRUE, multithread=TRUE)
dadaRs.lrn <- dada(derepRs, err=NULL, selfConsist = TRUE, multithread=TRUE)

pdf("plotErrors.pdf", width=10, height=10, pointsize=12)
plotErrors(dadaFs.lrn[[1]], nominalQ=TRUE)
plotErrors(dadaRs.lrn[[1]], nominalQ=TRUE)
dev.off()

# Pooled sample Inference
errF <- dadaFs.lrn[[1]]$err_out
dadaFs <- dada(derepFs, err=errF, pool=TRUE, multithread=TRUE)
errR <- dadaRs.lrn[[1]]$err_out
dadaRs <- dada(derepRs, err=errR, pool=TRUE, multithread=TRUE)

# Merge paired reads (Default: minOverlap = 20; maxMismatch = 0)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Construct a sample-by-sequence observation matrix
seqtab <- makeSequenceTable(mergers)

# Save sequence table
saveRDS(seqtab, "seqtab.rds") # or as an example, use seqtab[c(1:5),] to save data for a subset of 5 samples

# Save current workspace
# save.image(file="image.RData")
