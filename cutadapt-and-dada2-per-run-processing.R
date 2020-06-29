#

library("dada2")
library("ShortRead")
library("Biostrings")
library("ggplot2")
options(width = 190)

# Location of raw fastq files
# Note to change the R1 & R2 fastq file pattern accordingly
path.raw = file.path("raw")
fnFs = sort(list.files(path.raw, pattern = ".1.fastq.gz", full.names = TRUE))
fnRs = sort(list.files(path.raw, pattern = ".2.fastq.gz", full.names = TRUE))

# Set directory to place N-filterd files
path.filtN = file.path("filtN")
fnFs.filtN = file.path(path.filtN, basename(fnFs))
fnRs.filtN = file.path(path.filtN, basename(fnRs))

# Pre-filter to remove those with Ns, but perform no other filtering
out1 = filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)
out1 = as.data.frame(out1)
out1$Perc = out1$reads.out/out1$reads.in * 100
#head(out1)

# Set directory to place cutadapt-filterd files
path.cut = file.path("cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)

fnFs.cut = file.path(path.cut, basename(fnFs))
fnRs.cut = file.path(path.cut, basename(fnRs))
log.cut = gsub(".1.fastq.gz", ".out", fnFs.cut)
sample.names = gsub(".1.fastq.gz", "", basename(fnFs.cut))

# Define the primer set used to perform PCR
FWD = "CCTACGGRAGGCAGCAG"	# 341F
REV = "GACTACHVGGGTATCTAATCC"	# 805R

# Get reverse complement DNA sequences
FWD.RC = dada2::rc(FWD)
REV.RC = dada2::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags = paste("-g", FWD, "-a", REV.RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags = paste("-G", REV, "-A", FWD.RC)

# Run cutadapt to remove primers
# Note to change the PATH to cutadapt accordingly
cutadapt = "/path-to-bin/cutadapt"

for(i in seq_along(fnFs)) {
	print(sample.names[i])

	system2(cutadapt,
		stdout = log.cut[i], stderr = log.cut[i],	# log file
		args = c(R1.flags, R2.flags, 
		"-n 2",					# -n 2 required to remove FWD and REV from reads
		"-m 150",				# discard reads shorter than LEN (avoid length zero sequences)
#		"-j 0",					# auto-detection of CPU cores, only available on Python 3
		"-o", fnFs.cut[i], "-p", fnRs.cut[i],	# trimmed files
		fnFs.filtN[i], fnRs.filtN[i])		# input files
	)
}

# Plot quality profile of fastq files
ii = 1:length(sample.names)
pdf("plotQualityProfile.pdf", width = 8, height = 8, pointsize = 12)
for(i in ii) {
	print(sample.names[i])
	print(plotQualityProfile(fnFs.cut[i]) + ggtitle("Fwd"))
	print(plotQualityProfile(fnRs.cut[i]) + ggtitle("Rev"))
}
dev.off()

# Set directory to place dada2-filterd files
path.filt = file.path("filtered")
filtFs = file.path(path.filt, basename(fnFs.cut))
filtRs = file.path(path.filt, basename(fnRs.cut))

# Perform filtering and trimming
# Review "plotQualityProfile.pdf" to select the best paramters for truncLen
out2 = filterAndTrim(fnFs.cut, filtFs, fnRs.cut, filtRs,
	# Need to keep paramters consistent between runs of the same study
	truncLen = c(225,225), minLen = 150, maxN = 0, truncQ = 2, maxEE = c(2,5),
	rm.phix = TRUE, compress = TRUE, verbose = TRUE, multithread = TRUE)
out2 = as.data.frame(out2)
out2$Perc = out2$reads.out/out2$reads.in * 100
#head(out2)

# Dereplication and learn the error rates (Default: nbases = 1e8)
# derepFastq() has been intergrated into learnErrors()
errF = learnErrors(filtFs, multithread = TRUE)
errR = learnErrors(filtRs, multithread = TRUE)

pdf("plotErrors.pdf", width = 10, height = 10, pointsize = 12)
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
dev.off()

# Sample Inference
dadaFs = dada(filtFs, err = errF, multithread = TRUE)
dadaRs = dada(filtRs, err = errR, multithread = TRUE)

# Merge paired reads (Default: minOverlap = 12; maxMismatch = 0)
mergers = mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

# Construct sequence table
seqtab = makeSequenceTable(mergers)

# Save sequence table
saveRDS(seqtab, "seqtab.rds") # or as an example, use seqtab[c(1:5),] to save data for a subset of the first 5 samples

# Save current workspace
# save.image(file="image.RData")
