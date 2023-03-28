# copyright Pedro Martinez Arbizu & Sahar Khodami (2023)
# pmartinez@senckenberg.de
#  
# script for batch processing
# Feb 2020
# R script to process data using dada2
# dada2 commands based on https://benjjneb.github.io/dada2/tutorial.html

# USE like this:
# nohup Rscript --vanilla SGN_dada2_batch.r &
 

##########################################
# Edit this variables
##########################################
#working directory
setwd('/Dataexchange/dada2batch/')
#path to fastq files
path <- "/Dataexchange/dada2batch/rawseq"

patternF <- "_R1_tr.fastq"
patternR <- "_R2_tr.fastq"

#max number expected errors
maxEF <- 1
maxER <- 2

#truncate length
#> (60+306)/2
#[1] 183  for COI
#> (60+364)/2
#[1] 212 for V1V2

#tlenF <- 212
#tlenR <- 212

#mRNA
tlenF <- 150
tlenR <- 150

###########################################
# Do not edit beyond this line
###########################################

#load packages
library(dada2)

# here files names contain '_tr.' change to '_001.' if needed
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern=patternF, full.names = TRUE))
fnRs <- sort(list.files(path, pattern=patternR, full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#save the workspace
save.image('dada2.RData')

#if needed
# check quality of reads
pdf(file='quality_profileF.pdf')
plotQualityProfile(fnFs[1:2])
dev.off()

pdf(file='quality_profileR.pdf')
plotQualityProfile(fnRs[1:2])
dev.off()


# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(tlenF,tlenR), maxN=0, maxEE=c(maxEF,maxER), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) 

#learning errors on forward
errF <- learnErrors(filtFs, multithread=TRUE)
#learning errors on reverse
errR <- learnErrors(filtRs, multithread=TRUE)

#plot the trained error model
pdf(file='errorModelF.pdf')
plotErrors(errF, nominalQ=TRUE)
dev.off()
#plot the trained error model
pdf(file='errorModelR.pdf')
plotErrors(errR, nominalQ=TRUE)
dev.off()

# dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)


#merge paired reads 
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#build seqeunce table
seqtab <- makeSequenceTable(mergers)

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

#export the track data
write.csv(track,'trackdada2.cvs')

#export the dada2 sequences
write.table(colnames(seqtab.nochim),'nonchim.dada2.txt')

#export the community table
write.table(t(seqtab.nochim),'full.nonchim.dada2.txt',col.names=FALSE)

#save the workspace
save.image('dada2.RData')
