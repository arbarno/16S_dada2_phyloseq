z## This script provides the workflow to raw demultiplexed Illumina Reads using DADA2 with cutadapt primer removal
## This requires downloading the required packages before running, as well as downloading cutadapt
## DADA2 tutorials can be found here: https://benjjneb.github.io/dada2/
## Cutadapt documentation: https://cutadapt.readthedocs.io/en/stable/
### cutadapt installation defers significantly if using a Mac vs Windows machine.. read the documentation
### for using cutadapt, the primer sequence must be known
### if primer sequence is not known, remove primers using BBDuk (not described here)
## BBDuk documentation: https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/

## Also, special shout out to Nauras Daraghmeh (https://github.com/naurasd) who created the skeleton much of this is adapted from

# Install necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dada2") # if this does not work, try to install via devtools (requires prior installation of devtools)
BiocManager::install("ShortRead")
BiocManager::install("Biostrings")

# If dada2 installation failed via BiocManager, try to install via devtools (requires prior installation of devtools)
install.packages("devtools") # If not installed prior
library("devtools")
devtools::install_github("benjjneb/dada2")

# We will also use some functions of ggplot2 to plot read quality profiles.
# If not done prior, install ggplot2 via:
install.packages("ggplot2")
# Or if this fails:
install.packages("devtools") # If not installed prior
library("devtools")
devtools::install_github("tidyverse/ggplot2")

## DADA2 16S workflow with cutadapt primer removal ##

# Load required packages
library(dada2) # v1.36.0
library(ShortRead) # v1.67.0
library(Biostrings) # v2.77.0
library(ggplot2) # v3.5.2

# Read in the input files containing the sequence reads
# dada2 can deal with fasta / fastq / fastq.gz files as input
# Set the directory path containing files 
# Make sure path contains / instead of \

path <- "/path/to/fastq/folder"
list.files(path)

# Generate matched lists of the forward and reverse read files
#sometimes written as R1 and R2 fastq: SAMPLENAME_R1_001.fastq, SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fq.gz", full.names = TRUE))

# Check number of reads in the files
#countFastq(fnFs)
#countFastq(fnRs)

## Remove primers with cutadapt ##

# Designate primer sequences [including ambiguous nucleotides (base = N, Y, W, etc.) if present) 
# I included several common primers we use in the MMs lab
#FWD <- "ACTCCTACGGGAGGCAGCA"  # 338F 
#FWD <- "CCTACGGGNGGCWGCAG"  # 341F
FWD <- "CCTAYGGGRBGCASCAG"  # 341F BMK
#REV <- "GACTACHVGGGTATCTAATCC"  ## 785R
#REV <- "GGACTACHVGGGTWTCTAAT"  ## 806R
REV <- "GGACTACNNGGGTATCTAAT"  ## 806R BMK


# Verify the presence and orientation of these primers in the data
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

# Calculate number of reads containing forward and reverse primer sequences (Only exact matches are found.).
# Only one set of paired end fastq files will be checked (first sample in this case).
# This is is sufficient, assuming all the files were created using the same library preparation.
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), # [...] defines number of sample to be checked
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]])) # [...] defines number of sample to be checked

# Output interpretation:
# FWD primer is found in the forward reads in its forward orientation
# REV primer is found in the reverse reads in its forward orientation

# Define the path to cutadapt.exe file:
cutadapt <- "/path/to/cutadapt"

# To see if this worked and cutadapt has indeed been found, check installed version of cutadapt:
system2(cutadapt, args = "--version") # Run shell commands from R. See if R recognizes cutadapt and shows its version

# Create output filenames for the cutadapt-ed files.
# Define the parameters for the cutadapt command.
# See here for a detailed explanation of parameter settings: https://cutadapt.readthedocs.io/en/stable/guide.html#
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

# Trim FWD off of R1 (forward reads) - 
R1.flags <- paste0("-g", " ", FWD)
# I have removed the anchored primer identification from past experiences where anchoring has not worked. To anchor the primers, replace " " with " ^"

# Trim REV off of R2 (reverse reads)
R2.flags <- paste0("-G", " ", REV)
# I have removed the anchored primer identification from past experiences where anchoring has not worked. To anchor the primers, replace " " with " ^"

# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c("-e 0.07",R1.flags, R2.flags, # "-e 0.07" for an error rate of 7%.
                             "-m 1", # -m 1 discards reads having a length of zero bp after cutadapting to avoid problems in quality profile plotting later on
                             "--discard-untrimmed",
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs[i], fnRs[i])) # input files
}

# see here for a detailed explanation of the output:
# https://cutadapt.readthedocs.io/en/stable/guide.html#cutadapt-s-output
# Sometimes, you will see this: "WARNING: One or more of your adapter sequences may be incomplete. Please see the detailed output above."
# This usually refers to: "WARNING: The adapter is preceded by "T" (or any other base) extremely often. The provided adapter sequence could be incomplete at COI 3' end."
# The amplified regions and primer binding sites are usually highly conserved, so primer sequences are often preceded by the same base.
# Cutadapt just warns us that this is the case and tells us to check if the preceding base is indeed not part of the primer. 

# Count the presence of primers in the first cutadapt-ed sample to check if cutadapt worked as intended:
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

## The primer-devoid sequence files are now ready to be analyzed. ##

##-----------------------------------------------------------------------------------------------------------

# Similar to the earlier steps of reading in FASTQ files, read in the names of the cutadapt-ed FASTQ files. 
# Get the matched lists of forward and reverse fastq files.

# Forward and reverse fastq filenames have the format (make sure file name patterns are correct)
cutFs <- sort(list.files(path.cut, pattern = "_1.fq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_2.fq.gz", full.names = TRUE))

# Check number of reads in the files
#countFastq(cutFs)
#countFastq(cutRs)

# Check if forward and reverse files match:
if(length(cutFs) == length(cutRs)) print("Forward and reverse files match. Go forth and explore")
if (length(cutFs) != length(cutRs)) stop("Forward and reverse files do not match. Better go back and have a check")

# Generate quality profile plots for our reads.
# In case we have more than 20 fastq files, the following command will randomly choose 20 files to be plotted, otherwise all files will be plotted.
if(length(cutFs) <= 20) {
  fwd_qual_plots<-plotQualityProfile(cutFs) + 
    scale_x_continuous(breaks=seq(0,300,10)) + 
    scale_y_continuous(breaks=seq(0,40,5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    geom_hline(yintercept = 30)
  rev_qual_plots<-plotQualityProfile(cutRs) + 
    scale_x_continuous(breaks=seq(0,300,10)) + 
    scale_y_continuous(breaks=seq(0,40,5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    geom_hline(yintercept = 30)
} else {
  rand_samples <- sample(size = 20, 1:length(cutFs)) # grab 20 random samples to plot
  fwd_qual_plots <- plotQualityProfile(cutFs[rand_samples]) + 
    scale_x_continuous(breaks=seq(0,300,10)) + 
    scale_y_continuous(breaks=seq(0,40,5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    geom_hline(yintercept = 30)
  rev_qual_plots <- plotQualityProfile(cutRs[rand_samples]) + 
    scale_x_continuous(breaks=seq(0,300,10)) + 
    scale_y_continuous(breaks=seq(0,40,5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    geom_hline(yintercept = 30)
}
fwd_qual_plots
rev_qual_plots


# In case you want to plot reads of a particular sample, run:
#fwd_qual_plots1<-plotQualityProfile(cutFs[1]) + # Plots quality profile of first sample. You may also use cutFs[1:5], cutFs[c(2,4,7)] etc
#  scale_x_continuous(breaks=seq(0,300,10)) + 
#  scale_y_continuous(breaks=seq(0,40,5)) +
#  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
#  geom_hline(yintercept = 30)
#rev_qual_plots1<-plotQualityProfile(cutRs[1]) + # Plots quality profile of first sample. You may also use cutRs[1:5], cutRs[c(2,4,7)] etc
#  scale_x_continuous(breaks=seq(0,300,10)) + 
#  scale_y_continuous(breaks=seq(0,40,5)) +
#  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
#  geom_hline(yintercept = 30)
#fwd_qual_plots1
#rev_qual_plots1

# Print out the forward quality plot
jpeg(file="16S_run1_quality_forward.jpg",res=300, width=15, height=8, units="in")
fwd_qual_plots
dev.off()

# Print out the reverse quality plot
jpeg(file="16S_run1_quality_reverse.jpg",res=300, width=15, height=8, units="in")
rev_qual_plots
dev.off()

## Filter and trim ##

# Assigning the directory for the filtered reads to be stored in.
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

# Set filtering parameters #

# maxN = 0: dada2 does not allow reads with any ambiguous bases for downstream analysis, so has to be set to zero.
# truncLen: Truncate reads after n bases. Reads shorter than this are discarded.
# This truncates the low-quality ends off of our reads. Given that reverse reads show lower overall quality, they probably need to be truncated further.
# When specifying truncLen, remember that forward and reverse reads will be merged later and sufficient overlap needs to remain.
# After truncation, reads with higher than maxEE "expected errors" will be discarded.
# maxEE: After truncation, reads with higher than maxEE "expected errors" will be discarded. Expected errors are calculated from the nominal definition of the quality score: EE = sum(10^(-Q/10)).
# You may set a higher maxEE for reverse reads due to their lower quality.
# truncQ: Truncate reads at the first instance of a quality score less than or equal to truncQ.
# minLen: Remove reads with length less than minLen. It is enforced after trimming and truncation and filters spurious short reads.
# rm.phix=TRUE: Removes any reads that match the PhiX bacteriophage genome, which is typically added to Illumina sequencing runs for quality monitoring.
# Multithreading allows for parallel computation on several cores.
# On windows, multithreading is not possible for the filterAndTrim function. Setting multithread = TRUE is not a problem though, a simple warning will occur. 
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, truncLen=c(220,220), maxEE = c(2,2), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE) 

# check how many reads remain after filtering
out

# Save this output as RDS file:
saveRDS(out, "filter_and_trim_out_16S_run1.rds")

# Check number of reads in the files
#countFastq(filtFs)
#countFastq(filtRs)

# Dereplicate reads (this is a new addition after emailing Antony - https://github.com/cpantony)
# does not appear absolutely necessary with newer dada2 versions, but will include here
derepFs <- derepFastq(filtFs, verbose=T)
derepRs <- derepFastq(filtRs, verbose=T)

## Error model generation ##

# dada2 learns the specific error-signature of our dataset.
# This is why files from separate sequencing runs have to be processed separately till later on.
set.seed(1) # set seed to ensure that randomized steps are replicable
errF <- learnErrors(filtFs, multithread=T)
errR <- learnErrors(filtRs, multithread=T)

# save error calculation as RDS files:
saveRDS(errF, "errF_16S_run1.rds")
saveRDS(errR, "errR_16S_run1.rds")

# Visualize the estimated error rates:
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

# Apply the dada2's core sequence-variant inference algorithm:
dadaFs <- dada(derepFs, err=errF, multithread=T,pool="pseudo")
dadaRs <- dada(derepRs, err=errR, multithread=T,pool="pseudo")

# Save sequence-variant inference as RDS files which may be uploaded in case R crashes: 
saveRDS(dadaFs, "dadaFs_16S_run1.rds")
saveRDS(dadaRs, "dadaRs_16S_run1.rds")

# Inspecting the returned dada-class object of the first sample:
dadaFs[[1]]

# Merge the forward and reverse reads together to obtain the full denoised sequences.
# Adjust the minimum overlap (default = 12) and maximum mismatch allowed (e.g. = 1) if necessary.
# The minimum overlap is adjusted to 10
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap = 10, maxMismatch = 1, verbose=T)

# Inspect the merger data.frame of the first sample
# See here for output explanation: https://www.rdocumentation.org/packages/dada2/versions/1.0.3/topics/mergePairs
head(mergers[[1]])

saveRDS(mergers,"mergers_16S_run1.rds")

# Construct an amplicon sequence variant table (ASV) table
# If maxMismatch > 0 has been allowed in the mergePairs step,
# "Duplicate sequences detected and merged" may appear as output during the sequence table creation. This is not a problem.
seqtab <- makeSequenceTable(mergers)

# How many sequence variants were inferred?
dim(seqtab)

# Save sequence table
saveRDS(seqtab, "16S_seqtab_run1.rds")

## If you have read datasets of various sequencing runs, continue here:
## If samples were  obtained using a single sequencing run, use the dada2_raw_pipeline.R pipeline

# Create a table to track read numbers throughout the pipeline:

# There will be a problem creating this table if we have a sample with zero reads after the filter and trim step in our out object and omitted this sample during the downstream analysis.
# We create a new out object by subsetting all rows that contain no zeros.
# Later on, we may add the omitted samples manually to our tracking table in excel.
row_sub = apply(out, 1, function(row) all(row !=0 ))
out<-out[row_sub,]

# If no samples were omitted during the steps above, ignore the two previous lines and continue as follows:
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,getN))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged_run")
rownames(track) <- sample.names
track
write.table(track,"track_16S_run1.txt",sep="\t",col.names = NA)

##-----------------------------------------------------------------------------------------------------------

# Load the sequence table of the different sequence runs
run1<-readRDS("16S_seqtab_run1.rds")
run2<-readRDS("16S_seqtab_run2.rds")

# Inspect distribution of read lengths
hist(nchar(getSequences(run1)), main="Distribution of sequence lengths")
hist(nchar(getSequences(run2)), main="Distribution of sequence lengths")

# Merge the sequence tables of the separate sequencing runs
merged<-mergeSequenceTables(run1,run2)

# Save the merged sequence table
saveRDS(merged,"merged_seqtab.rds")

# Show overall length distribution of reads
table(nchar(getSequences(merged)))
hist(nchar(getSequences(merged)), main="Distribution of sequence lengths")

# The 341F - 805R 16S primer pair will generate an amplicon set with a bimodal length distribution.
# We will only keep sequence reads with a length of 400-430 bp.
seqtab.filtered <- merged[,nchar(colnames(merged)) %in% seq(400,430)]

## Remove chimeras ##

# Chimeras are PCR artefacts and may inflate diversity estimates.
# Chimeric sequences are identified if they can be exactly reconstructed by combining a left-segment and a right-segment from two more abundant "parent" sequences.
# If pooling was used during the sequence variant inference step above, method = "pool" should be specified here.
# In this case, the default setting minFoldParentOverabundance = 2 may be too stringent for pooled chimera removal, consider setting this to 8 (or maybe 4, 6, etc.).
seqtab.nochim <- removeBimeraDenovo(seqtab.filtered, multithread=T, verbose=TRUE)

# Save table with the non-chimeric sequences as rds-file:
saveRDS(seqtab.nochim, "16S_seqtab_nochim.rds")

# It is possible that a large fraction of the total number of UNIQUE SEQUENCES will be chimeras.
# However, this is usually not the case for the majority of the READS.
# Calculate percentage of the reads that were non-chimeric:
sum(seqtab.nochim)/sum(seqtab.filtered)

## If you want to taxonomically classify your ASVs without clustering, go ahead as follows:

# For 16S, pre-arranged reference databases in dada2 format exist.
# Check dada2's dedicated website: https://benjjneb.github.io/dada2/training.html
# We will use Silva v138: https://zenodo.org/records/14169026 (silva_nr99_v138.2_train_set.fa.gz)

# dada2 uses the RDP classifier applying a bootstrap method.
# Assign taxonomy with minimum Bootstrap confidence of 70
# tryRC can be set to TRUE to assign taxonomy using the reverse-complement of each sequence if it is a better match to the reference sequences than the forward sequence. This will increase computation time, though.
# outputBootstraps was set to TRUE to output bootstrap values in an integer matrix. A named list containing the assigned taxonomies (named "taxa") and the bootstrap values (named "boot") will be returned. Minimum bootstrap confidence filtering still takes place.
# Set seed for replicability (already happened for error rate estimation, repreat if pipeline has been terminated in between).
# boostrapping has not been necessary, but can still be included
taxa.70 <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.2_train_set.fa.gz", minBoot=70, multithread=T,outputBootstraps = T)

# Save the taxa70 object as rds file
saveRDS(taxa.70$tax,"taxa70.rds")

# Check the structure of the taxa.70 object
# The taxa.70$tax list contains the taxonomy
# The taxa.70$boot list contains the bootstrap values for each taxonomic level
str(taxa.70)

# Inspect the taxonomic assignments:
taxa.print <- taxa.70$tax  
# Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# Write a taxa table with short "ASV_..."-names instead of complete sequences.
rownames(taxa.70$tax) <- paste0("ASV", seq(nrow(taxa.70$tax)))

write.table(taxa.70$tax,file="tax_table_16S.txt",sep="\t", quote=F,col.names=NA)

# Write a bootstrap value table with short "ASV_..."-names instead of complete sequences.
rownames(taxa.70$boot) <- paste0("ASV", seq(nrow(taxa.70$boot)))

write.table(taxa.70$tax,file="boot_table_16S.txt",sep="\t", quote=F,col.names=NA)

#-------------------------------------------------------------------------------------------------------

### Continue here after taxonomy assignment ###

# Write a fasta file of the final, non-chimeric sequences with short >ASV... type headers:
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="")
}

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "16S_nochim_ASVs.fa")

# Write an ASV count table of the final, non-chimeric sequences with short >ASV... type names
seqtab.nochim.counts <- seqtab.nochim
colnames(seqtab.nochim.counts) <- paste0("ASV", seq(ncol(seqtab.nochim)))

#ASV_counts<-t(seqtab.nochim) # transposing table
write.table(seqtab.nochim.counts,file="16S_ASV_counts.txt",sep="\t", quote=F,col.names=NA)

# If you would like to have an ASV set excluding singletons as well, continue as follows:

# Remove singletons from the non-chimeric ASVs:
#Transform counts to numeric (as they will most likely be integers)
mode(seqtab.nochim) = "numeric"

# Subset columns with counts of > 1
seqtab.nochim.nosingle<-seqtab.nochim[,colSums(seqtab.nochim) > 1]

# Write a fasta file of the final, non-chimeric , non-singleton sequences with short >ASV... type headers
asv_seqs <- colnames(seqtab.nochim.nosingle)
asv_headers <- vector(dim(seqtab.nochim.nosingle)[2], mode="character")
for (i in 1:dim(seqtab.nochim.nosingle)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="")
}

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "16S_nochim_nosingle_ASVs.fa")

# Write an ASV count table of the final, non-chimeric, non-singleton sequences with short >ASV... type names
colnames(seqtab.nochim.nosingle) <- paste0("ASV", seq(ncol(seqtab.nochim.nosingle)))

#ASV_counts<-t(seqtab.nochim.nosingle) # transposing table
write.table(seqtab.nochim.nosingle,file="ASV_counts_nosingle.txt",sep="\t", quote=F,col.names=NA)

# You can then remove singleton ASVs from your taxonomy table outside of or within in R 
rownames(taxa.70$tax) <- paste0("ASV", seq_len(nrow(taxa.70$tax)))

# Subset the taxonomy table to only the ASVs in seqtab.nochim.nosingle
taxa.nosingle <- taxa.70$tax[colnames(seqtab.nochim.nosingle), , drop=FALSE]

write.table(taxa.nosingle, file = "tax_table_16S_nosingle.txt", sep = "\t", quote = F, col.names = NA)

# Make a table to track reads through the chimera as well as singleton removal step.
# Manually add this table outside of R to the track tables of the individual sequencing datasets produced with the previous script.
getN <- function(x) sum(getUniques(x))
track_nochim<-cbind(rowSums(merged),rowSums(seqtab.nochim),rowSums(seqtab.nochim.nosingle))
colnames(track_nochim) <- c("merged", "nonchim","nosingle")
track_nochim
write.table(track_nochim,"track_16S_nochim_nosingle.txt",sep="\t",col.names = NA)
