# Process raw sequences into phyloseq object for analyses

# Load dada2 and prep ####
library(dada2); packageVersion("dada2")
library(vegan)
library(ggplot2)
library(dplyr)

# File parsing - For this, we will use only the forward illumina reads - make sure to move fwd reads into their own directory for simplest processing
path <- "." # CHANGE to the directory containing your demultiplexed fastq files
filtpath <- file.path(path, "filtered") # Filtered files go into the filtered/ subdirectory
if(!file_test("-d", filtpath)) dir.create(filtpath) # make directory for filtered fqs if not already present
fns <- list.files(path)
fastqs <- fns[grepl(".fastq.gz$", fns)] # CHANGE if different file extensions or to target only certain sequences
rm(fns)
list.files()
fnFs <- sort(list.files(path, pattern="_R1_001", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "-"), `[`, 1)

# visualize a couple of fwd read quality profiles to help you decide reasonable filtration parameters
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Filter and trim ####
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

# learn error rates ####
errF <- learnErrors(filtFs, multithread=TRUE, MAX_CONSIST = 20, nbases = 1e10)
errR <- learnErrors(filtRs, multithread=TRUE, MAX_CONSIST = 20, nbases = 1e10)

# sanity check
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


# Dereplication ####
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


# Sample inference ####
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)


# Merge paired reads ####
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


# Make a sequence table ####
seqtab <- makeSequenceTable(mergers)


# Remove Chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)


# Track Reads through pipeline ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track = as.data.frame(track)
track$filter.loss = (track[,1]-track[,2])/track[,1]
write.csv(track, file = "../Output/read_counts_at_each_step.csv", row.names = TRUE)



# remove contaminants using controls ####
library(decontam)


# import metadata ####
meta = read.csv("../metadata.csv")[,1:5]
row.names(meta) <- meta$SampleName
row.names(seqtab.nochim)

# reorder metadata
meta = meta[order(row.names(meta)),]

# Find controlsamples (extraction negatives) ####
meta$controls <- meta$Island.Collected.From == "Blank"

# find contaminants
contams = isContaminant(seqtab.nochim, neg = meta$controls, normalize = TRUE)
table(contams$contaminant)
write.csv(contams, file = "../Output/likely_contaminants.csv", row.names = TRUE)


# remove them
seqtab.nochim = seqtab.nochim[,(which(contams$contaminant != TRUE))]
seqtab.nochim = seqtab.nochim[meta$controls == FALSE,]
meta = meta[meta$controls == FALSE,]

# Assign Taxonomy ####
# taxa <- assignTaxonomy(seqtab.nochim, "../tax/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "../tax/silva_species_assignment_v132.fa.gz")

write.csv(as.data.frame(seqtab.nochim), file = "../Output/SeqTable_no-chimera_no-contams.csv", row.names = TRUE, quote = FALSE)
saveRDS(seqtab.nochim, file = "../Output/clean_dada2_seqtable.RDS")
saveRDS(taxa, file = "../Output/Silva_Taxonomy_from_dada2.RDS")

seqtab.nochim <- readRDS(file = "../Output/clean_dada2_seqtable.RDS")
taxa <- readRDS(file = "../Output/Silva_Taxonomy_from_dada2.RDS")



# Build phylogeny ####
    # source("https://bioconductor.org/biocLite.R")
    # biocLite("DECIPHER")
library(DECIPHER)
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)


library(phangorn)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
saveRDS(treeNJ, file = "../Output/treeNJ.RDS")


## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
saveRDS(fitGTR, file = "../Output/FitGTR1.RDS")
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

fitGTR <- readRDS(file = "../Output/FitGTR1.RDS")

# Hand off to Phyloseq ####

library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(meta), 
               tax_table(taxa),
               phy_tree(fitGTR$tree))


names(meta)
plot_richness(ps, x="Island.Collected.From", measures = "Shannon")

saveRDS(ps, file = "../Output/clean_phyloseq_object.RDS")
ps = readRDS(file = "../Output/clean_phyloseq_object.RDS")


