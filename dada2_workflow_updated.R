

# library("devtools")
# devtools::install_github("benjjneb/dada2")
# devtools::install_github("joey711/phyloseq")

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

# source('http://bioconductor.org/biocLite.R')
# biocLite('phyloseq')
# library(devtools)
# devtools::install_github("benjjneb/decontam")

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

# quick peek at ps object
sample_names(ps)
rank_names(ps)
sample_variables(ps)
otu_table(ps)[1:5, 1:5]
tax_table(ps)[1:5, 1:7]


# Change sequences to unique IDs to make viewing easier
seqs_16S = taxa_names(ps)
names(seqs_16S) <- 1:length(seqs_16S) # save seqs and IDs combination
taxa_names(ps) <- names(seqs_16S)
tax_table(ps)[1:5, 1:7]




# quick exploratory look at top families by island
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
ps.islands = prune_samples(ps.top20@sam_data$Island.Collected.From != "",ps.top20)
plot_bar(ps.top20, x="Island.Collected.From", fill="Family")


table(tax_table(ps)[, "Phylum"], exclude = NULL)
table(tax_table(ps)[, "Family"], exclude = NULL)


# subset ps to remove mitochondrial and chloroplast sequences, missing phyla and those of low <5% prevalence ####
ps0 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps0),
               MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps0),
                    tax_table(ps0))

# compute mean and total prevalence
phylum_prev = plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
family_prev = plyr::ddply(prevdf, "Family", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# Remove "mitochondria" taxa
filterFamily = c("Mitochondria")
ps1 = subset_taxa(ps0, !Family %in% filterFamily)

# Remove "Chloroplast" taxa
filterOrder = c("Chloroplast")
ps2 = subset_taxa(ps1, !Order %in% filterOrder)

# Look at abundance vs prevalence
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps2, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps0),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.025, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
ggsave(filename = "../Output/abundance_vs_prevalence_Phylum.png", dpi = 300)

#  Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps2)
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps3 = prune_taxa(keepTaxa, ps2)
ps3

# Remove samples not associated with an island
ps3 = prune_samples(ps3@sam_data$Island.Collected.From != "", ps3)




# # Look at tree ####
# library(ape)
# p2tree = plot_tree(ps3, method = "treeonly",
#                    ladderize = "left",
#                    title = "Before Agglomeration")
# 
# plot_tree(ps0)

# ggsave(filename = "../Output/Glommed_Tree_by_Island.png", dpi=300, width = 10, height = 8)

# Normalize (relative abundance) ####
ps3ra = transform_sample_counts(ps3, function(x){x / sum(x)})

# Plot glommed tree
ps4 = tip_glom(ps3ra, h = .01)
plot_tree(ps4, color = "Island.Collected.From") + aes(color = "Island.Collected.From") + scale_color_discrete(name = "Island")
plot_tree(ps4, color = "Phylum") 


# Plot abundance function ####
plot_abundance = function(physeq,title = "", 
                          Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = physeq
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "Orientation",y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}

library(gridExtra)

plotBefore = plot_abundance(ps3,"") + ggtitle("Before Normalization") + theme(axis.text.x = element_text(angle = 75, hjust = 1))
ggsave(plotBefore, filename = "../Output/Abundance_before_norm.png", dpi= 300)
plotAfter = plot_abundance(ps3ra,"") + ggtitle("After Relative") + theme(axis.text.x = element_text(angle = 75, hjust = 1))
ggsave(plotAfter, filename = "../Output/Abundance_after_norm.png", dpi=300)
# grid.arrange(plotBefore,plotAfter)

# Make final figure of relative abundance of phyla by island ####
  # Arbitrary subset, based on Phylum, for plotting
  
  mphyseq = psmelt(ps3ra)
  mphyseq <- subset(mphyseq, Abundance > 0)
  
  relabundplot = ggplot(data = mphyseq, mapping = aes_string(x = "Orientation",y = "Abundance",
                                              color = "Phylum")) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = "Phylum") + scale_y_log10()+
    theme_bw() + theme(legend.position="none") + 
    labs(x="Island",y="Relative Abundance") +
    theme(axis.text.x = element_text(angle = 75, hjust = 1, size = 8)) 
ggsave(relabundplot, filename = "../Output/Phylum_relabund_by_Orientation.png",dpi=300)

  
# Ordination(s) ####
NMDS = ordinate(ps3ra, method = "NMDS", distance = "bray")
DPCoA = ordinate(ps3ra, method = "DPCoA", distance = "bray")
PCoA = ordinate(ps3ra, method = "PCoA", distance = "bray")
DCA = ordinate(ps3ra, method = "DCA", distance = "bray")
CCA = ordinate(ps3ra, method = "CCA", distance = "bray")
RDA = ordinate(ps3ra, method = "RDA", distance = "bray")

p_NMDS = plot_ordination(ps3ra, NMDS, color = "Island.Collected.From", shape = "Orientation") + scale_color_discrete(name = "Island") + 
  scale_color_discrete(name = "Site") + 
    theme_bw()
p_DPCoA = plot_ordination(ps3ra, DPCoA, color = "Island.Collected.From", shape = "Orientation") + ggtitle("DPCoA") + scale_color_discrete(name = "Site")+theme_bw()
p_PCoA = plot_ordination(ps3ra, PCoA, color = "Island.Collected.From", shape = "Orientation") + ggtitle("PCoA")  + scale_color_discrete(name = "Site")+theme_bw()
p_DCA = plot_ordination(ps3ra, DCA, color = "Island.Collected.From", shape = "Orientation") + ggtitle("DCA")  + scale_color_discrete(name = "Site")+theme_bw()
p_CCA = plot_ordination(ps3ra, CCA, color = "Island.Collected.From", shape = "Orientation") + ggtitle("CCA")  + scale_color_discrete(name = "Site")+theme_bw()
p_RDA = plot_ordination(ps3ra, RDA, color = "Island.Collected.From", shape = "Orientation") + ggtitle("RDA")  + scale_color_discrete(name = "Site")+theme_bw()

ggsave(p_NMDS, filename = "../Output/NMDS_w_corrected_names.png", dpi = 300)

library(ggpubr)

ord_plots = ggarrange(p_NMDS,p_DPCoA,p_PCoA,p_DCA,p_CCA,p_RDA,common.legend = TRUE, legend = "right")
ggsave(ord_plots, filename = "../Output/Ordination_Plots_by_Orientation.png", dpi=300,width = 10,height = 6,units = "in")


# Heatmap ####

?plot_heatmap

plot_heatmap(ps3ra, method = "NMDS", distance="bray", 
             sample.label="Orientation", 
             taxa.label="Family", 
             sample.order = "Orientation",
             taxa.order = "Family",
             low="#FFFFCC", high="#000033", na.value="white")

# Extract data frames from physeq object
otu = as.data.frame(otu_table(ps3ra))
tax = as.data.frame(tax_table(ps3ra))
samp = as.data.frame(sample_data(ps3ra))

# add color vector to samp
island.names = levels(samp$Island.Collected.From)
color.names = colors()[c(121,20,96,40,80,60,200,15,105)]
cols = plyr::mapvalues(samp$Island.Collected.From, from = island.names, to = color.names)
samp$island.colors = cols
rm(cols)

# remove empty taxa
emptytax = which(colSums(otu) == 0)
tax = tax[-emptytax,]
otu = otu[,-emptytax]

# plot and save heatmap
png(filename = "../Output/Heatmap_Family_by_Island", width = 960, height = 960, pointsize = 14)
heatmap(t(as.matrix(otu)), ColSideColors = as.character(samp$island.colors), Rowv = NA,
        labRow = tax$Family, labCol = samp$Island.Collected.From, margins = c(10,5), cexRow = 1, cexCol = 1,
        col=gray.colors(100))
dev.off()

# PERMANOVA TEST ####
ps = readRDS("../Output/clean_phyloseq_object.RDS")
ps = prune_samples(ps@sam_data$Orientation != "", ps)

otu = as.data.frame(otu_table(ps))
meta = as.data.frame(sample_data(ps))
df = data.frame(SampleID = meta$SampleID, Island = meta$Island.Collected.From, Species = meta$Species, Orientation = meta$Orientation)
library(vegan)
sink("../Output/adonis_table.txt")
adonis(otu ~ Orientation, data = df)
sink(NULL)

# Test for species driving Island patterns ####
# Remove islands with only one observation

table(samp$Island.Collected.From) > 1
samp$Island.Collected.From


samp1 = filter(samp, Island.Collected.From %in% c("Hantu","Jong","Kusu","Raffles Light House","Semakau",
                                          "Sisters","TPT") )
otu1 = otu[row.names(otu) %in% samp1$SampleName,]

# Simper function
simp = simper(otu1,samp1[,5])
summary(simp)
data("dune")
data("dune.env")

sink(file = "../Output/Simper_Summary.txt")
summary(simp)
sink(NULL)


simp$`Raffles Light House_Kusu`$average > .1
simp$TPT_Kusu$average > .1
simp$Semakau_Kusu$average > .1
simp$Jong_Kusu$average > .1
simp$Hantu_Kusu$average > .1


simp$TPT_Sisters$average
simp$Semakau_Sisters$average > .1
simp$Jong_Sisters$average > .1
simp$`Raffles Light House_Sisters`$average > .1
simp$Hantu_Sisters$average > .1

simp$Sisters_Kusu

tax[1,]

dd2 = readRDS(file = "../Output/clean_dada2_seqtable.RDS")
colnames(dd2)[4]
library(dplyr)
# NMDS without taxon4
newotu = select(otu,-4)
mnds1 = metaMDS(otu)
MDS1 = mnds1$points[,1]
MDS2 = mnds1$points[,2]
island = samp$Island.Collected.From
library(ggplot2)

ggplot(mapping = aes(x=MDS1,y=MDS2,color=meta$Orientation)) +
  geom_point() + theme_bw() +
  scale_color_discrete(name = "Orientation")
ggsave("../Output/NMDS_wo_Taxon-4_colored-by-Orientation.png", dpi = 300, height = 6, width = 6)

ggplot(mapping = aes(x=MDS1,y=MDS2,color=meta$Island, shape = meta$Orientation)) +
  geom_point() + theme_bw() +
  scale_color_discrete(name = "Site") +
  scale_shape_discrete(name = "Orientation")
ggsave("../Output/NMDS_wo_Taxon-4_Site-and-Orientation.png", dpi = 300, height = 6, width = 6)


# Shannon Diversity plot

div = data.frame(Site = ps3ra@sam_data$Island.Collected.From, 
                 Shannon = diversity(otu_table(ps3ra))) 
                 # Richness = colSums(decostand(otu_table(ps3ra), method = "pa"))
                 

length(Shannon)
length(Site)

?diversity



write.csv(div, file = "../Output/Diversity_table.csv", quote = FALSE)



div %>% group_by(Site) %>%
  summarise(N = n(), Mean = mean(Shannon))

ggplot(div, aes(x=Site, y=Shannon, fill = Site)) +
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  labs(y="Shannon Diversity")

ggsave(filename = "../Output/Shannon_Diversity.png", dpi = 300)  




# Plot diversity bar charts ####
bact_barplot = plot_bar(ps3ra, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "identity") + labs(y="Relative Abundance",x = "Sample ID") + 
  coord_flip() + 
ggsave(bact_barplot, filename = "../Output/RelAbund_Barplot_by_Sample.png", dpi = 300, height = 10, width = 8)  

bact_barplot_order = plot_bar(ps3, fill = "Order") + 
  geom_bar(aes(color=Order, fill=Order), stat = "identity") + labs(y="Relative Abundance",x = "Sample ID") + 
  coord_flip()
ggsave(bact_barplot_order, filename = "../Output/RelAbund_Barplot_Order_by_Sample.png", dpi = 300, height = 10, width = 20)  

# Same, but mean for each island
phylum_summary = summarize_taxa(ps3, Rank = "Phylum", GroupBy = "Island.Collected.From")
ggplot(phylum_summary, aes(y=meanRA, x = Island.Collected.From, fill = Phylum)) +
  geom_bar(stat="identity") + coord_flip() + theme_bw()



for(island in levels(ps3@sam_data$Island.Collected.From)){
  df = (summarize_taxa(subset_samples(ps3, Island.Collected.From == island), Rank = "Phylum"))
  df$Island = island
  df$RA = df$meanRA / sum(df$meanRA)
  assign(island, df, envir = .GlobalEnv)   
  print(island)
}
Sisters

phylum_island = rbind(Hantu,Jong,Kusu,`Raffles Light House`,Semakau,Sisters,`St John`,`Sultan Shoals`,TPT)

for(island in levels(ps3@sam_data$Island.Collected.From)){
  df = (summarize_taxa(subset_samples(ps3, Island.Collected.From == island), Rank = "Order"))
  df$Island = island
  df$RA = df$meanRA / sum(df$meanRA)
  assign(island, df, envir = .GlobalEnv)   
  print(island)
}

order_island = rbind(Hantu,Jong,Kusu,`Raffles Light House`,Semakau,Sisters,`St John`,`Sultan Shoals`,TPT)

ggplot(phylum_island, aes(y=((RA)), x = Island, fill = Phylum)) +
  geom_bar(stat="identity") + coord_flip() + theme_bw() + labs(x= "Island", y = "Relative Abundance")
ggsave(filename = "../Output/Relabund_by_Island.png", dpi = 300, height = 6.07, width = 10)

ggplot(order_island, aes(y=((RA)), x = Island, fill = Order)) +
  geom_bar(stat="identity") + coord_flip() + theme_bw() + labs(x= "Island", y = "Relative Abundance")
ggsave(filename = "../Output/Relabund_order_by_Island.png", dpi = 300, height = 6.07, width = 16)


