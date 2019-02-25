# Community analyses

# Load dada2 and prep ####
library(vegan)
library(ggplot2)
library(dplyr)
library(phyloseq)
library(gridExtra)
library(ggpubr)
library(broom)
library(stargazer)
library(sjPlot)
library(RColorBrewer)
source("plot_bar2.R")
source("summarize_taxa.R")

# Color palette
pal = c("#6b5456","#ec8d1b","#6abf2a","#8b53b7","#70acbe","#01c95b","#c00014","#31332f","#f7d000","#abba00")

# load processed data
ps = readRDS(file = "./Output/clean_phyloseq_object.RDS")

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
ggsave(filename = "./Output/abundance_vs_prevalence_Phylum.png", dpi = 300, width = 12, height = 12)

#  Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps2)

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps3 = prune_taxa(keepTaxa, ps2)

# Remove samples not associated with an island
ps3 = prune_samples(ps3@sam_data$Island.Collected.From != "", ps3)

# Normalize (relative abundance) ####
ps3ra = transform_sample_counts(ps3, function(x){x / sum(x)})

# Plot glommed tree
ps4 = tip_glom(ps3ra, h = .01)
plot_tree(ps4, color = "Island.Collected.From")  + scale_color_discrete(name = "Island")
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


# abundance before and after normalization
plotBefore = plot_abundance(ps3,"") + ggtitle("Before Normalization") + theme(axis.text.x = element_text(angle = 75, hjust = 1))
ggsave(plotBefore, filename = "./Output/Abundance_before_norm.png", dpi= 300, width = 12, height = 10)
plotAfter = plot_abundance(ps3ra,"") + ggtitle("After Relative") + theme(axis.text.x = element_text(angle = 75, hjust = 1))
ggsave(plotAfter, filename = "./Output/Abundance_after_norm.png", dpi=300, width = 12, height = 10)



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
    labs(x="Orientation",y="Relative Abundance") +
    theme(axis.text.x = element_text(angle = 75, hjust = 1, size = 8)) +
    scale_color_manual(values = pal)
  
ggsave(relabundplot, filename = "./Output/Phylum_relabund_by_Orientation.png",dpi=300)



# Ordination(s) ####
NMDS = ordinate(ps3ra, method = "NMDS", distance = "bray")
DPCoA = ordinate(ps3ra, method = "DPCoA", distance = "bray")
PCoA = ordinate(ps3ra, method = "PCoA", distance = "bray")
DCA = ordinate(ps3ra, method = "DCA", distance = "bray")
CCA = ordinate(ps3ra, method = "CCA", distance = "bray")
RDA = ordinate(ps3ra, method = "RDA", distance = "bray")

p_NMDS = plot_ordination(ps3ra, NMDS, color = "Island.Collected.From", shape = "Orientation") + scale_color_discrete(name = "Island") + 
  scale_color_discrete(name = "Site") + theme_bw() + scale_color_manual(values = pal)

p_DPCoA = plot_ordination(ps3ra, DPCoA, color = "Island.Collected.From", shape = "Orientation") + ggtitle("DPCoA") + scale_color_discrete(name = "Site")+theme_bw()+scale_color_manual(values = pal)
p_PCoA = plot_ordination(ps3ra, PCoA, color = "Island.Collected.From", shape = "Orientation") + ggtitle("PCoA")  + scale_color_discrete(name = "Site")+theme_bw()+scale_color_manual(values = pal)
p_DCA = plot_ordination(ps3ra, DCA, color = "Island.Collected.From", shape = "Orientation") + ggtitle("DCA")  + scale_color_discrete(name = "Site")+theme_bw()+scale_color_manual(values = pal)
p_CCA = plot_ordination(ps3ra, CCA, color = "Island.Collected.From", shape = "Orientation") + ggtitle("CCA")  + scale_color_discrete(name = "Site")+theme_bw()+scale_color_manual(values = pal)
p_RDA = plot_ordination(ps3ra, RDA, color = "Island.Collected.From", shape = "Orientation") + ggtitle("RDA")  + scale_color_discrete(name = "Site")+theme_bw()+scale_color_manual(values = pal)

ggsave(p_NMDS, filename = "./Output/NMDS_w_corrected_names.png", dpi = 300)


ord_plots = ggarrange(p_NMDS,p_DPCoA,p_PCoA,p_DCA,p_CCA,p_RDA,common.legend = TRUE, legend = "right")
ggsave(ord_plots, filename = "./Output/Ordination_Plots_by_Orientation.png", dpi=300,width = 10,height = 8,units = "in")


# Heatmap ####
plot_heatmap(ps3ra, method = "NMDS", distance="bray", 
             sample.label="Orientation", 
             taxa.label="Family", 
             sample.order = "Orientation",
             taxa.order = "Family",
             low="#FFFFCC", high="#000033", na.value="white")
ggsave("./Output/Heatmap_of_Family_by_Orientation.png", dpi=300,width = 10, height = 12)

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
# otu = otu[,-c(emptytax)]

# plot and save heatmap
png(filename = "./Output/Heatmap_Family_by_Island_test.png", width = 960, height = 960, pointsize = 14)
heatmap(t(as.matrix(otu)), ColSideColors = as.character(samp$island.colors), Rowv = NA,
        labRow = tax$Family, labCol = samp$Island.Collected.From, margins = c(10,5), cexRow = 1, cexCol = 1,
        col=gray.colors(100))
dev.off()

# PERMANOVA TEST ####
ps = readRDS("./Output/clean_phyloseq_object.RDS")
ps = prune_samples(ps@sam_data$Orientation != "", ps)

otu = as.data.frame(otu_table(ps))
meta = as.data.frame(sample_data(ps))
df = data.frame(SampleID = meta$SampleID, Island = meta$Island.Collected.From, Species = meta$Species, Orientation = meta$Orientation)

permanova = adonis(otu ~ Orientation * Island, data = df)


sink("./Output/adonis_table.txt")
noquote(print("PermANOVA Table:"))
permanova
sink(NULL)

source("pairwise_adonis.R")
padonis = pairwise.adonis(otu,as.character(meta$Island.Collected))
sink("./Output/adonis_table.txt",append = TRUE)
noquote(print("Pairwise adonis between islands (Bonferroni corrected Pvalues:"))
padonis
sink(NULL)

dd2 = readRDS(file = "./Output/clean_dada2_seqtable.RDS")


# NMDS without taxon4
newotu = select(otu,-4)
mnds1 = metaMDS(otu)
MDS1 = mnds1$points[,1]
MDS2 = mnds1$points[,2]
island = samp$Island.Collected.From


ggplot(mapping = aes(x=MDS1,y=MDS2,color=meta$Orientation)) +
  geom_point() + theme_bw() +
  scale_color_discrete(name = "Orientation")
ggsave("./Output/NMDS_colored-by-Orientation.png", dpi = 300, height = 6, width = 6)

ggplot(mapping = aes(x=MDS1,y=MDS2,color=meta$Island, shape = meta$Orientation)) +
  geom_point() + theme_bw() +
  scale_color_discrete(name = "Site") +
  scale_shape_discrete(name = "Orientation") + scale_color_manual(values = pal)
ggsave("./Output/NMDS_wo_Taxon-4_Site-and-Orientation.png", dpi = 300, height = 6, width = 6)


# Shannon Diversity plot
div = data.frame(Site = ps3ra@sam_data$Island.Collected.From, 
                 Shannon = diversity(otu_table(ps3ra))) 
                 # Richness = colSums(decostand(otu_table(ps3ra), method = "pa"))
                 
write.csv(div, file = "../Output/Diversity_table.csv", quote = FALSE)



div %>% group_by(Site) %>%
  summarise(N = n(), Mean = mean(Shannon))

ggplot(div, aes(x=Site, y=Shannon, fill = Site)) +
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  labs(y="Shannon Diversity") + scale_color_manual(values = pal)
ggsave(filename = "./Output/Shannon_Diversity.png", dpi = 300)  


# Plot diversity bar charts ####
bact_barplot = plot_bar2(ps3ra, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "identity") + labs(y="Relative Abundance",x = "Sample ID") + 
  scale_fill_manual(values = pal) + scale_color_manual(values = pal) +
  coord_flip()
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

phylum_island = rbind(Hantu,Jong,Kusu,`Raffles Lighthouse`,Semakau,Sisters,`St John`,`Sultan Shoal`,TPT)

for(island in levels(ps3@sam_data$Island.Collected.From)){
  df = (summarize_taxa(subset_samples(ps3, Island.Collected.From == island), Rank = "Order"))
  df$Island = island
  df$RA = df$meanRA / sum(df$meanRA)
  assign(island, df, envir = .GlobalEnv)   
  print(island)
}

order_island = rbind(Hantu,Jong,Kusu,`Raffles Lighthouse`,Semakau,Sisters,`St John`,`Sultan Shoal`,TPT)

ggplot(phylum_island, aes(y=((RA)), x = Island, fill = Phylum)) +
  geom_bar(stat="identity") + coord_flip() + theme_bw() + labs(x= "Island", y = "Relative Abundance") +
  scale_fill_manual(values = pal) 
ggsave(filename = "./Output/Relabund_by_Island.png", dpi = 300, height = 6.07, width = 10)

ggplot(order_island, aes(y=((RA)), x = Island, fill = Order)) +
  geom_bar(stat="identity") + coord_flip() + theme_bw() + labs(x= "Island", y = "Relative Abundance") +
ggsave(filename = "./Output/Relabund_order_by_Island.png", dpi = 300, height = 6.07, width = 16)


