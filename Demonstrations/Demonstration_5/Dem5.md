# Chapter 6: Community Diversity Measures and Calculations                      
Complete the following demonstration in RStudio. Create a markdown file of your script. You can read about diversity indices in Xia et al. (2018), Chapter 6: Community Diversity Measures and Calculations

# Estimate alpha-diversity

# Use different R packages to estimate several alpha-diversity and plot the results

library(BiocManager)
#source("https://bioconductor.org/install")
#useDevel()
BiocManager::install(version='devel')
BiocManager::install("microbiome")

library(microbiome)
library(phyloseq)
library(knitr)

data(GlobalPatterns)

# Rename GlobalPatterns data (which is a phyloseq object)
ps1 <- GlobalPatterns
ps1

# Remove taxa represented by <1000 reads to create a smaller dataset
physeq = prune_taxa(taxa_sums(ps1) > 1000, ps1)
physeq

# microbiome: The following R commands calculate Chao 1 richness
chao1 <- alpha(physeq, index = "chao1")
head(chao1,3)

# microbiome: The following R commands calculate Chao 1 richness and the observed taxa
rich <- richness(physeq)
head(rich,3)

# microbiome: The following R commands calculate pielou index
pielou <- alpha(physeq, index = "pielou")
head(pielou,3)

# microbiome: The following R commands calculate all available evenness measures 
even <- evenness(physeq, "all")
head(even,3)

# phyloseq: The following R commands estimate multiple indices
diver<-estimate_richness(physeq, measures=c("Observed", "Chao1"))
head(diver,3)

diver<-estimate_richness(physeq)
head(diver,3)

# Estimate phylogenetic diversity (PD) using picante R package (and observed taxa)
install.packages("picante")
library(picante)

# Extract a transposed OTU table
otuD<-as.data.frame(t(otu_table(physeq)))

# Estimates phylogenetic diversity including a rooted tree via midpoint rooting
PD<-pd(otuD, phy_tree(physeq), include.root=TRUE)
head(PD)

# Plot alpha-diversity estimates

# Generating Boxplot
library(ggpubr)

# Add diversity indices to metadata
diver_all<-cbind(sample_data(physeq),PD,diver)
head(diver_all,3)

# generate boxplots using the ggboxplot() function from the ggpubr package
p <- ggboxplot(diver_all, x = "SampleType", y = "Shannon", color = "SampleType", add = "jitter", shape = "SampleType")
p

# generate boxplots using ggplot
chao <- ggplot(diver_all, aes(factor(SampleType), Chao1)) +
        geom_boxplot(aes(fill = factor(SampleType)),outlier.colour = "black", outlier.size = 1)+ geom_jitter(size=1,shape=1)+ ggtitle("Chao1 richness")+labs(y = "Chao1 richness")
chao

# Summarize the Diversity Measures per group using FSA R package
library(FSA)
Summarize(Shannon ~ SampleType, data = diver_all)

# Plot Histogram of the Diversity Distributions
library(lattice)
histogram(~ Shannon|SampleType, data=diver_all,layout=c(3,3))

# Estimate beta-diversity

# Use phyloseq to estimate beta-diversity indices

# Bray-Curtis Dissimilarity
brayd<-phyloseq::distance(physeq, method="bray")

# Jaccard Dissimilarity
jaccd<-phyloseq::distance(physeq, method="jaccard")

# Phylogenetic Beta Diversity: Unweighted UniFrac
uniun<-phyloseq::distance(physeq, method="unifrac")

# # Phylogenetic Beta Diversity: Weighted UniFrac
uniweigh<-phyloseq::distance(physeq, method="wunifrac")

# Save distance matrix as table
library(reshape)
brayd_m <- as.matrix(brayd)
head(brayd_m)
write.csv(brayd_m, file = "brayd_m.csv")
