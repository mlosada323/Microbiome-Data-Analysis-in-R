# Demonstration 3
Complete the following demonstration in RStudio. Xia et al. (2018), Chapter 7: Exploratory Analysis of Microbiome Data and Beyond provides some directions, but more information can be found at the https://bioconductor.org/packages/release/bioc/html/phyloseq.html. Review them to interpret scripts and outcomes of your analyses

# Introduction to phyloseq and microbiome packages

## Phyloseq package
phyloseq provides a set of classes and tools to facilitate the import, storage, analysis, and graphical display of microbiome census data

```r
# Here, we use some datasets to illustrate how to create a phyloseq object and some operations on phyloseq object. First, we need to install this package using BiocManager

# install physloseq
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")
packageVersion("phyloseq")
library(phyloseq)

# Set the R working directory to the folder where you'll store the datasets and load them
setwd("your directory")

tax_tab <- read.delim("new_tax_tab.txt")
otu_tab <- read.delim("new_otu_tab.txt", row.names = 1)
meta_tab <- read.csv("new_meta_tab.csv", header=TRUE, row.names = 1)

# Check and make sure that the datasets have been correctly loaded
head(otu_tab,3) # are taxa rows or not? see below
head(tax_tab)
head(meta_tab,3)

# Check and make sure that the classes of datasets are all data.frame
class(otu_tab)
class(tax_tab)
class(meta_tab)

# Since the above codes show that both otu_table and tax_table components are data frames, we need to convert them into matrices
otumat<-as.matrix(otu_tab)
taxmat<-as.matrix(tax_tab)

class(otumat)
class(taxmat)

# Create phyloseq object using phyloseq () 
otu<-otu_table(otumat,taxa_are_rows = TRUE)
tax<-tax_table(taxmat)
sam<-sample_data(meta_tab)
physeq <- phyloseq(otu, tax, sam)
physeq

ps1 <- prune_taxa(taxa_sums(physeq) > 10000, physeq);ps1
new_metadata <- meta(ps1)
write.csv(new_metadata,file="new_meta_tab.csv")
table_otu<-otu_table(ps1)
write.csv(table_otu,file="new_otu_tab.csv")
table_tax<-tax_table(ps1)
write.csv(table_tax,file="new_tax_tab.csv")

# check names
sample_names(otu)
sample_names(sam)
sample_variables


# Add a random phylogenetic tree component
library("ape")
random_tree = ape::rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
plot(random_tree)

# Now merge the tree data to the phyloseq object we already have by using the merge_phyloseq()
ps = merge_phyloseq(physeq, random_tree)
ps

# create phyloseq object in one step
ps <- phyloseq(tax_table(as.matrix(tax_tab)), phy_tree(tree), sample_data(meta_tab), otu_table(otu_tab, taxa_are_rows = TRUE))
ps

# phyloseq object information
str(ps)
nsamples(ps)
sample_names(ps)
min(colSums(otu_table(ps)))
sample_sums(ps)
sample_variables(ps)
phy_tree(ps)
taxa_names(ps)[1:10]
head(sample_data(ps))
head(otu_table(ps))
head(tax_table(ps))
rank_names(ps)
get_taxa_unique(ps, "Genus")

# manipulate phyloseq object

# separate phylum Proteobacteria
ps_prot <- subset_taxa(ps, Phylum=="Proteobacteria")
ps_prot

# separate male samples (see gender)
ps_male <- subset_samples(ps, gender == "M")
ps_male

# agglomerate taxa by genus
ps_gen <- tax_glom(ps, taxrank = "Genus")
ps_gen

# prune taxa with <50000 reads
ps1 <- prune_taxa(taxa_sums(ps) > 50000, ps)
ps1

# prune samples with <10000 reads
ps2 = prune_samples(sample_sums(ps) >= 10000, ps)
ps2

# transform otu counts to otu proportions
head(otu_table(ps))
ps_prop <- transform_sample_counts(ps, function (x) x/sum(x))
head(otu_table(ps_prop))

# Export phyloseq data into CSV Files

#Export taxonomy
taxonomy<-tax_table(ps)
write.csv(taxonomy,file="tax.csv")

#Export table of OTUs
table_otu<-otu_table(ps)
write.csv(table_otu,file="otu.csv")

# Export table of OTUs with taxonomy
table_all<-cbind(tax_table(ps),otu_table(ps))
write.csv(table_all,file="otu_tax.csv")
```

## Microbiome package
The microbiome package is built on the phyloseq objects and extends some functions of the phyloseq package in order to facilitate manipulation and processing
microbiome datasets

```r
# install microbiome
remotes::install_github("microbiome/microbiome")
library(microbiome)

# Use the GlobalPatterns datasets included in the phyloseq package to illustrate how to use the microbiome
library(microbiome)
library(phyloseq)
library(knitr)

data(GlobalPatterns)

# Rename GlobalPatterns data (which is a phyloseq object)
physeq <- GlobalPatterns
physeq

# Summarize the Contents of phyloseq Object
summarize_phyloseq(physeq)

# display absolute abundances
head(otu_abs <- abundances(physeq),3)

# Relative abundances
head(otu_rel <- abundances(physeq, "compositional"),3)

# retrieve the total read counts
read_tot <- readcount(physeq)
head(read_tot)

# retrieve the taxonomy table
tax <- tax_table(physeq)
head((tax),3)

# retrieve metadata table
meta <- meta(physeq)
head((meta),3)

# melt the phyloseq data as a data frame table for easier plotting and downstream statistical analysis
# Melt phyloseq data as a data frame table
df <- psmelt(physeq)
kable(head(df,4))

# Number of taxa
num_tax <- ntaxa(physeq)
num_tax

# Most abundant taxa
top10tax <- top_taxa(physeq, n = 10)
top10tax

# Data Transformations
# Transform absolute abundances to relative abundances
head(otu_table(physeq),3)
physeq_comp <- microbiome::transform(physeq, "compositional")
head(otu_table(physeq_comp),3)

# The Hellinger transformation, which is square root of the relative abundance given at the scale [0,1].
physeq_hellinger <- microbiome::transform(physeq_comp, "hellinger")
head(otu_table(physeq_hellinger),3)

# log10 transformation
physeq_log10 <- microbiome::transform(physeq, "log10")
head(otu_table(physeq_log10),3)

# the central log-ratio (CLR) transformation
physeq_clr <- microbiome::transform(physeq, "clr")
head(otu_table(physeq_clr),3)

# Export phyloseq Data into CSV Files

# check files in phyloseq object
slotNames(physeq)

# convert files to data frames
otu_df = as.data.frame(physeq@otu_table)
tax_df = as.data.frame(physeq@tax_table)
sam_df = as.data.frame(physeq@sam_data)

# export otu_table, tax_table, and sam_data using the readr package.
library(readr)

write_csv(otu_df, le = "otu_tab_GlobalPatterns.csv")
write_csv(tax_df, le = "tax_tab_GlobalPatterns.csv")
write_csv(sam_df, le = "sam_tab_GlobalPatterns.csv")
```
