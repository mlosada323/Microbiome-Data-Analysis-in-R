# Demonstration 10
Complete the following demonstration in RStudio. You can follow detail instructions in Xia et al. (2018), Chapter 11: Modeling Over-Dispersed Microbiome Data. All the sections below match the sections in the book

# Modeling Over-Dispersed Microbiome Data
## 11.5. The DESeq2 Package                                              
```r
install.packages("GUniFrac")
library(GUniFrac)
data(throat.otu.tab)
head(throat.otu.tab) 
otu_tab<-throat.otu.tab
head(otu_tab)

# prepare dataset

countData<-as(otu_tab, "matrix")
head(countData)

#DESeq2 need taxa(genes=rows) by samples(=columns)format
countData<-(t(countData))
head(countData)

data(throat.meta)
head(throat.meta)

group<-throat.meta$SmokingStatus
head(group)

metaData<-data.frame(row.names=colnames(countData),group=group)
head(metaData)

# create DESeq2 object and run

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = metaData,
                              design = ~ group)
# filter data
dds <- dds[rowSums(counts(dds)) > 0,]
dds

# Normalize the Count Data 
dds <- estimateSizeFactors(dds) 
sizeFactors(dds)

# Estimate the Dispersion
dds<- estimateDispersions(dds)

# Test the Differential Abundance
dds$group <- relevel(dds$group, "NonSmoker")

#or by

dds$group <- factor(dds$group, levels = c("NonSmoker", "Smoker"))

dds <- DESeq(dds)

# Extract the Results Table
res <- results(dds)
res

mcols(res, use.names=TRUE)
mcols(dds,use.names=TRUE)[1:4,1:4]

substr(names(mcols(dds)),1,10)
head(assays(dds)[["mu"]])
head(dispersions(dds))
head(mcols(dds)$dispersion)

sizeFactors(dds)
head(coef(dds))

# Compare Differential Abundance Between Groups Using Contrast
res <- results(dds, contrast = c("group", "Smoker", "NonSmoker") )
res

# Adjust p-Values Using FDR
sum(res$pvalue < 0.01, na.rm=TRUE )
table(is.na(res$pvalue))

table(res[,"padj"] < 0.1)  
sum(res$padj < 0.1, na.rm=TRUE )

res_Sig <- res[which(res$padj < 0.1 ),]
head(res_Sig[order(res_Sig$log2FoldChange),])
```
![Alt text](Rplot7.png)
```r
tail(res_Sig[order( res_Sig$log2FoldChange ),])

##Diagnostic Plots Using plotMA
plotMA(res)
```
![Alt text](Rplot1.png)
```r
##Diagnostic Plots Using plotDispEsts
plotDispEsts(dds, ylim = c(1e-2, 1e3))
```
![Alt text](Rplot2.png)
```r
##Clustering with Heatmap
rld <- rlog(dds)
vst <-varianceStabilizingTransformation(dds)

par(mfrow = c(1, 3))
plot(log2( 1+counts(dds, normalized=TRUE)[,1:2] ), main="Ordinary log2",col="#00000020", pch=20, cex=0.3 )
plot(assay(rld)[,1:2], main="Regularized-logarithm", col="#00000020", pch=20, cex=0.3 )
plot(assay(vst)[,1:2], main="Variance stabilizing",col="#00000020", pch=20, cex=0.3 )
```
![Alt text](Rplot3.png)
```r
head(assay(rld))[,1-3]

par(mfrow = c(1, 1))

install.packages("gplots")
library("gplots" )
library("RColorBrewer" )
library("genefilter" )
library(SummarizedExperiment)
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 10 )
heatmap.2(assay(rld)[ topVarGenes, ], scale="row",
          trace="none", dendrogram="column",
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
```
![Alt text](Rplot4.png)
```r

##Histogram of p-Values
hist(res$pvalue, breaks=20, col="grey",
     main = "Smoker vs. NonSmoker", xlab = "p-values")
```
![Alt text](Rplot5.png)
```r
##Independent Filtering
metadata(res)

metadata(res)$alpha
metadata(res)$filterThreshold

plot(metadata(res)$filterNumRej,
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)


##Re-estimate the p-Values
res   

# remove filtered out OTUs by independent filtering 
# they have NA adj. pvals
res <- res[ !is.na(res$padj),]
# with NA pvals (outliers)
res <- res[ !is.na(res$pvalue),]

res <- res[, -which(names(res) == "padj")]

install.packages("fdrtool")
library(fdrtool)
res_fdr <- fdrtool(res$stat, statistic= "normal", plot = T)

head(res_fdr)

res_fdr$param[1, "sd"]
sd

res[,"padj"] <- p.adjust(res_fdr$pval, method = "BH")

hist(res_fdr$pval, col = "gray",
     main = "Smoker vs. NonSkoer, correct null model", xlab = "Corrected p-values")
```
![Alt text](Rplot6.png)
```r
# Extract Differentially Abundant OTUs and Export Results Table
table(res[,"padj"] < 0.1)

res[1:2,]
```
## use Negative Binomial distribution on phyloseq object to normalize counts and correct over-dispersion
```r
ps<-readRDS(file="Demo10.RDS")
ps
sample_variables(ps)
head(otu_table(ps))
```
![Alt text](Rplot8.png)
```r
# create DESeQ object
diagdds = phyloseq_to_deseq2(ps, ~region_c) 
# Any variable of the metadata would work to create the DESeQ object

# Calculate geometric means
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate size factors
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

# Get Normalized read counts
normcounts <- counts(diagdds, normalized = TRUE)

# Round read counts
round(normcounts, digits = 0) -> normcountsrd

# Transform matrix of normalized counts to phyloseq object
otu_table(normcountsrd, taxa_are_rows = TRUE) -> ncr

# Replace otu_table in original phyloseq object
otu_table(ps) <- ncr
head(otu_table(ps))
```
![Alt text](Rplot10.png)

Now you can apply all the statistical analyses we have revised in previous chapters