#source("https://bioconductor.org/biocLite.R")
library("DESeq2")
library("edgeR")

setwd("/Users/fayezheng/Dropbox/Side Projects/Prostate Side Project 2015")

###---------- Getting dataset ready ----------###

d.full <- read.csv("GeneCountMatrix.AllBatches_withannot_.csv")
rownames(d.full) <- d.full[,1]
d.full <- d.full[,-c(1,2)]

samples.full <- names(d.full)

# Remove the control and repeat samples (for now)
ctrl.idx <- which(substr(samples.full, nchar(samples.full)-3+1, nchar(samples.full))=="trl", arr.ind=TRUE)
repeat.idx <- which(substr(samples.full, nchar(samples.full)-3+1, nchar(samples.full))=="eat", arr.ind=TRUE)
remove.idx <- c(ctrl.idx, repeat.idx)
ctrls <- d.full[, ctrl.idx]
repeats <- d.full[, repeat.idx]

d <- d.full[, -remove.idx] # 36135 genes, 399 samples
G <- dim(d)[1]
N <- dim(d)[2]

samples <- samples.full[-remove.idx]
groups <- substr(samples,1,1)
batch <- substr(samples,2,2)

sample.info <- data.frame(sample=samples, group=groups, batch=batch)

design <- model.matrix(~batch+groups)
rownames(design) <- colnames(d)

###---------- Make into DESeq object and filter ----------###

dds <- DESeqDataSetFromMatrix(countData = d,
                            colData = sample.info,
                            design = ~group)
dds
design(dds)

###---------- Differential expression testing ----------###

# Need to do the testing in parallel, otherwise super slow
library("BiocParallel")
register(MulticoreParam(8)) #register 4 cores

# In one function, estimates size factors, dispersions, gene-wise dispersions
dds <- DESeq(dds, parallel=TRUE)
res <- results(dds, parallel=TRUE, alpha=0.05)
res

# check out different contrasts
resMF.batch12 <- results(ddsMF, contrast=c("batch","1","2"))
resMF.batch23 <- results(ddsMF, contrast=c("batch","2","3"))
resMF.batch13 <- results(ddsMF, contrast=c("batch","1","3"))
resMF.group <- results(ddsMF, contrast=c("group","C","S"))

# Order results table by smallest adjusted p-value
resOrdered <- res[order(res$padj),]
head(resOrdered)

# Subsetting results for adjusted p-value < 0.05
resSig <- subset(resOrdered, padj < 0.05)

# Write results to file
write.csv(as.data.frame(resSig), file="DESeq2results.csv")
write.table(as.data.frame(resSig), file="DESeq2results.txt")

# Summarize basic tallies
summary(res)

# More info on results columns
mcols(res)$description

###---------- Visualizations ----------###

# MA plot shows the log2 fold changes attributable to a given variable over the mean of normalized counts.
# Points will be colored red if the adjusted p value is less than 0.05
plotMA(res, main="DESeq2", ylim=c(-2,2))

# Plot of counts for smallest p-value gene across samples
g <- plotCounts(dds, gene=which.min(res$padj), intgroup="group", returnData=TRUE)
library("ggplot2")
ggplot(g, aes(x=condition, y=count)) +
  geom_point(position=position_jitter(w=0.1,h=0)) +
  scale_y_log10(breaks=c(25,100,400))


###---------- PCA plot ----------###

# Transform the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size
rld <- rlog(dds)

# Overall PCA plot with all combinations of batch x group
plotPCA(rld, intgroup=c("batch", "group"))

# PCA plot with points colored for batch
plotPCA(rld, intgroup="batch")

# PCA plot with points colored for group
plotPCA(rld, intgroup="group")

# Example of how to customize the PCA plot (the overall plot)
data <- plotPCA(rld, intgroup=c("condition", "type"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))







