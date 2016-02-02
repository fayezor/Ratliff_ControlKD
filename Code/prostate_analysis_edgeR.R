#source("https://bioconductor.org/biocLite.R")
library("edgeR")

setwd("/Users/fayezheng/Dropbox/PhD/Side Projects/Prostate Side Project 2015")

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
groups <- ifelse(groups=="C","C","KD")
batch <- substr(samples,2,2)

sample.info <- data.frame(sample=samples, group=groups, batch=batch)

design <- model.matrix(~batch+groups)
rownames(design) <- colnames(d)

###---------- Make into DGEList object and filter ----------###

d <- DGEList(counts=d,group=groups)
d.filt <- d

# Filter out low expression genes, and readjust lib sizes
keep <- which((rowSums(d$counts)/N) >= 5, arr.ind=TRUE) #keep 10854 out of 36135 genes
d.filt <- d.filt[keep, , keep.lib.sizes=FALSE]

###---------- Normalize, estimate disps, and BCV plot ----------###

# Calculate normalization factors to scale the library 
d.filt <- calcNormFactors(d.filt)

# Estimate dispersions
d.filt <- estimateGLMCommonDisp(d.filt, design, verbose=TRUE) 
# Disp = 3.74311, BCV = 1.9347
d.filt <- estimateGLMTrendedDisp(d.filt, design)
d.filt <- estimateGLMTagwiseDisp(d.filt, design)

# MDS plot: labeled/colored by batch
pdf("MDSplot_batch.pdf")
plotMDS(d.filt, main="Multidimensional Scaling Plot of Batches", xlab="LogFC Distance, Dim 1", ylab="LogFC Distance, Dim 2", col=ifelse(batch=="1", "blue", ifelse(batch=="2", "black", "red")), labels=batch, cex=0.75)
legend(x=1.5, y=2.8, legend=c("Batch 1", "Batch 2", "Batch 3"), col=c("blue", "black", "red"), pch = 20, cex=0.75)
dev.off()

# MDS plot: labeled/colored by group
pdf("MDSplot_group.pdf")
plotMDS(d.filt, main="Multidimensional Scaling Plot of Treatments", xlab="LogFC Distance, Dim 1", ylab="LogFC Distance, Dim 2", col=ifelse(groups=="C", "blue", "red"), labels=groups, cex=0.75)
legend(x=1.5, y=2.8, legend=unique(groups), col=c("blue", "red"), pch = 20)
dev.off()

# Plot tagwise dispersions against log2-CPM
pdf("BCVplot.pdf")
plotBCV(d.filt) 
dev.off()

###---------- Test batch effect ----------###

# Fit genewise GLMs
fit <- glmFit(d.filt, design)

# Conduct LRTs
#lrt <- glmLRT(fit, coef=2:3)
#FDR <- p.adjust(lrt$table$PValue, method="BH")
#sum(FDR < 0.5) #3880 genes DE between batches; it's good that we controlled for this

###---------- Test treatment effect ----------###

lrt <- glmLRT(fit)

# Total number of DE genes at FDR<0.05
summary(de <- decideTestsDGE(lrt, p=0.05))

# Get table for only DE genes
table <- topTags(lrt, n=dim(d.filt)[1])$table
table.de <- table[table$FDR<0.05,] #2029 DE genes

# Write to file
write.table(table.de, "DEresults.txt", row.names=TRUE)
write.csv(table.de, "DEresults.csv", row.names=TRUE)
write.table(table, "DEresults_full.txt", row.names=TRUE)
write.csv(table, "DEresults_full.csv", row.names=TRUE)

# Smearplot displays the log-fold changes with the DE genes highlighted:
pdf("smearplot.pdf")
DEnames <- rownames(d.filt)[as.logical(de)]
plotSmear(lrt, de.tags=DEnames)
abline(h = c(-2, 2), col = "blue")
dev.off()

# Gene ontology analysis
#go <- goana(lrt, species="Hs")
#topGO(go,ont="BP",sort="Up",n=30)

