# http://web.stanford.edu/~cengel/cgi-bin/anthrospace/building-my-first-shiny-application-with-ggplot

#source("https://bioconductor.org/biocLite.R")
library("edgeR")
library("ggplot2")
library("lattice")
library("grid")

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
groups <- ifelse(substr(samples,1,1)=="C","Control","Knockdown")
batch <- substr(samples,2,2)

sample.info <- data.frame(sample=samples, group=groups, batch=batch)

###---------- DESCRIPTIVES ----------###

# AR - ENSG00000169083
# PSA (officially named KLK3) - ENSG00000142515
# SULT2B1b - ENSG00000088002
# TNFα - ENSG00000232810
# SREBP1c - ENSG00000072310
# ABCG1 - ENSG00000160179
# NFKB1 - ENSG00000109320

# Just the sult2b1b expression values
sult2b1b.df <- data.frame(sult2b1b=as.numeric(d[rownames(d)=="ENSG00000088002",]), group=groups)
head(sult2b1b)

# Boxplots, but they're useless (min, 1st quartile, median, 3rd quartile all 0)
summary(sult2b1b.df[sult2b1b.df$group=="Control",]$sult2b1b)
summary(sult2b1b.df[sult2b1b.df$group=="Knockdown",]$sult2b1b)

# Summary of the nonzero values - these are a little more interpretable
summary(sult2b1b.df[sult2b1b.df$group=="Control" & sult2b1b.df$sult2b1b!=0,]$sult2b1b)
summary(sult2b1b.df[sult2b1b.df$group=="Knockdown" & sult2b1b.df$sult2b1b!=0,]$sult2b1b)

# Number of 0's in each group, potentially more useful number
C.0s <- length(which(sult2b1b.df[sult2b1b.df$group=="Control",]$sult2b1b==0)) #166 / 209
S.0s <- length(which(sult2b1b.df[sult2b1b.df$group=="Knockdown",]$sult2b1b==0)) #186 / 190

# Table with proportion of zeros in each group
zeros <- data.frame("pZeros"= c(C.0s/length(which(groups=="Control")),S.0s/length(which(groups=="Knockdown"))), "Group" = c("Control","Knockdown"))

###---------- EXPRESSION OF A GENE IN BOTH GROUPS ----------###

# Set seed to reproduce jitter
set.seed(12345)
p <- ggplot(sult2b1b.df, aes(x = group, y = sult2b1b)) + geom_point(position="jitter", aes(color=group))
plot <- p + 
  scale_x_discrete("Group") +
  scale_y_continuous("SULT2B1b Expression") + 
  theme(legend.position="none") + 
  annotate("text", label = paste("% Zeros \n C: ", round(C.0s/length(which(groups=="Control")),2), "%\n KD: ", round(S.0s/length(which(groups=="Knockdown")),2), "%", sep=""), x = 2.6, y = max(sult2b1b.df$sult2b1b), size = 3.5, fontface="bold") + 
  annotate("rect", xmin = 2.4, xmax = 2.8, ymin = max(sult2b1b.df$sult2b1b)-7, ymax = max(sult2b1b.df$sult2b1b)+6, alpha = .2)
plot

# Same plot, with log2 values
set.seed(12345)
p <- ggplot(sult2b1b.df, aes(x = group, y = log2(sult2b1b+1))) + geom_point(position="jitter", aes(color=group))
plot <- p + 
  scale_x_discrete("Group") +
  scale_y_continuous("SULT2B1b log2(expression)") + 
  theme(legend.position="none") + 
  annotate("text", label = paste("% Zeros \n C: ", round(C.0s/length(which(groups=="Control")),2), "%\n KD: ", round(S.0s/length(which(groups=="Knockdown")),2), "%", sep=""), x = 2.6, y = max(log2(sult2b1b.df$sult2b1b)), size = 3.5, fontface="bold") + 
  annotate("rect", xmin = 2.4, xmax = 2.8, ymin = max(log2(sult2b1b.df$sult2b1b))-0.8, ymax = max(log2(sult2b1b.df$sult2b1b))+0.8, alpha = .2)
plot

#

###---------- EXPRESSION OF A GENE, THRESHOLDED TOP/BOTTOM ----------###

sult2b1b.df <- data.frame(sult2b1b=as.numeric(d[rownames(d)=="ENSG00000088002",]), x = rep(1, length(groups)), group=groups)

# Given low and high threshold, calculate % C's in each group
low=1 
high=40 
high.df = sult2b1b.df[sult2b1b.df$sult2b1b>high,]
low.df = sult2b1b.df[sult2b1b.df$sult2b1b<low,]
high.n = dim(high.df)[1]
low.n = dim(low.df)[1]
high.C = length(which(high.df[high.df$group=="Control",]$group=="Control"))
low.C = length(which(low.df[low.df$group=="Control",]$group=="Control"))

set.seed(37291)
p <- ggplot(sult2b1b.df, aes(x = x, y = sult2b1b)) + geom_jitter(position=position_jitter(width=0.05), aes(color=group)) +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), legend.title=element_blank()) +
  scale_y_continuous("SULT2B1b Expression", breaks=sort(c(seq(0, max(sult2b1b.df$sult2b1b), by=10), low,high)) ) +
  geom_hline(yintercept=low, linetype="dashed") + 
  geom_hline(yintercept=high, linetype="dashed") + 
  annotate("text", label = paste("High Group Cells \n C: ", round(high.C*100/high.n,2), "% \n KD: ", round((high.n-high.C)*100/high.n,2), "%", sep=""), x = 1.1, y = high+8, size = 3.5, fontface="bold") + 
  annotate("rect", xmin = 1.06, xmax = 1.14, ymin = high+2, ymax = high+13, alpha = .2) + 
  annotate("text", label = paste("Low Group Cells \n C: ", round(low.C*100/low.n,2), "% \n KD: ", round((low.n-low.C)*100/low.n,2), "%", sep=""), x = 1.1, y = low+8, size = 3.5, fontface="bold") + 
  annotate("rect", xmin = 1.06, xmax = 1.14, ymin = low+2, ymax = low+13, alpha = .2)
p

###---------- CORRELATION BETWEEN TWO GENES ----------###

# AR - ENSG00000169083
# PSA (officially named KLK3) - ENSG00000142515
# SULT2B1b - ENSG00000088002
# TNFα - ENSG00000232810
# SREBP1c - ENSG00000072310
# ABCG1 - ENSG00000160179
# NFKB1 - ENSG00000109320

sult2b1b.AR.df <- data.frame(sult2b1b=as.numeric(d[rownames(d)=="ENSG00000088002",]), AR=as.numeric(d[rownames(d)=="ENSG00000169083",]), group=groups)

set.seed(5749)
p <- ggplot(sult2b1b.AR.df, aes(x = sult2b1b, y = AR)) + geom_point(position="jitter", aes(color=group)) + 
  scale_x_continuous("SULT2B1b Expression") +
  scale_y_continuous("AR Expression") + 
  theme(legend.title=element_blank()) + 
  annotate("text", label = paste("Correlation:", round(cor(sult2b1b.AR.df$sult2b1b, sult2b1b.AR.df$AR),2)), x = 45, y = 1750, size = 3.5, fontface="bold")
p
