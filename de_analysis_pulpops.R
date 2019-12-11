#creates a series of plots relating to differential expression analysis
# uses the bioconductor package 'limma'

## loading required packages (requires prior installation of each of the following)
require("limma")
require("edgeR")
require("Rsubread")
require("Biobase")
require("gplots")

### set your path to the scripts directory
WD="/home/ssnyde11/scratch/Daphnia_RNAseq_090719/scripts"
setwd(WD)

#obtaining directory paths for both groups
bamDir <- "/home/ssnyde11/scratch/Daphnia_RNAseq_090719/BWAdir"
pulexAnnotation <- "/home/ssnyde11/scratch/genomes/PA42_new_v5.saf"

#obtaining list of file names for both groups
POV_2014_12 <- list.files(bamDir, pattern="\\POV", full.names=TRUE)
NFL3 <- list.files(bamDir, pattern="\\NFL3", full.names=TRUE)
#LPB_2014_32 <- list.files(bamDir, pattern="\\LPB.bam", full.names=TRUE)
#LPA_2014_16 <- list.files(bamDir, pattern="\\LPA201416", full.names=TRUE)
#LPA_2014_32 <- list.files(bamDir, pattern="\\LPA201432", full.names=TRUE)
#KAP_2014_114 <-list.files(bamDir, pattern="\\KAP.bam", full.names=TRUE)

#pulex_pops <- c(POV_2014_12, NFL3, LPB_2014_32, LPA_2014_16, LPA_2014_32, KAP_2014_114)
pulex_pops <- c(NFL3, POV_2014_12)

#creating a count table
pulexpops_fc <- featureCounts(pulex_pops, annot.ext=pulexAnnotation, useMetaFeatures=TRUE, strandSpecific=1, isPairedEnd=FALSE, nthreads=16, isGTFAnnotationFile=FALSE, primaryOnly=TRUE)

save(pulexpops_fc, file="pulex_pops_DE.RData") #saving our featureCounts data to an R binary

## end of read counts section ##

#load("pulex_pops_DE.RData") #starting from an R binary containing the featureCounts list created using the commands above. To run the above commands simply uncomment them (remove the leading '#' from each individual command), and commend out this line.

dge <- DGEList(counts = pulexpops_fc$counts,
               group = c(rep("KAP-2014-16",4), rep("LPA-2014-16",4)),
               genes = pulexpops_fc$annotation$GeneID)

### Now we will apply TMM normalization

dge <- calcNormFactors(dge)

### Let's take a look at the normalization factors

dge$samples

### making a plot of the library counts data
#barplot(dge$samples$lib.size, names=c("KAP-2013-114_1","KAP-2013-114-2","KAP-2013-114-3", "KAP-2013-114-4", "LPA-2014-16-1","LPA-2014-16-2", "LPA-2014-16-3", "LPA-2014-16-4"), las=2, ylim=c(0,30000000))


### making a plot of the counts value
logcounts <- cpm(dge,log=TRUE)
boxplot(logcounts, xlab="", ylab="(log2) counts per million",las=2)
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

# Filtering out genes that have a low number of counts (i.e. are lowly-expressed)

keep <- rowSums(cpm(dge)>1) >= 2
dge <- dge[keep, keep.lib.sizes=FALSE]

# Creating a design matrix to model our experiment

design <- model.matrix(~dge$samples$group)

colnames(design) <- c("KAP-2014-114", "LPA-2014-15")
design #what does this object look like?

#estimate the dispersion

dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

# evaluate the common dispersion

sqrt(dge$common.disp)

# Now we plot the tagwise dispersions against the log2-scaled counts-per million (CPM) values

plotBCV(dge)

# Now we performed the differential expression calculation

fit <- glmFit(dge, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)

summary(de <- decideTestsDGE(lrt, p=0.01, adjust="BH"))
de_tags <- rownames(decideTestsDGE(lrt, p=0.01, adjust="BH"))
de_tags <- rownames(dge)[as.logical(de)]

#mkaing a smear (i.e. a mean-difference) plot of our data

plotSmear(lrt, de.tags=de_tags)
abline(h=c(-2,2), col="blue")

# We can also make the (classic) volcano plot from our data
volcanoData <- cbind(lrt$table$logFC, -log10(lrt$table$PValue))
plot(volcanoData, pch=19)
abline(v=c(-2,2), col="red")

save(dge, file="pulexpopDGE.RData") #saving the updated dge object to our working directory

pulexpops_top_tags <- topTags(lrt, adjust.method="BH", sort.by="PValue", p.value=0.01)
head(pulexpops_top_tags[[1]]) #shows the top results on the screen
write.csv(pulexpops_top_tags[[1]], file="pulexpops_top_tags.csv", row.names=FALSE) #writes a csv file to your working directory
save(pulexpops_top_tags, file= "pulexpops_top_tags.RData") #saves the prist_top_tags file as a p-value

##############
#Making a heatmap with the differentially-expressed genes
library(Biobase) #load this required package if you haven't already done so
library(gplots)
de_data <- dge$counts
colnames(de_data) <- c("KAP-2013-114_1","KAP-2013-114-2","KAP-2013-114-3", "KAP-2013-114-4", "LPA-2014-16-1","L\
PA-2014-16-2", "LPA-2014-16-3", "LPA-2014-16-4")
head(de_data)

top_tags <- topTags(lrt, n= 18146, sort.by="none")

#differential analysis results
de_data <- cbind(de_data, top_tags[[1]])
head(de_data)

diff.genes = rownames(de_data[de_data$FDR<0.01, ])
head(diff.genes)
length(diff.genes)

dge.subset = dge[diff.genes, ]
colnames(dge.subset$counts) <- c("KAP-2013-114_1","KAP-2013-114-2","KAP-2013-114-3", "KAP-2013-114-4", "LPA-2014-16-1","LPA-2014-16-2", "LPA-2014-16-3", "LPA-2014-16-4")
rownames(dge.subset$counts) <- NULL

# plotting the heatmap
heatmap.2(dge.subset$counts,symm=FALSE,symkey=FALSE, scale="row",
          density.info="none",trace="none", key=TRUE,margins=c(10,10))

dev.off()

# plotting and saving the heatmap to a file
pdf("Prist_dge_heatmap.pdf")
heatmap.2(dge.subset$counts,symm=FALSE,symkey=FALSE, scale="row", density.info="none",trace="none",
          key=TRUE,margins=c(10,10))
dev.off()

plotMDS(v, labels=c("KAP-2013-114_1","KAP-2013-114-2","KAP-2013-114-3", "KAP-2013-114-4", "LPA-2014-16-1","L\
PA-2014-16-2", "LPA-2014-16-3", "LPA-2014-16-4"), main="MDS plot for all eight libraries")
#### Done! ######