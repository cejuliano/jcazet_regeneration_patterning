# this script will use the pre-generated chromVar read counts object and find 
# HOMER motifs associated with significant variability

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(chromVAR)
library(motifmatchr)
library(GenomicRanges)
library(BiocParallel)
library(BSgenome.Hymag.NIH.Hymag2)
library(Biostrings)
library(SummarizedExperiment)
library(chromVARmotifs)
library(ggplot2)
library(plyr)
library(gplots)

RNGkind(sample.kind = "default")

# load pre-calculated fragment counts data
load("resources/untreated_chromVar_Counts.RData")

# load HOMER PWMs
data("homer_pwms")

# find PWM hits in the genome
set.seed(123456)
Motifs <- matchMotifs(homer_pwms, fragment_counts, genome = BSgenome.Hymag.NIH.Hymag2)

set.seed(123456)
dev <- computeDeviations(object = fragment_counts, annotations = Motifs)
variability <- computeVariability(dev)

#save deviation scores
write.csv(deviationScores(dev), file = "Analysis_Output/ATAC/untreated_deviations.csv")

#set a variability cutoff of 2
variability <- variability[variability$variability >= 2,]

#load motif cluster information
motClust <- read.csv("resources/HOMER_motif_clusters.csv", row.names = 1)

#select a representative motif from each cluster based on the one that best accounts for
#variability across sample/treatment conditions

variability$ID <- rownames(variability)

difdev <- differentialDeviations(dev, "sampleType")

significanceTest <- difdev
significanceTest$ID <- rownames(significanceTest)

difdev <- differentialDeviations(dev, "structure")
difdev$ID <- rownames(difdev)

significanceTest <- merge(significanceTest[,c(2,3)], difdev[,c(2,3)], by = "ID", suffixes = c("_type","_structure"))

difdev <- differentialDeviations(dev, "timePoint")
difdev$ID <- rownames(difdev)

significanceTest <- merge(significanceTest, difdev[,c(2,3)], by = "ID")

significanceTest <- merge(significanceTest, motClust, by = "ID")

significanceTest$minP <- apply(significanceTest[,2:4],1,min)

sig.mot <- merge(variability,significanceTest, by = "ID")

sig.mot <- sig.mot[order(sig.mot$minP),]

sig.mot <- sig.mot[!duplicated(sig.mot$clust),]

#redo dev and variability calc using only the nr motif list
homer_pwms <- homer_pwms[which(names(homer_pwms) %in% sig.mot$ID)]

set.seed(123456)
Motifs <- matchMotifs(homer_pwms, fragment_counts, genome = BSgenome.Hymag.NIH.Hymag2)

set.seed(123456)
dev <- computeDeviations(object = fragment_counts, annotations = Motifs)
variability <- computeVariability(dev)

write.csv(variability, file = "Analysis_Output/ATAC/HOMER_chromVar_variability.csv")

motifAccessibility <- function(x) {
  
  getScores <- as.data.frame(deviationScores(dev))
  # IDs <- strsplit(rownames(getScores), split = "/")
  # IDs <- vapply(IDs, function(x) x[[1]],"")
  # getScores$ID <- IDs
  # rm(IDs)
  
  getScores <- as.data.frame(t(getScores))
  
  getScores$sample <- as.factor(gsub("\\d_final_shift.bam","",rownames(getScores)))
  
  getScores.mean <- aggregate(getScores[,1:(ncol(getScores) - 1)], by = list(getScores$sample), FUN=mean)
  
  getScores.mean$type <- gsub("^\\d*","",getScores.mean$Group.1)
  getScores.mean$time <- as.factor(as.numeric(gsub("[HF].*","",getScores.mean$Group.1)))
  getScores.mean$Group.1 <- NULL
  
  x <- gsub("[(]","[(]",x)
  x <- gsub("[)]","[)]",x)
  x <- gsub("[+].*$","",x)
  x <- gsub("[,].*$","",x)
  x <- gsub("[?].*$","",x)
  
  
  print(x)
  
  
  motifPlot <- colnames(getScores.mean)[grep(x,colnames(getScores.mean),ignore.case = T)][1]
  
  getScores.mean <- getScores.mean[,c(motifPlot,"type","time")]
  
  colnames(getScores.mean) <- c("Accessibility","type","time")
  
  getScores.all <- getScores
  getScores.all$type <- gsub("^\\d*","",rownames(getScores.all))
  getScores.all$type <- as.factor(gsub("\\d_final_shift.bam","",getScores.all$type))
  
  getScores.all$time <- as.factor(gsub("[HF].*","",getScores.all$sample))

  getScores.all <- getScores.all[,c(motifPlot,"type","time")]
  
  colnames(getScores.all) <- c("Accessibility","type","time")
  
  gg <- ggplot(data = getScores.mean, aes(x = time, y = Accessibility, group = type, color = type)) + geom_line(size = 1.5)
  gg <- gg + theme_bw()
  gg <- gg + geom_point(data = getScores.all, aes(x = time, y = Accessibility, group = type, color = type), size = 2, alpha = 0.75)
  gg <- gg + scale_color_manual(values=c("blue", "red"))
  gg
  
  ggsave(paste0("plots/ATAC/HOMER_motifs/chromVar_",x,"_accessibility_HOMER.pdf"), width = 5, height = 5)
}  

unlink("plots/ATAC/HOMER_motifs", recursive = T)
dir.create("plots/ATAC/HOMER_motifs", showWarnings = F)

lapply(as.character(variability[variability$variability >= 2, "name"]), function(x) motifAccessibility(x))

#heatmap of motif accessibility changes
getScores <- as.data.frame(deviationScores(dev))
getScores <- getScores[rownames(getScores) %in% sig.mot$ID,]

#average across replicates
getScores <- as.data.frame(t(getScores))

getScores$ID <- as.factor(sub("\\d*_final_shift.bam","",rownames(getScores)))

getScores <- aggregate(getScores[,1:(ncol(getScores)-1)], by = list(getScores$ID), FUN = mean)

rownames(getScores) <- getScores$Group.1
getScores$Group.1 <- NULL

getScores <- as.data.frame(t(getScores))

getScores <- getScores[,c(1,5,7,3,2,6,8,4)]

getPalette = colorRampPalette(c("white","red"))

pdf(file = "plots/ATAC/HOMER_motifs/HOMER_motif_heatmap.pdf", height = 8, width = 5)
heatmap.2(as.matrix(getScores), Rowv = T, Colv = F, dendrogram = "none", scale = "row",
          hclustfun = function(x) hclust(x, method = "ward.D2"),
          distfun = function(x) dist(x, method = "euclidean"),
          col = getPalette(90), trace = "none", key = FALSE, sepcolor = "white",
          keysize = 0.1, margins = c(8,12), colsep = c(4),
          labRow = sub("[()].*$","",rownames(getScores)), labCol = sub("_final.bam","",colnames(getScores)))
dev.off()

