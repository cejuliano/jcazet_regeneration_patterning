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
library(gridExtra)
library(RColorBrewer)
library(plyr)
library(ggmotif)
library(TFBSTools)
library(gplots)

RNGkind(sample.kind = "default")

# load pre-calculated fragment counts data
load("resources/chromVar_iCounts.RData")

# load HOMER PWMs
data("homer_pwms")

# find PWM hits in the genome
set.seed(123456)
Motifs <- matchMotifs(homer_pwms, fragment_counts, genome = BSgenome.Hymag.NIH.Hymag2)

set.seed(123456)
dev <- computeDeviations(object = fragment_counts, annotations = Motifs)
variability <- computeVariability(dev)

#save deviation scores
write.csv(deviationScores(dev), file = "Analysis_Output/ATAC/full_deviations.csv")

#focus on the injury responsive TFBMs identified in the untreated dataset
motifList <- list.files("plots/ATAC/HOMER_motifs", pattern = "^chromVar_")
motifList <- gsub("\\[","",motifList)
motifList <- gsub("\\]","",motifList)
motifList <- gsub("chromVar_","",motifList)
motifList <- gsub("chromVar_","",motifList)
motifList <- gsub("_accessibility_HOMER.pdf","",motifList)

motifAccessibility <- function(x) {
  
  getScores <- as.data.frame(deviationScores(dev))
  # IDs <- strsplit(rownames(getScores), split = "/")
  # IDs <- vapply(IDs, function(x) x[[1]],"")
  # getScores$ID <- IDs
  # rm(IDs)
  
  getScores <- as.data.frame(t(getScores))
  
  getScores$sample <- as.factor(gsub("\\d_final_shift.bam","",rownames(getScores)))
  
  getScores.mean <- aggregate(getScores[,1:(ncol(getScores) - 1)], by = list(getScores$sample), FUN=mean)
  
  getScores.mean$type <- gsub("^\\d+","",getScores.mean$Group.1)
  getScores.mean$time <- as.factor(as.numeric(gsub("i?[HF]","",getScores.mean$Group.1)))
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
  
  getScores.mean$treatment <- grepl("i",getScores.mean$type)
  
  getScores.mean$treatment[getScores.mean$treatment == T] <- "iCRT14"
  getScores.mean$treatment[getScores.mean$treatment == F] <- "Untreated"
  
  getScores.mean$treatment <- factor(getScores.mean$treatment, levels = c("Untreated","iCRT14"))
  
  
  getScores.all <- getScores
  getScores.all$type <- gsub("^\\d+","",rownames(getScores.all))
  getScores.all$type <- as.factor(gsub("\\d_final_shift.bam","",getScores.all$type))
  
  getScores.all$time <- as.factor(gsub("i?[HF].*","",getScores.all$sample))

  getScores.all <- getScores.all[,c(motifPlot,"type","time")]
  
  colnames(getScores.all) <- c("Accessibility","type","time")
  
  getScores.all$treatment <- grepl("i",getScores.all$type)
  
  getScores.all$treatment[getScores.all$treatment == T] <- "iCRT14"
  getScores.all$treatment[getScores.all$treatment == F] <- "Untreated"
  
  getScores.all$treatment <- factor(getScores.all$treatment, levels = c("Untreated","iCRT14"))
  
  getScores.mean$type <- gsub("i","",getScores.mean$type)
  getScores.all$type <- gsub("i","",getScores.all$type)
  
  getScores.mean$treatment <- paste(getScores.mean$type, getScores.mean$treatment, sep = "_")
  getScores.all$treatment <- paste(getScores.all$type, getScores.all$treatment, sep = "_")
  
  gg <- ggplot(data = getScores.mean, aes(x = time, y = Accessibility, group = treatment, color = treatment)) + geom_line(size = 1.5)
  gg <- gg + theme_bw()
  gg <- gg + geom_point(data = getScores.all, aes(x = time, y = Accessibility, group = treatment, color = treatment), size = 2, alpha = 0.75)
  gg <- gg + scale_color_manual(values=c("turquoise", "blue", "purple", "red"))
  gg <- gg + facet_grid(. ~ type)
  gg <- gg + theme(legend.position = "none")
  gg
  
  ggsave(paste0("plots/ATAC/HOMER_motifs_icrt/chromVar_",x,"_accessibility_HOMER.pdf"), width = 9, height = 5)
}  

unlink("plots/ATAC/HOMER_motifs_icrt", recursive = T)
dir.create("plots/ATAC/HOMER_motifs_icrt", showWarnings = F)

lapply(motifList, function(x) motifAccessibility(x))


