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

#save unclustered significance results
sig.mot.unclust <- sig.mot

sig.mot <- sig.mot[!duplicated(sig.mot$clust),]

#save initial deviation and variability results
variability.orig <- variability
dev.orig <- dev

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

####Identify candidate targets of injury responsive motifs (0-3hpa)####

#pull all motifs associated with the injury responsive group of motifs
#recreate clustering used for heatmap
scoreDist <- dist(getScores, method = "euclidean")
scoreClust <- hclust(scoreDist, method = "ward.D2")

#identify cluster associated with increase from 0 to 3hpa
plot(scoreClust, labels = colnames(scoreDist))
abline(a = 20, b = 0, col = "red")
clusts <- cutree(scoreClust, h = 20)
clusts <- as.data.frame(clusts)

#pull 0-3hpa cluster (1)
earlyI <- rownames(clusts)[clusts$clusts == 1]

#for each representative motif, pull all motifs that are in the same cluster
earlyI <- motClust[motClust$clust %in% motClust[motClust$ID %in% earlyI,"clust"],]

#verify that all of the motifs we consider show meaningful variability and increase in accessibility from 0 to 3 hpa
getScores <- as.data.frame(deviationScores(dev.orig))
getScores <- getScores[rownames(getScores) %in% rownames(variability.orig),]

#average across replicates
getScores <- as.data.frame(t(getScores))

getScores$ID <- as.factor(sub("\\d*_final_shift.bam","",rownames(getScores)))

getScores <- aggregate(getScores[,1:(ncol(getScores)-1)], by = list(getScores$ID), FUN = mean)

rownames(getScores) <- getScores$Group.1
getScores$Group.1 <- NULL

getScores <- as.data.frame(t(getScores))

getScores <- getScores[,c(1,5,7,3,2,6,8,4)]

#calculate the change in accessibility from 0 to 3 hpa
#scale it to the overall range of accessibility for that motif
change.range <- apply(getScores, 1, function(x)max(x) - min(x))
change.H <- (getScores$`3H` - getScores$`0H`)/change.range
change.F <- (getScores$`3F` - getScores$`0F`)/change.range

#look for motifs that incraese from 0 to 3 by a factor greater than 25% 
#of the largest change in accessibility for that motif
change.H <- change.H[change.H > 0.25]
change.F <- change.H[change.H > 0.25]

change.B <- unique(c(names(change.H),names(change.F)))

#filter earlyI to include only injury response cluster motifs
#that actually show an increase in the chromVar analysis
earlyI <- earlyI[earlyI$ID %in% change.B,"ID"]

#reformat to match target filenames
earlyI <- gsub("/Homer.*",".txt",earlyI)
earlyI <- gsub("/","_",earlyI)

unlink("resources/chromVar_HOMER_sub/injuryClusts", recursive = T)
dir.create("resources/chromVar_HOMER_sub/injuryClusts", showWarnings = F)

#generate individual HOMER motif files from multi-motif file
system('cd resources/chromVar_HOMER_sub;bash homerSplit.sh')

#concatenate individual motif files for a given cluster
command <- paste(earlyI,collapse = " ")
command <- paste0("cd resources/chromVar_HOMER_sub; cat ", command, " > injuryClusts/earlyI.motif")
command <- gsub("([(){}:])","\\\\\\1",command)
system(command)

#find all hits for the injury response motifs genome wide using homer
system("cd resources/chromVar_HOMER_sub; ./pullMotifHits.sh")

#read in peaks with injury responsive motifs
injPeaks <- read.table("resources/chromVar_HOMER_sub/injuryClusts/earlyI.hits.peaks.bed", sep = "\t")

#bring in uropa annotation info
peakLink <- read.table("resources/untreated_consensus_finalhits.txt", header = T, stringsAsFactors = F)

annot <- read.csv("resources/expanded_dovetail_SP_annot.csv", row.names = 1)

injPeaks$ID <- mapvalues(injPeaks$V5, from = peakLink$peak_id, to = peakLink$ID, warn_missing = F)

#drop peaks without gene annotation
injPeaks <- injPeaks[grepl("^g",injPeaks$ID),]

injPeaks <- merge(injPeaks, annot, by = "ID", all.x = T)

#bring in RNA-seq data to only focus on genes with an increase in expression from 0 to 3 hpa
load("Analysis_Output/RNA/RNA_DGE.RData")

injPeaks <- injPeaks[injPeaks$ID %in% c(H3vH0.DG.up$ID, F3vF0.DG.up$ID),]

#bring in ATAC-seq data to only focus on peaks with an increase in expression from 0 to 3 hpa
load("Analysis_Output/ATAC/ATAC_DGE.RData")

injPeaks <- injPeaks[injPeaks$V5 %in% c(H3H0.DG.up$ID, F3F0.DG.up$ID),]

#bring in info on Wnt pathway genes
wntGenes <- read.table("Wnt_Survey/wntGeneIDs.txt", stringsAsFactors = F)[,1]

#focus on Wnt pathway genes that are up from 0 to 3
wntGenes <- injPeaks[injPeaks$ID %in% wntGenes,]

#import motif enrichment data for 3hpa
enrichment.3 <- read.table("Analysis_Output/ATAC/Up3Enrichment/knownResults.txt", sep = "\t", skip = 1)

#reformat motif name
enrichment.3$V1 <- gsub("[/].*","",enrichment.3$V1)

#drop all hits above an FDR of 0.05
enrichment.3 <- enrichment.3[enrichment.3$V5 <= 0.05,]

injIndex <- logical(nrow(wntGenes))
eachMotifRes <- character(0)

for (i in 1:nrow(wntGenes)) {
  res <- unlist(strsplit(wntGenes[i,5], split = ","))
  resV <- character(0)
  for (j in res) {
    if(j %in% enrichment.3$V1) {
      injIndex[i] <- T
      resV <- c(resV, j)
    }
  }
  if(length(resV) != 0) {
    eachMotifRes <- c(eachMotifRes,paste(resV, collapse = ","))
  }
}

wntGenes <- wntGenes[injIndex,]
wntGenes$enrichedMotifs <- eachMotifRes

#look for candidate regulators based on interpro annotations of genes that increase from 0-3hpa
#the interproscan results were generated by running interProScan on the Dovetail protein models using standard settings.
ipro <-read.csv("resources/Dovetail_Interproscan.csv")

searchTerms <- paste(c("homeo","nuclear receptor", "c2h2", "ets-domain", "bhlh","bzip","leucine zipper","irf","duf3446","EGR_N","EGR_C"), collapse = "|")

tf.index <- grepl(searchTerms, apply(ipro,1,function(x) paste(x, collapse = ", ")), ignore.case = T)

ipro.tf <- ipro[tf.index,]

ipro.tf <- ipro.tf[ipro.tf$ID %in% c(H3vH0.DG.up$ID, F3vF0.DG.up$ID),]

ipro.tf <- merge(ipro.tf, annot, by = "ID", all.x = T)

