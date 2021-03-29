#This script will generate plots based on the results from the edgeR DGE of the ATAC-seq dataset

####Setup####

library(rstudioapi)
library(ggplot2)
library(Seurat)
library(gplots)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(plyr)

options(stringsAsFactors = F)
setwd(dirname(getActiveDocumentContext()$path))

#load in the edgeR DGE results
load("Analysis_Output/ATAC/ATAC_DGE.RData")
load("Analysis_Output/ATAC/iATAC_DGE.RData")

#load in scRNA-seq data
ds <- readRDS("resources/Hydra_Seurat_Whole_Genome_updated.rds")

hFind <-function (y,x) { 
  return (y@assays$RNA@data@Dimnames[[1]][grep(x,y@assays$RNA@data@Dimnames[[1]],ignore.case = T)])
}

fcCompare <-function(FR,HR,intDFs) {
  
  fc.comp <- merge(FR, HR, by = "ID", suffixes = c("_F", "_H"))
  
  intIDs <- lapply(intDFs, function(x) x$ID)
  intIDs <- do.call(c, intIDs)
  
  fc.comp <- fc.comp[fc.comp$ID %in% intIDs,]
  
  fc.comp$Interesting <- "Not Interesting"
  fc.comp$Interesting[fc.comp$ID %in% intDFs[[length(intDFs)]]$ID] <- "Interesting"
  
  for (i in 1:nrow(fc.comp)) {
    if(fc.comp[i,"logFC_F"] > fc.comp[i,"logFC_H"] & fc.comp[i,"Interesting"] == "Interesting") {
      fc.comp[i,"Interesting"] <- "Foot"
    } else if (fc.comp[i,"logFC_F"] < fc.comp[i,"logFC_H"] & fc.comp[i,"Interesting"] == "Interesting"){
      fc.comp[i,"Interesting"] <- "Head"
    }
  }
  
  return(fc.comp)
}

scatterTheme <- theme_bw() + 
  theme(legend.position="none") +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(axis.title = element_blank()) +
  theme(axis.text = element_blank())

scatterPlot <- function(df, plotName, x = "logFC_F", y = "logFC_H", colorBy = "Interesting", corPlot = T) {
  
  #set the colors and sizes
  if(length(unique(df[,colorBy])) == 3) {
    colors.use <- c("#0000FF","#FF0000","darkgrey")
    sizes.use <- c(2,2,1)
  } else if(length(unique(df[,colorBy])) == 2) {
    colors.use <- c("#FF0000","darkgrey")
    sizes.use <- c(2,1)
  } else {
    colors.use <- c("darkgrey")
    sizes.use <- c(1)
  }
  
  #set the axis limits
  axis.ulimit <- ceiling(max(c(df[,x],df[,y])))
  axis.llimit <- floor(min(c(df[,x],df[,y])))
  
  #initialize ggplot object and set basic aesthetics
  gg <- ggplot(df) 
  gg <- gg + scatterTheme
  
  #make custom axis lines that run through the center of the plot
  gg <- gg + geom_abline(intercept = 0, slope = 0, size = 1)
  gg <- gg + geom_vline(xintercept = 0, size = 1)
  
  
  #add scatterplot
  gg <- gg + geom_point(data = df, aes(x = !!sym(x), y = !!sym(y), color = !!sym(colorBy), size = !!sym(colorBy), name = ID))
  gg <- gg + scale_color_manual(values = colors.use)
  gg <- gg + scale_size_manual(values = sizes.use)
  
  if(corPlot) {
    gg <- gg + geom_abline(intercept = 0, slope = 1, linetype = "longdash", size = 0.8)
    gg <- gg + scale_y_continuous(limits = c(axis.llimit,axis.ulimit))
    gg <- gg + scale_x_continuous(limits = c(axis.llimit,axis.ulimit))
  }
  
  gg
  
  ggsave(filename = paste0("plots/ATAC/",plotName,".pdf"), plot = gg, width = 4, height = 4, useDingbats=FALSE)
  
  return(gg)
}

getDropIDs <- function(x) {
  ds.IDs <- lapply(x, function(y) hFind(ds,paste0(y,".t")))
  
  ds.IDs <- Filter(length, ds.IDs)  
  ds.IDs <- unlist(ds.IDs)  
  ds.IDs <- unique(ds.IDs)
  return(ds.IDs)
}

strEnrichment <- function(x) {
  IDs <- getDropIDs(x)
  
  enrichmentRes <- FindMarkers(ds, ident.1 = c("ecEp_head/hypostome","ecEp_battery(mp)","enEp_head","enEp_tentacle"),
                               ident.2 = c("enEp_foot","ecEp_basal_disk","ecEp_peduncle"),
                               features = IDs, logfc.threshold = 0, min.pct = 0, min.cells.group = 0, test.use = "wilcox")
  
  enrichmentRes$ID <- gsub("[.]t.*","",rownames(enrichmentRes))
  
  enrichmentRes <- enrichmentRes[order(enrichmentRes$p_val_adj),]
  enrichmentRes <- enrichmentRes[!duplicated(enrichmentRes$ID),]  
  
  # partition into head and foot markers
  enrichmentRes.H <- enrichmentRes[enrichmentRes$avg_logFC > 0 & enrichmentRes$p_val_adj < 1e-6,]
  enrichmentRes.F <- enrichmentRes[enrichmentRes$avg_logFC < 0 & enrichmentRes$p_val_adj < 1e-6,]
  
  structureRes <- rep("NoEnrichment", length(x))
  structureRes[x %in% enrichmentRes.H$ID] <- "Head"
  structureRes[x %in% enrichmentRes.F$ID] <- "Foot"
  
  structureRes <- data.frame(ID = x,
                             structureRes = structureRes, 
                             strFC = mapvalues(x, from = enrichmentRes$ID, to = enrichmentRes$avg_logFC, warn_missing = F))
  
  structureRes$strFC <- as.numeric(gsub("^g.*",NA, structureRes$strFC))
  
  return(list(structureRes,enrichmentRes))
  #return(structureRes)
}

plotATAC <- function(chr,start,stop) {
  
  ROI <- GRanges(seqnames = chr,
                 ranges = IRanges(start = start, end = stop))
  
  IDs <- c("0H","3H","8H",
           "0F","3F","8F")
  
  regType <- c(rep("H",3), rep("F", 3))
  
  files <- paste0(IDs, "_MG_final_shift.bw")
  
  bigWigs <- data.frame(ID = IDs,
                        files = files,
                        type = regType)
  
  bwValues <- data.frame(position = numeric(0), scores = numeric(0), 
                         ID = character(0), type = character(0))
  
  for (i in 1:nrow(bigWigs)) {
    
    path <- paste0("resources/bigwigs/", bigWigs[i,"files"])
    bwf <- BigWigFile(path)
    BW <- import(bwf, which = ROI)
    
    BW.df <- data.frame(position = start(BW), scores = score(BW), ID = bigWigs[i,"ID"], 
                        type = bigWigs[i,"type"])
    bwValues <- rbind(bwValues, BW.df)
    
  }
  
  
  gg <- ggplot(bwValues, aes(x = position, y = scores, color = ID))
  gg <- gg + geom_line(size = 2)
  gg <- gg + theme_bw() + facet_grid(. ~ type)
  gg <- gg + scale_color_manual(values = c("#dbdbff","#ffdbdb",
                                           "#8484ff","#ff8484",
                                           "#0000ff","#ff0000"))
  gg <- gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  gg <- gg + theme(legend.position="none")
  gg
}

####export normalized counts matrix####
inormalizedCounts$ID <- rownames(inormalizedCounts)
normalizedCounts$ID <- rownames(normalizedCounts)
normalizedCounts <- merge(normalizedCounts,inormalizedCounts, by = "ID")

rownames(normalizedCounts) <- normalizedCounts$ID

normalizedCounts$ID <- NULL

write.csv(normalizedCounts, file = "Analysis_Output/ATAC/fullNormCounts.csv")

#need to reload results
load("Analysis_Output/ATAC/ATAC_DGE.RData")
load("Analysis_Output/ATAC/iATAC_DGE.RData")

####Visualize differentially activated peaks####

FC.comp.3 <- fcCompare(F3F0,H3H0,list(F3F0.DG,H3H0.DG,HRFR3c0.DG))

FC.comp.3.exp <- FC.comp.3[,c(1,6,10,2,16,30)]
colnames(FC.comp.3.exp) <- c("Peak ID","Gene ID", "Swissprot Annotation",
                             "Foot Reg logFC","Head Reg logFC","Structural Enrichment")
write.csv(FC.comp.3.exp, file = "Analysis_Output/ATAC/fc.comp.3.csv", row.names = F)

scatterPlot(FC.comp.3, "FC_Comp_3")

FC.comp.8 <- fcCompare(F8F0,H8H0,list(F8F0.DG,H8H0.DG,HRFR8c0.DG))

FC.comp.8.exp <- FC.comp.8[,c(1,6,10,2,16,30)]
colnames(FC.comp.8.exp) <- c("Peak ID","Gene ID", "Swissprot Annotation",
                             "Foot Reg logFC","Head Reg logFC","Structural Enrichment")
write.csv(FC.comp.8.exp, file = "Analysis_Output/ATAC/fc.comp.8.csv", row.names = F)

scatterPlot(FC.comp.8, "FC_Comp_8")

FC.comp.8i <- fcCompare(Fi8Fi0,Hi8Hi0,list(Fi8Fi0.DG,Hi8Hi0.DG,inhib8c0.DG))

FC.comp.8i.exp <- FC.comp.8i[,c(1,6,10,2,16,30)]
colnames(FC.comp.8i.exp) <- c("Peak ID","Gene ID", "Swissprot Annotation",
                             "Foot Reg logFC","Head Reg logFC","Structural Enrichment")
write.csv(FC.comp.8i.exp, file = "Analysis_Output/ATAC/fc.comp.8i.csv", row.names = F)

scatterPlot(FC.comp.8i, "FC_Comp_8i")

FC.comp.12 <- fcCompare(F12F0,H12H0,list(F12F0.DG,H12H0.DG,HRFR12c0.DG))

FC.comp.12.exp <- FC.comp.12[,c(1,6,10,2,16,30)]
colnames(FC.comp.12.exp) <- c("Peak ID","Gene ID", "Swissprot Annotation",
                              "Foot Reg logFC","Head Reg logFC","Structural Enrichment")
write.csv(FC.comp.12.exp, file = "Analysis_Output/ATAC/fc.comp.12.csv", row.names = F)

scatterPlot(FC.comp.12, "FC_Comp_12")

FC.comp.12i <- fcCompare(Fi12Fi0,Hi12Hi0,list(Fi12Fi0.DG,Hi12Hi0.DG,inhib12c0.DG))

FC.comp.12i.exp <- FC.comp.12i[,c(1,6,10,2,16,30)]
colnames(FC.comp.12i.exp) <- c("Peak ID","Gene ID", "Swissprot Annotation",
                              "Foot Reg logFC","Head Reg logFC","Structural Enrichment")
write.csv(FC.comp.12i.exp, file = "Analysis_Output/ATAC/fc.comp.12i.csv", row.names = F)

scatterPlot(FC.comp.12i, "FC_Comp_12i")

####Identify head and foot marker genes in differentially activated peaks####

# 8hpa
FC.comp.8.str <- HRFR8c0[HRFR8c0$ID %in% FC.comp.8$ID,]

str.8 <- unique(FC.comp.8.str[complete.cases(FC.comp.8.str$G_Name),"gene_ID"])

str.8 <- strEnrichment(str.8)

str.8.exp <- str.8[[2]]

str.8 <- str.8[[1]]

FC.comp.8.str$str <- mapvalues(FC.comp.8$gene_ID_F, from = str.8$ID, to = str.8$structureRes, warn_missing = F)

FC.comp.8.str$strFC <- mapvalues(FC.comp.8$gene_ID_F, from = str.8$ID, to = str.8$strFC, warn_missing = F)

FC.comp.8.str$str[grepl("^g",FC.comp.8.str$str) | is.na(FC.comp.8.str$str)] <- "NoEnrichment"

FC.comp.8.str$strFC[is.na(FC.comp.8.str$strFC) | grepl("^g",FC.comp.8.str$strFC)] <- 0

FC.comp.8.str$strFC <- as.numeric(FC.comp.8.str$strFC)

FC.comp.8.str$str[sign(FC.comp.8.str$strFC) != sign(FC.comp.8.str$logFC) | !(FC.comp.8.str$ID %in% HRFR8c0.DG$ID)] <- "NoEnrichment"

scatterPlot(FC.comp.8.str, x = "logFC", y = "strFC", colorBy = "str", plotName = "Struct_8hpa", corPlot = F)

FC.comp.8.str.exp <- FC.comp.8.str[,c(1,6,2,17,16,10)]
colnames(FC.comp.8.str.exp) <- c("Peak ID","Gene ID","Regen LogFC", "SC LogFC", "Structure Enrichment", "Swissprot Annotation")

str.8.exp$ID <- rownames(str.8.exp)
str.8.exp <- str.8.exp[,c(6,2,3,4,5)]
colnames(str.8.exp) <- c("ID","LogFC","Percent Positive, Head Cells", "Percent Positive, Foot Cells", "Adjusted P-value")

write.csv(FC.comp.8.str.exp,"Analysis_Output/ATAC/str.comp.8.csv", row.names = F)
write.csv(str.8.exp,"Analysis_Output/ATAC/str.8.csv", row.names = F)

# 12hpa

FC.comp.12.str <- HRFR12c0[HRFR12c0$ID %in% FC.comp.12$ID,]

str.12 <- unique(FC.comp.12.str[complete.cases(FC.comp.12.str$G_Name),"gene_ID"])

str.12 <- strEnrichment(str.12)

str.12.exp <- str.12[[2]]

str.12 <- str.12[[1]]

FC.comp.12.str$str <- mapvalues(FC.comp.12$gene_ID_F, from = str.12$ID, to = str.12$structureRes, warn_missing = F)

FC.comp.12.str$strFC <- mapvalues(FC.comp.12$gene_ID_F, from = str.12$ID, to = str.12$strFC, warn_missing = F)

FC.comp.12.str$str[grepl("^g",FC.comp.12.str$str) | is.na(FC.comp.12.str$str)] <- "NoEnrichment"

FC.comp.12.str$strFC[is.na(FC.comp.12.str$strFC) | grepl("^g",FC.comp.12.str$strFC)] <- 0

FC.comp.12.str$strFC <- as.numeric(FC.comp.12.str$strFC)

FC.comp.12.str$str[sign(FC.comp.12.str$strFC) != sign(FC.comp.12.str$logFC) | !(FC.comp.12.str$ID %in% HRFR12c0.DG$ID)] <- "NoEnrichment"

scatterPlot(FC.comp.12.str, x = "logFC", y = "strFC", colorBy = "str", plotName = "Struct_12hpa", corPlot = F)

FC.comp.12.str.exp <- FC.comp.12.str[,c(1,6,2,17,16,10)]
colnames(FC.comp.12.str.exp) <- c("Peak ID","Gene ID","Regen LogFC", "SC LogFC", "Structure Enrichment", "Swissprot Annotation")

str.12.exp$ID <- rownames(str.12.exp)
str.12.exp <- str.12.exp[,c(6,2,3,4,5)]
colnames(str.12.exp) <- c("ID","LogFC","Percent Positive, Head Cells", "Percent Positive, Foot Cells", "Adjusted P-value")

write.csv(FC.comp.12.str.exp,"Analysis_Output/ATAC/str.comp.12.csv", row.names = F)
write.csv(str.12.exp,"Analysis_Output/ATAC/str.12.csv", row.names = F)

####Visualize ATAC Tracks ####

plotATAC(chr = "Sc4wPfr_224.1", start  = 65750, stop = 67500)
ggsave("plots/ATAC/Wnt910c_prom.pdf", height = 6, width = 12)

plotATAC(chr = "Sc4wPfr_1061", start  = 1489500, stop = 1490750)
ggsave("plots/ATAC/Wntless_prom.pdf", height = 6, width = 12)

plotATAC(chr = "Sc4wPfr_399", start  = 437500, stop = 440000)
ggsave("plots/ATAC/Wnt3_prom.pdf", height = 6, width = 12)

####Export Peak Lists####

write.bed <- function(x,p) {
  bedFile <- peakLink[peakLink$peak_id %in% x$ID,]
  bedFile <- bedFile[,2:5]
  bedFile$score <- 0
  bedFile$strand <- "."
  write.table(bedFile, file = p, quote = F, sep = "\t", row.names = F, col.names = F)
}

#get peaks that increase in accessibility from 0 to 3 hpa
Up.3 <- H3H0[H3H0$ID %in% c(H3H0.DG.up$ID,F3F0.DG.up$ID),]

#get control peaks that don't go up at 3
Up.3.c <- H3H0[!(H3H0$ID %in% Up.3$ID),]

write.bed(Up.3, "resources/Up3.bed")
write.bed(Up.3.c, "resources/Up3c.bed")

system("./homerEnrichment.sh resources/Up3.bed resources/Up3c.bed Up3Enrichment")


####Visualize ATAC Tracks using Gviz####
library(Gviz)
genomeInfo <- read.table("resources/Dovetail.genome", stringsAsFactors = F)

chrominfo <- data.frame(chrom=genomeInfo$V1, length=genomeInfo$V2, is_circular=rep(FALSE,length(genomeInfo$V1)))
txDb.mi <- makeTxDbFromGFF(file="resources/hydra.augustus.gff", format="gff", dataSource="Gene Models", organism ="Hydra vulgaris", chrominfo=chrominfo)

options(ucscChromosomeNames=F)

AtacPlot <- function(chr, height, left, right, buffer,fpath = "plot.pdf") {
  chrID <- chr
  
  #This specifies the maximum value on the y-axis
  z <- height

  #chromosome map
  gtrack <- GenomeAxisTrack(name = chrID,add53=T, add35 = T, fontsize = 13, fontcolor.title = "black", 
                            fontsize.title = 13, showTitle = F, rotation.title = 0, grid = T,
                            cex = 0.6, labelPos = "below")
  #gene models
  grtrack <- GeneRegionTrack(txDb.mi, chromosome=chrID, name="Genes", transcriptAnnotation = "gene", 
                             col = "black", fontsize.group = 13, fontcolor.group = "black", fill = "black", 
                             fontsize=25, rotation.title = 0, background.title = "white", col.line = "black",
                             just.group = "below", collapseTranscripts = "longest")
  
  H0_1 <- DataTrack(range = "resources/bigwigs/0H1_final_shift.bw", 
                  type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#7f0000","#7f0000"), 
                  chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                  background.title = "white", fontcolor.title = "black", col.axis = "black", 
                  span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  H0_2 <- DataTrack(range = "resources/bigwigs/0H2_final_shift.bw", 
                  type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#7f0000","#7f0000"), 
                  chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                  background.title = "white", fontcolor.title = "black", col.axis = "black", 
                  span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  H0_3 <- DataTrack(range = "resources/bigwigs/0H3_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#7f0000","#7f0000"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  H0_4 <- DataTrack(range = "resources/bigwigs/0H4_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#7f0000","#7f0000"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  H0_5 <- DataTrack(range = "resources/bigwigs/0H5_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#7f0000","#7f0000"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  F0_1 <- DataTrack(range = "resources/bigwigs/0F1_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#00007f","#00007f"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  F0_2 <- DataTrack(range = "resources/bigwigs/0F2_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#00007f","#00007f"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  F0_3 <- DataTrack(range = "resources/bigwigs/0F3_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#00007f","#00007f"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  F0_4 <- DataTrack(range = "resources/bigwigs/0F4_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#00007f","#00007f"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  F0_5 <- DataTrack(range = "resources/bigwigs/0F5_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#00007f","#00007f"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  H3_1 <- DataTrack(range = "resources/bigwigs/3H1_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#ff0000","#ff0000"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  H3_2 <- DataTrack(range = "resources/bigwigs/3H2_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#ff0000","#ff0000"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  H3_3 <- DataTrack(range = "resources/bigwigs/3H3_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#ff0000","#ff0000"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)

  F3_1 <- DataTrack(range = "resources/bigwigs/3F1_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#0000ff","#0000ff"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  F3_2 <- DataTrack(range = "resources/bigwigs/3F2_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#0000ff","#0000ff"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  F3_3 <- DataTrack(range = "resources/bigwigs/3F3_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#0000ff","#0000ff"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  H8_1 <- DataTrack(range = "resources/bigwigs/8H1_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#ff4c4c","#ff4c4c"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  H8_2 <- DataTrack(range = "resources/bigwigs/8H2_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#ff4c4c","#ff4c4c"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  H8_3 <- DataTrack(range = "resources/bigwigs/8H3_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#ff4c4c","#ff4c4c"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  H8_4 <- DataTrack(range = "resources/bigwigs/8H4_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#ff4c4c","#ff4c4c"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  H8_5 <- DataTrack(range = "resources/bigwigs/8H5_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#ff4c4c","#ff4c4c"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  F8_1 <- DataTrack(range = "resources/bigwigs/8F1_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#4c4cff","#4c4cff"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  F8_2 <- DataTrack(range = "resources/bigwigs/8F2_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#4c4cff","#4c4cff"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  F8_3 <- DataTrack(range = "resources/bigwigs/8F3_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#4c4cff","#4c4cff"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  F8_4 <- DataTrack(range = "resources/bigwigs/8F4_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#4c4cff","#4c4cff"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  F8_5 <- DataTrack(range = "resources/bigwigs/8F5_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#4c4cff","#4c4cff"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  H12_1 <- DataTrack(range = "resources/bigwigs/12H1_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#ff9999","#ff9999"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  H12_2 <- DataTrack(range = "resources/bigwigs/12H2_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#ff9999","#ff9999"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  H12_3 <- DataTrack(range = "resources/bigwigs/12H3_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#ff9999","#ff9999"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  H12_4 <- DataTrack(range = "resources/bigwigs/12H4_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#ff9999","#ff9999"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  H12_5 <- DataTrack(range = "resources/bigwigs/12H5_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#ff9999","#ff9999"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  F12_1 <- DataTrack(range = "resources/bigwigs/12F1_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#9999ff","#9999ff"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  F12_2 <- DataTrack(range = "resources/bigwigs/12F2_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#9999ff","#9999ff"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  F12_3 <- DataTrack(range = "resources/bigwigs/12F3_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#9999ff","#9999ff"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  F12_4 <- DataTrack(range = "resources/bigwigs/12F4_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#9999ff","#9999ff"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  F12_5 <- DataTrack(range = "resources/bigwigs/12F5_final_shift.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c("#9999ff","#9999ff"), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  pdf(file = fpath, width=15, height=20)
  plotTracks(list(F0_1,F0_2,F0_3,F0_4,F0_5,
                  F3_1,F3_2,F3_3,
                  F8_1,F8_2,F8_3,F8_4,F8_5,
                  F12_1,F12_2,F12_3,F12_4,F12_5,
                  H0_1,H0_2,H0_3,H0_4,H0_5,
                  H3_1,H3_2,H3_3,
                  H8_1,H8_2,H8_3,H8_4,H8_5,
                  H12_1,H12_2,H12_3,H12_4,H12_5,
                  grtrack,gtrack), 
             from=(left - buffer), to=(right + buffer), title.width = 0.7, margin = 0, sizes = c(rep(4,36),4,3))
  dev.off()
}

AtacPlot("Sc4wPfr_399",12,350000,550000,0, "plots/ATAC/wnt3_broad.pdf")

AtacPlot("Sc4wPfr_224.1",10,10,200000,0, "plots/ATAC/wnt910c_broad.pdf")

AtacPlot("Sc4wPfr_390",7,1380000,1575000,0, "plots/ATAC/pitx_broad.pdf")



