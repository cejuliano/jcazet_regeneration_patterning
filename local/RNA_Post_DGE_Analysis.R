####Setup####

library(rstudioapi)
library(ggplot2)
library(Seurat)
library(plyr)
library(gsubfn)
library(ggfortify)
library(gt)

#set working directory to folder containing this script
setwd(dirname(getActiveDocumentContext()$path))

#load in the edgeR DGE results
load("Analysis_Output/RNA/RNA_DGE.RData")
load("Analysis_Output/RNA/icrt_RNA_DGE.RData")

#load in scRNA-seq data
ds <- readRDS("resources/Hydra_Seurat_Whole_Genome_updated.rds")

#defining functions
plotGeneExpression <- function(x, drug.include=F) {
  
  inormalizedCounts$ID <- rownames(inormalizedCounts)
  normalizedCounts$ID <- rownames(normalizedCounts)
  normalizedCounts <- merge(normalizedCounts,inormalizedCounts, by = "ID")
  
  rownames(normalizedCounts) <- normalizedCounts$ID
  
  #extract expression data for the module peaks
  moduleExpression <- normalizedCounts[normalizedCounts$ID %in% x,]
  
  #moduleExpression <- moduleExpression[,!grepl("i",colnames(moduleExpression))]
  
  #drop ID column
  moduleExpression$ID <- NULL

  #transpose expression
  moduleExpression <- as.data.frame(t(moduleExpression))
  
  moduleExpression.full <- moduleExpression
  
  #pool reps by treatment
  
  moduleExpression$Type <- as.factor(gsub("\\d","",rownames(moduleExpression)))
  
  moduleExpression$Time <- as.factor(as.numeric(gsub("\\D*R\\d$","",rownames(moduleExpression))))
  
  moduleExpression <- aggregate(moduleExpression[,1], by = list(moduleExpression$Time, moduleExpression$Type), mean)
  
  #also provide IDs for full data
  moduleExpression.full$Type <- as.factor(gsub("\\d","",rownames(moduleExpression.full)))
  
  moduleExpression.full$Time <- as.factor(as.numeric(gsub("\\D*R\\d$","",rownames(moduleExpression.full))))
  
  moduleExpression.full$drug <- grepl("i",moduleExpression.full$Type)
  
  colnames(moduleExpression) <- c("Time", "Type", "Expression")
  
  moduleExpression$drug <- grepl("i",moduleExpression$Type)
  
  moduleExpression[moduleExpression$drug == T,"drug"] <- "5µM iCRT14"
  moduleExpression[moduleExpression$drug == F,"drug"] <- "Untreated"
  
  moduleExpression.full[moduleExpression.full$drug == T,"drug"] <- "5µM iCRT14"
  moduleExpression.full[moduleExpression.full$drug == F,"drug"] <- "Untreated"
  
  moduleExpression.full$drug <- factor(moduleExpression.full$drug, levels = c("Untreated", "5µM iCRT14"))

  moduleExpression$drug <- factor(moduleExpression$drug, levels = c("Untreated", "5µM iCRT14"))
  
  if(drug.include) {
    moduleExpression$Type <- gsub("i","", moduleExpression$Type)
    moduleExpression.full$Type <- gsub("i","", moduleExpression.full$Type)
    
    moduleExpression$drug <- paste(moduleExpression$Type, moduleExpression$drug, sep = "_")
    moduleExpression.full$drug <- paste(moduleExpression.full$Type, moduleExpression.full$drug, sep = "_")
    
    gg <- ggplot(data=moduleExpression, aes(x = Time, y = Expression, group = drug, color = drug)) + geom_line(size = 1.5)
    gg <- gg + facet_grid(. ~ Type)
    gg <- gg + geom_point(data = moduleExpression.full, aes(x = Time, y = eval(parse(text = x)), group = drug, color = drug), size = 2, alpha = 0.75)
    gg <- gg + scale_color_manual(values=c("#40E0D0", "blue", "#d040e0", "red"))
  } else {
    moduleExpression <- moduleExpression[moduleExpression$drug != "5µM iCRT14",]
    moduleExpression.full <- moduleExpression.full[moduleExpression.full$drug != "5µM iCRT14",]
    gg <- ggplot(data=moduleExpression, aes(x = Time, y = Expression, group = Type, color = Type)) + geom_line(size = 1.5)
    gg <- gg + geom_point(data = moduleExpression.full, aes(x = Time, y = eval(parse(text = x)), group = Type, color = Type), size = 2, alpha = 0.75)
    gg <- gg + scale_color_manual(values = c("blue","red"))
    }
  
  gg <- gg + theme_bw()
  
  
  gg <- gg + ggtitle(paste0(x, " ", gsub("sp[|].*[|]","",Annotations[which(Annotations$ID == x),3])))
  
  gg <- gg + labs(y = "log2(CPM)")
  
  gg <- gg + theme(legend.position = "none")
  
  gg
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
  gg <- gg + geom_point(data = df, aes(x = !!sym(x), y = !!sym(y), color = !!sym(colorBy), 
                                       size = !!sym(colorBy), name = ID))
  gg <- gg + scale_color_manual(values = colors.use)
  gg <- gg + scale_size_manual(values = sizes.use)
  
  if(corPlot) {
    gg <- gg + geom_abline(intercept = 0, slope = 1, linetype = "longdash", size = 0.8)
    gg <- gg + scale_y_continuous(limits = c(axis.llimit,axis.ulimit))
    gg <- gg + scale_x_continuous(limits = c(axis.llimit,axis.ulimit))
  }
  
  gg
  
  ggsave(filename = paste0("plots/RNA/",plotName,".pdf"), plot = gg, width = 4, height = 4, useDingbats=FALSE)
  
  return(gg)
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

hFind <-function (y,x) { 
  return (y@assays$RNA@data@Dimnames[[1]][grep(x,y@assays$RNA@data@Dimnames[[1]],ignore.case = T)])
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
  
  structureRes <- data.frame(structureRes = structureRes, 
                             strFC = mapvalues(x, from = enrichmentRes$ID, to = enrichmentRes$avg_logFC, warn_missing = F))
  
  structureRes$strFC <- as.numeric(gsub("^g.*",NA, structureRes$strFC))
  
  return(structureRes)
}

####Visualize differentially activated genes####

FC.comp.3 <- fcCompare(F3vF0,H3vH0,list(F3vF0.DG,H3vH0.DG,HR3vFR3.DG))

scatterPlot(FC.comp.3, "FC_Comp_3")

FC.comp.8 <- fcCompare(F8vF0,H8vH0,list(F8vF0.DG,H8vH0.DG,HR8vFR8c0.DG))

scatterPlot(FC.comp.8, "FC_Comp_8")

FC.comp.8i <- fcCompare(Fi8vFi0,Hi8vHi0,list(Fi8vFi0.DG,Hi8vHi0.DG,Hi8vFi8c0.DG))

scatterPlot(FC.comp.8i, "FC_Comp_8i")

FC.comp.12 <- fcCompare(F12vF0,H12vH0,list(F12vF0.DG,H12vH0.DG,HR12vFR12c0.DG))

scatterPlot(FC.comp.12,"FC_Comp_12")

FC.comp.12i <- fcCompare(Fi12vFi0,Hi12vHi0,list(Fi12vFi0,Hi12vHi0,Hi12vFi12c0.DG))

scatterPlot(FC.comp.12i,"FC_Comp_12i")


####Identify head and foot marker genes in differentially activated genes####
FC.comp.8.str <- HR8vFR8c0[HR8vFR8c0$ID %in% FC.comp.8$ID,]

str.8 <- strEnrichment(FC.comp.8.str$ID)

FC.comp.8.str$str <- str.8[,1]
FC.comp.8.str$strFC <- str.8[,2]

FC.comp.8.str <- FC.comp.8.str[!is.na(FC.comp.8.str$strFC),]

FC.comp.8.str$str[!(FC.comp.8.str$ID %in% HR8vFR8c0.DG$ID)] <- "NoEnrichment"

FC.comp.8.str$str[sign(FC.comp.8.str$strFC) != sign(FC.comp.8.str$logFC)] <- "NoEnrichment"

scatterPlot(FC.comp.8.str, x = "logFC", y = "strFC", colorBy = "str", plotName = "Struct_8hpa", corPlot = F)

#repeat for 12hpa
FC.comp.12.str <- HR12vFR12c0[HR12vFR12c0$ID %in% FC.comp.12$ID,]

str.12 <- strEnrichment(FC.comp.12.str$ID)

FC.comp.12.str$str <- str.12[,1]
FC.comp.12.str$strFC <- str.12[,2]

FC.comp.12.str <- FC.comp.12.str[!is.na(FC.comp.12.str$strFC),]

FC.comp.12.str$str[!(FC.comp.12.str$ID %in% HR12vFR12c0.DG$ID)] <- "NoEnrichment"

FC.comp.12.str$str[sign(FC.comp.12.str$strFC) != sign(FC.comp.12.str$logFC)] <- "NoEnrichment"

scatterPlot(FC.comp.12.str, x = "logFC", y = "strFC", colorBy = "str", plotName = "Struct_12hpa", corPlot = F)
                
####plot noteworthy genes####
unlink("plots/RNA/Noteworthy_Genes", recursive = T)
dir.create("plots/RNA/Noteworthy_Genes", showWarnings = F)

#wnt3
plotGeneExpression("g28064")
ggsave("plots/RNA/Noteworthy_Genes/Wnt3.pdf", height = 6, width = 6)

plotGeneExpression("g28064", drug.include = T)
ggsave("plots/RNA/Noteworthy_Genes/Wnt3_i.pdf", height = 6, width = 11)

#wnt910c
plotGeneExpression("g33373")
ggsave("plots/RNA/Noteworthy_Genes/Wnt910c.pdf", height = 6, width = 6)

plotGeneExpression("g33373", drug.include = T)
ggsave("plots/RNA/Noteworthy_Genes/Wnt910c_i.pdf", height = 6, width = 11)

#dvl
plotGeneExpression("g29781")
ggsave("plots/RNA/Noteworthy_Genes/dvl.pdf", height = 6, width = 6)

#wntless
plotGeneExpression("g18842")
ggsave("plots/RNA/Noteworthy_Genes/Wntless.pdf", height = 6, width = 6)

plotGeneExpression("g18842", drug.include = T)
ggsave("plots/RNA/Noteworthy_Genes/Wntless_i.pdf", height = 6, width = 11)

#bcat
plotGeneExpression("g7262")
ggsave("plots/RNA/Noteworthy_Genes/bCat.pdf", height = 6, width = 6)

#wnt7
plotGeneExpression("g29750")
ggsave("plots/RNA/Noteworthy_Genes/wnt7.pdf", height = 6, width = 6)

#pitx
plotGeneExpression("g5621")
ggsave("plots/RNA/Noteworthy_Genes/pitx.pdf", height = 6, width = 6)

#otx
plotGeneExpression("g33396")
ggsave("plots/RNA/Noteworthy_Genes/otx.pdf", height = 6, width = 6)

#jun
plotGeneExpression("g1450")
ggsave("plots/RNA/Noteworthy_Genes/jun.pdf", height = 6, width = 6)

plotGeneExpression("g1450", drug.include = T)
ggsave("plots/RNA/Noteworthy_Genes/jun_i.pdf", height = 6, width = 11)

#tcf
plotGeneExpression("g27364")
ggsave("plots/RNA/Noteworthy_Genes/tcf.pdf", height = 6, width = 6)

#sp5
plotGeneExpression("g33422")
ggsave("plots/RNA/Noteworthy_Genes/sp5.pdf", height = 6, width = 6)

#sp5
plotGeneExpression("g33422", drug.include = T)
ggsave("plots/RNA/Noteworthy_Genes/sp5_i.pdf", height = 6, width = 11)

#fos
plotGeneExpression("g23720", drug.include = T)
ggsave("plots/RNA/Noteworthy_Genes/fos_i.pdf", height = 6, width = 11)

plotGeneExpression("g23720")
ggsave("plots/RNA/Noteworthy_Genes/fos.pdf", height = 6, width = 6)

#creb
plotGeneExpression("g16491", drug.include = T)
ggsave("plots/RNA/Noteworthy_Genes/creb_i.pdf", height = 6, width = 11)

plotGeneExpression("g16491")
ggsave("plots/RNA/Noteworthy_Genes/creb.pdf", height = 6, width = 6)

#cr3l
plotGeneExpression("g32585")
ggsave("plots/RNA/Noteworthy_Genes/cr3l.pdf", height = 6, width = 6)

#notum single cell plot
FeaturePlot(ds, hFind(ds, "g26902.t"), order = T, pt.size = 0.8) + NoLegend() + NoAxes()
ggsave("plots/RNA/Noteworthy_Genes/notum_sc.png", device = "png", height = 6, width = 6, dpi = 300)

####Early head and foot regeneration genes at 3hpa####
GvS <- merge(HR12vFR12c0.DG, H3vH0, by = "ID", suffixes = c("_12","_H3"))

GvS <- merge(GvS, F3vF0, by = "ID")

GvS$aveFC <- apply(GvS[,c("logFC_H3","logFC")], 1, function(x) mean(x))

#GvS <-GvS[GvS$ID %in% FC.comp.12.str[FC.comp.12.str$str != "NoEnrichment","ID"],]

GvS$Interesting <- "noEnrichment"
GvS$Interesting[(GvS$ID %in% HR12vFR12c0.DG.up$ID) & (GvS$ID %in% c(H3vH0.DG$ID,F3vF0.DG$ID))] <- "Head"
GvS$Interesting[(GvS$ID %in% HR12vFR12c0.DG.down$ID) & (GvS$ID %in% c(H3vH0.DG$ID,F3vF0.DG$ID))] <- "Foot"

scatterPlot(GvS,plotName = "EarlyWoundStructure", x = "logFC_12", y = "aveFC", colorBy = "Interesting", corPlot = F)

####icrt effect on structure specific program####

#divergently expressed genes, head regeneration at 8hpa
H8vHi8c0.headChange <- H8vHi8c0[H8vHi8c0$ID %in% HR8vFR8c0.DG$ID,]

H8vHi8c0.headChange$Interesting <- "noChange"
H8vHi8c0.headChange$Interesting[H8vHi8c0.headChange$ID %in% H8vHi8c0.DG.down$ID] <- "Less"
H8vHi8c0.headChange$Interesting[H8vHi8c0.headChange$ID %in% H8vHi8c0.DG.up$ID] <- "More"
H8vHi8c0.headChange$sigScore <- -log10(H8vHi8c0.headChange$FDR)

gg <- scatterPlot(H8vHi8c0.headChange, "drug_effect_head_8hpa", x = "logFC", y = "sigScore", colorBy = "Interesting", corPlot = F)
gg <- gg + geom_abline(intercept = 3, slope = 0, size = 0.8, linetype = "longdash")
gg <- gg + scale_color_manual(values = c("#008000","#e59400","darkgrey"))
ggsave(filename = "plots/RNA/drug_effect_head_8hpa.pdf", plot = gg, width = 4, height = 4, useDingbats=FALSE)

#divergently expressed genes, foot regeneration at 8hpa
F8vFi8c0.footChange <- F8vFi8c0[F8vFi8c0$ID %in% HR8vFR8c0.DG$ID,]

F8vFi8c0.footChange$Interesting <- "Not Interesting"
F8vFi8c0.footChange$Interesting[F8vFi8c0.footChange$ID %in% F8vFi8c0.DG.down$ID] <- "Less"
F8vFi8c0.footChange$Interesting[F8vFi8c0.footChange$ID %in% F8vFi8c0.DG.up$ID] <- "More"
F8vFi8c0.footChange$sigScore <- -log10(F8vFi8c0.footChange$FDR)

gg <- scatterPlot(F8vFi8c0.footChange, "drug_effect_foot_8hpa", x = "logFC", y = "sigScore", colorBy = "Interesting", corPlot = F)
gg <- gg + geom_abline(intercept = 3, slope = 0, size = 0.8, linetype = "longdash")
gg <- gg + scale_color_manual(values = c("#008000","#e59400","darkgrey"))
ggsave(filename = "plots/RNA/drug_effect_foot_8hpa.pdf", plot = gg, width = 4, height = 4, useDingbats=FALSE)

####check for correspondence between RNA and ATAC data####
#test for enrichment of injury induced peaks near injury induced genes

#load ATAC-seq data
attach("Analysis_Output/ATAC/ATAC_DGE.RData")
attach("Analysis_Output/RNA/wRNA_DGE.RData")


calcEnrichment <- function(sampleList,referenceList) {
  
  #define gene universe
  gU <- analysisList$genes$genes
  
  #hits in sample
  SinS <- length(intersect(sampleList,referenceList))
  
  #hits in background
  SinB <- length(referenceList) - SinS
  
  #non-hits in background
  FinB <- length(gU) - SinB
  
  #perform hypergeometric test, look for enrichment
  print('P-value')
  print(phyper(SinS,SinB,FinB,length(sampleList), lower.tail = F))
  
  #calclate fold enrichment
  print("Fold Enrichment:")
  FE <- (SinS/length(sampleList))/(SinB/(SinB + FinB))
  print(FE)
}

#test ATAC differential peak enrichment near DE transcripts
calcEnrichment(HR12vFR12c0.DG$ID, HRFR12c0.DG$gene_ID)

calcEnrichment(HR8vFR8c0.DG$ID, HRFR8c0.DG$gene_ID)

#test if our list of DE transcripts is enriched in transcripts
#that were identified as DE in the Wenger et al. dataset
calcEnrichment(HR12vFR12c0.DG$ID, HRw16vFRw16c0.DG$ID)

calcEnrichment(HR8vFR8c0.DG$ID, HRw8vFRw8c0.DG$ID)

