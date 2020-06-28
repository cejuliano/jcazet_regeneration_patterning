library(KEGGREST)
library(plyr)
library(Biostrings)
options(stringsAsFactors = F)

#set wd to folder containing this script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#get list of kegg wnt pathway genes 
#(curated from: https://www.genome.jp/dbget-bin/www_bget?pathway+hsa04310)
wntGenesList <- read.table("kegg_Wnt_Genes.txt", sep = "\t")

#pull kegg gene entries for each wnt component
results <- list()
for (i in 1:nrow(wntGenesList)) {
  results[[i]] <- keggGet(wntGenesList$V2[i])[[1]]
}
  
#pull the information on homologs across different references in the kegg database
results.genes <- lapply(results, function(x) x$GENES)

names(results.genes) <- paste0(wntGenesList$V1,".")

results.genes <- unlist(results.genes)

#function to pull wnt pathway homologs for specific species in kegg db
#and export a file contianing the sequences of those homologs for subsequent
#blast searches
getWntGenes <- function(x, results.genes) {
  
  results.genes <- results.genes[grepl(x, results.genes, ignore.case = T)]
  
  results.genes <- gsub(paste0(x,": "),"",results.genes, ignore.case = T)
  
  results.genes <- strsplit(results.genes, split = " ")
  
  results.genes <- unlist(results.genes)
  
  results.genes <- gsub("\\(.*$","",results.genes)
  
  tmp.names <- names(results.genes)
  
  results.genes <- paste0(x,":",results.genes)
  
  results.genes.save <- results.genes
  
  names(results.genes) <- paste0(tmp.names,"_",results.genes.save)
  
  wntGenes <- split(results.genes, ceiling(seq_along(results.genes)/10))
  
  results <- list()
  for (i in seq_along(wntGenes)) {
    workingReturn <- keggGet(wntGenes[[i]])
    names(workingReturn) <- names(wntGenes[[i]])
    results <- append(results, workingReturn)
  }
  
  results.seq <- lapply(results, function(x) x$AASEQ)
  
  for (i in seq_along(results.seq)) {
    names(results.seq[[i]]) <- names(results.seq)[i]
  }

  unlink(paste0(x,"WntGenes.fasta"))
  lapply(results.seq, function(y) 
    writeXStringSet(y,filepath = paste0(x,"WntGenes.fasta"), append = T))
  
  return(results.genes.save)
}

#get wnt genes for Hydra (old genome reference from ncbi)
getWntGenes("hmg", results.genes)

#Because the original Hydra genome is very fragmented, the kegg db
#may be missing a number of genes of interest
#to try and find additional candidates we'll use sequences from
#humans as well as other cnidarians to find matches in the 2.0 genome

#homo sapiens wnt genes
humanNames <- getWntGenes("hsa", results.genes)

#nematostella vectensis wnt genes 
nveNames <- getWntGenes("nve", results.genes)

#exaiptasia pallida wnt genes
epaNames <- getWntGenes("epa", results.genes)

#perform blast search to look for the wnt pathway homologs in hydra
system("./initial_blast.sh")

#we're using uniprot entries for our human reference, so we need to
#convert the kegg IDs to uniprot IDs
humanIDs <- data.frame(kg = character(0), matchName = character(0))
for (i in 1:length(humanNames)) {
  humanIDs[nrow(humanIDs)+1,] <- c(humanNames[i],keggConv("uniprot", humanNames[i])[1])
}

humanIDs$matchName <- gsub("up:","",humanIDs$matchName)

#our nematostella reference is from refseq, so we need to convert
#the kegg IDs to ncbi IDs
nveIDs <- data.frame(kg = character(0), matchName = character(0))
for (i in 1:length(nveNames)) {
  nveIDs[nrow(nveIDs)+1,] <- c(nveNames[i], keggConv("ncbi-proteinid", nveNames[i]))
}

nveIDs$matchName <- gsub("ncbi-proteinid:","",nveIDs$matchName)

#the exaiptasia reference is also from refseq
epaIDs <- data.frame(kg = character(0), matchName = character(0))
for (i in 1:length(epaNames)) {
  epaIDs[nrow(epaIDs)+1,] <- c(epaNames[i], keggConv("ncbi-proteinid", epaNames[i]))
}

epaIDs$matchName <- gsub("ncbi-proteinid:","",epaIDs$matchName)

#import the blast results containing the candidate wnt genes in Hydra
hm2hv <- read.table("hmToHv.txt")
nv2hv <- read.table("nvToHv.txt")
hs2hv <- read.table("hsToHv.txt")
ep2hv <- read.table("epToHv.txt")

#compile the candidate genes
wntCandidates <- unique(c(hm2hv$V2, nv2hv$V2, hs2hv$V2, ep2hv$V2))

#we next export these candidates and pull their sequence from the
#2.0 genome reference to perform a reciprocal blast search
write.table(wntCandidates, file = "hvCandidates.txt",
            sep = "\t", row.names = F, col.names = F, quote = F)

system("source ../resources/venv/bin/activate;./fetchHvSeq.sh")
system("./recipBlast.sh")

#now we import the results from the reciprocal search
#first for nematostella
hv2nv <- read.table("hvToNve.txt")

#reformat IDs (drop isoform label)
hv2nv$V2 <- gsub("[.].*","",hv2nv$V2)

#look for hits recovered in both the initial blast search and the
#reciprocal blast search
hv2nv <- hv2nv[hv2nv$V2 %in% nveIDs$matchName,]

#get the kegg ID for the nematostella hit
hv2nv$kg <- mapvalues(hv2nv$V2, from = nveIDs$matchName, to = nveIDs$kg, warn_missing = F)

nv2hv$kg <- gsub("^.*_","",nv2hv$V1)

#keep only hits that gave the same matches for the reciprocal blast searches
nv2hv <- nv2hv[paste0(nv2hv$kg,"_",nv2hv$V2) %in% 
                 paste0(hv2nv$kg,"_",hv2nv$V1),]

nv2hv <- nv2hv[!duplicated(nv2hv$V2),]

#repeat process for exaiptasia hits
hv2ep <- read.table("hvToEpa.txt")

hv2ep$V2 <- gsub("[.].*","",hv2ep$V2)

hv2ep <- hv2ep[hv2ep$V2 %in% epaIDs$matchName,]

hv2ep$kg <- mapvalues(hv2ep$V2, from = epaIDs$matchName, to = epaIDs$kg, warn_missing = F)

ep2hv$kg <- gsub("^.*_","",ep2hv$V1)

ep2hv <- ep2hv[paste0(ep2hv$kg,"_",ep2hv$V2) %in% 
                 paste0(hv2ep$kg,"_",hv2ep$V1),]

ep2hv <- ep2hv[!duplicated(ep2hv$V2),]

#repeat process for human hits
hv2hs <- read.table("hvTohs.txt")

hv2hs$V2 <- gsub("sp[|]","",hv2hs$V2)
hv2hs$V2 <- gsub("[|].*","",hv2hs$V2)

hv2hs <- hv2hs[hv2hs$V2 %in% humanIDs$matchName,]

hv2hs$kg <- mapvalues(hv2hs$V2, from = humanIDs$matchName, to = humanIDs$kg, warn_missing = F)

hs2hv$kg <- gsub("^.*_","",hs2hv$V1)

hs2hv <- hs2hv[paste0(hs2hv$kg,"_",hs2hv$V2) %in% 
                 paste0(hv2hs$kg,"_",hv2hs$V1),]

hs2hv <- hs2hv[!duplicated(hs2hv$V2),]

#combine all potential wnt genes 
pooledCandidates <- rbind(hs2hv,nv2hv, ep2hv)

pooledCandidates$kg <- NULL

#bring in the genes that were already annotated as hydra wnt genes from kegg
pooledCandidates <- rbind(pooledCandidates, hm2hv)

pooledCandidates <- pooledCandidates[!duplicated(pooledCandidates$V2),]

pooledCandidates$V1 <- gsub("[.].*","",pooledCandidates$V1)

#generate a heatmap of wnt genes
library(edgeR)
library(gplots)

#pull in regeneration RNA-seq data
load("../Analysis_Output/RNA/RNA_DGE.RData")

wntGenesRegen <- gsub("[.].*","",pooledCandidates$V2)

#manually curate list
wntGenesRegen <- wntGenesRegen[!(wntGenesRegen %in% c("g23666","g15494","g21815",
                                                    "g12374","g28262","g6790"))]

#adding sp5, sfrp, and wntless
wntGenesRegen <- c(wntGenesRegen, "g33422","g12274","g18842")

#get expression data for wnt genes
wntGenesRegen <- as.data.frame(t(normalizedCounts[rownames(normalizedCounts) %in% wntGenesRegen,]))

#because the single cell trajectory analysis used a transcriptome reference, we need to get
#a table that allows us to convert between genome and transcriptome references
lrIDs <- read.csv("crossReferenceIDs.csv")

#we will also pull in the precomputed trajectory analysis results 
load("ecto_endo_spline.RData")

#get genome ID equivalent for each transcriptome ID in the trajectory results
rownames(ms.both) <- mapvalues(rownames(ms.both), from = lrIDs$lrID, to = lrIDs$dvID, warn_missing = F)

ms.both <- ms.both[!duplicated(rownames(ms.both)),]

#drop transcripts with very low expression in the dataset
ms.both.cut <- which(apply(ms.both, 1, max) > 0.01)

ms.both <- ms.both[ms.both.cut,]

#normalize the results to have the maximum value be 1 and the minimum value be 0 
ms.both <- t(apply(ms.both, 1, function(x) x - min(x)))
ms.both <- t(apply(ms.both, 1, function(x) x/max(x)))

#subset to only include wnt genes
ms.both.plot <- ms.both[rownames(ms.both) %in% colnames(wntGenesRegen),]


#generate trajectory results heatmap
getPalette = colorRampPalette(c("white","blue"))

geneNames <- rownames(ms.both.plot)
geneNames.wnt <- mapvalues(geneNames, from = gsub("[.].*","",pooledCandidates$V2), 
                           to = pooledCandidates$V1, warn_missing = F)
geneNames <- paste0(geneNames, " ", geneNames.wnt)

pdf(file = "ds_spline_heatmap.pdf", height = 10, width = 12)
geneOrder <- heatmap.2(as.matrix(ms.both.plot), Colv = F, Rowv = T, dendrogram = "none", 
                       col = getPalette(90), trace = "none", density.info = "none",
                       key = F, cexRow = 0.35, margins = c(4,4), scale = "none", keysize = 0.2, 
                       colsep = c(16,40,46,71,117,133), labCol = NA, labRow = geneNames)
dev.off()

#save order from trajectory heatmap so the RNA-seq heatmap will have rows in the same order
geneOrder <- rownames(ms.both.plot)[geneOrder$rowInd]
geneOrder <- geneOrder[seq(length(geneOrder),1, by = -1)]

#next we generate the RNA-seq heatmap

#we'll to pool biological replicates
wntGenesRegen$treatment <- gsub("\\d+$","",rownames(wntGenesRegen))
wntGenesRegen <- aggregate(wntGenesRegen[,1:(ncol(wntGenesRegen)-1)],by = list(wntGenesRegen$treatment), mean)

rownames(wntGenesRegen) <- wntGenesRegen$Group.1

wntGenesRegen$Group.1 <- NULL

#reorder column so that they are sorted by time and regeneration type
wntGenesRegen <- wntGenesRegen[c(1,5,7,3,2,6,8,4),]

getPalette = colorRampPalette(c("white","red"))

#normalize the results to have the maximum value be 1 and the minimum value be 0 
wntGenesRegen <- apply(wntGenesRegen, 2, function(x) x - min(x))
wntGenesRegen <- as.data.frame(apply(wntGenesRegen, 2, function(x) x/max(x)))

wntGenesRegen <- wntGenesRegen[,geneOrder]

heatRowNames <- mapvalues(colnames(wntGenesRegen), from = gsub("[.].*","",pooledCandidates$V2), 
                          pooledCandidates$V1, warn_missing = F)

heatRowNames <- paste0(heatRowNames, " ", colnames(wntGenesRegen))

pdf(file = "wntGenesRegen.pdf", height = 10, width = 6)
heatmap.2(t(wntGenesRegen), Rowv = F, Colv = F, dendrogram = "none", scale = "none",
          col = getPalette(90), trace = "none", key = FALSE, sepcolor = "white",
          keysize = 0.1, margins = c(8,12), colsep = c(4), labRow = heatRowNames)
dev.off()


#check for siginificant changes in the RNA seq data for all the considered Wnt genes
differential.lists <- paste0(ls(pattern = "^[HF][1-9]+v[HF]0.DG."),"$ID")

checkSig <- function(x) {
  differentials <- lapply(differential.lists, function(y) {
    if(x %in% eval(parse(text=y))) {
      return(gsub("[$]ID","",y))
    }
  })
  differentials <- unlist(differentials)
  return(differentials)
}

wntGeneDifferentials <- lapply(colnames(wntGenesRegen), function(x) checkSig(x))
names(wntGeneDifferentials) <- colnames(wntGenesRegen)
