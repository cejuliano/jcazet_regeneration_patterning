args <- commandArgs(trailingOnly=TRUE)

prefix <- args[1]

print(prefix)

#function to calculate TSS enrichment score
TSScalc <- function(x) {
  #note: if you want to exclude rows (genes) that have no reads whatsoever, you can execute this line
  countMatrix <- read.table(x, skip = 3)
  y <- countMatrix
  countMatrix <- countMatrix[rowSums(countMatrix) > 0,]
  
  #calculate the average read density per bin for all genes
  readAve <- colMeans(countMatrix, na.rm = T)

  
  #calculate the averages for the 100bp in the flanking regions
  #this establishes the background levels
  end5 <- sum(readAve[1:10])/10
  end3 <- sum(readAve[190:200])/10
  normFact <- (end5+end3)/2 
  
  #normalize values to the background levels
  heatmapNorm <- readAve/normFact
  print(round(max(heatmapNorm), digits = 2))
  readAve <- readAve/(sum(readAve))
  
  #generate plot
  plot(1:length(readAve),readAve, type="n", ylab = "Average Read Density", xlab = "Position", xaxt = "n")
  lines(1:length(readAve),readAve)
  text(x = 150, y = 0.014, labels = round(max(heatmapNorm), digits = 2), cex = 1.5)
  axis(side = 1, at = c(0,100,200), labels = c("-1kb", "TSS", "+1kb"))

  #report the enrichment at TSS relative to background
  return(max(heatmapNorm))
}

#get matrix files
file.list <- list.files()

file.list <- file.list[grepl(paste0("^Values_",prefix),file.list)]

#get scores
filename <- paste0(prefix,"_TSS.pdf")

pdf(file = filename, width = (3*4), height = (ceiling(length(file.list)/3)) * 4)
par(mar=c(2,2,2,2))
par(mfrow=c(ceiling(length(file.list)/3),3))

for (i in file.list) {
	TSScalc(i)
	plotTitle <- i
	plotTitle <- gsub("Values_","",plotTitle)
	plotTitle <- gsub("_final_ATAC.txt","",plotTitle)
	title(main = plotTitle)
}

dev.off()
