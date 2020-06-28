# set working directory to folder containing this script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(stringsAsFactors = F)

#generate pairwise motif similarity matrix
system("perl ~/Homer/bin/compareMotifs.pl \\
       resources/chromVar_HOMER.motifs resources/HOMER_similarity_test \\
       -matrix resources/HOMER_similarity.txt")

#import similarity matrix for JASPAR animal PWMs
similarityMat <- read.table("resources/HOMER_similarity.txt", sep = "\t")

#for some reason there's a duplicate entry, so we need to drop that
keepThese <- which(!duplicated(similarityMat[,1]))
similarityMat <- similarityMat[keepThese,keepThese]

#reformat row and column names
labels <- similarityMat[2:nrow(similarityMat),1]
similarityMat <- similarityMat[2:nrow(similarityMat),2:ncol(similarityMat)]
rownames(similarityMat) <- labels
colnames(similarityMat) <- labels

#convert entries to numeric values (convert to matrix)
similarityMat <- apply(similarityMat,2,as.numeric)

#conversion to a matrix drops rownames for some reason
rownames(similarityMat) <- colnames(similarityMat)

#calculate distance matrix from correlation coefficient
d <- as.dist(1 - similarityMat)

#use hierarchical clustering to group motifs based on correlation
hc1 <- hclust(d, method = "average")

#plot tree (visualize proposed cutoff)
pdf(file = "plots/ATAC/HOMER_Motif_Dendrogram.pdf", height = 20, width = 80)
plot(hc1, labels = gsub("[/].*$","",colnames(similarityMat)))
abline(a = 0.2, b = 0, col = "red")
dev.off()

#cut the branches at 2 and group accordingly
clusts <- cutree(hc1, h = 0.2)
length(unique(clusts))

#for clusters with more than one motif, we will pick a representative motif by finding
#the PWM that has the highest enrichment score when comparing accessible regions to 
#non-accessible regions using HOMER (i.e. sequences within peaks against random genomic 
#background)

length(which(table(clusts) > 1))

#save cluster information
repMotifs <- data.frame(ID = names(clusts), clust = clusts)

write.csv(repMotifs, file = "resources/HOMER_motif_clusters.csv")
