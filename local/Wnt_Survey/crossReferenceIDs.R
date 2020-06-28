library(plyr)
options(stringsAsFactors = F)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

system("./crossRefBlast.sh")

#cross reference mapping from transcriptome to genome
dv2lr <- read.table("dvToLr.txt")
lr2dv <- read.table("lrToDv.txt")

#add tiny psuedovalue to evalue for log transformation
dv2lr$V11 <- dv2lr$V11 + 1e-181
lr2dv$V11 <- lr2dv$V11 + 1e-181

dv2lr$V1 <- gsub("[.].*","",dv2lr$V1)
lr2dv$V2 <- gsub("[.].*","",lr2dv$V2)

dv2lr$hit <- paste0(dv2lr$V1,"_",dv2lr$V2)
lr2dv$hit <- paste0(lr2dv$V2, "_", lr2dv$V1)

dv2lr <- dv2lr[!duplicated(dv2lr$hit),]
lr2dv <- lr2dv[!duplicated(lr2dv$hit),]

#only keep hits that are the best hit or whose evalue is
# ~5 orders of magnitude from that hit
dv2lr <- split(dv2lr, dv2lr$V1)

dv2lr <- lapply(dv2lr, function(x) {
  x[(log(min(x$V11)) - log(x$V11)) > -10,]
})

dv2lr <- do.call(rbind, dv2lr)


lr2dv <- split(lr2dv, lr2dv$V1)

lr2dv <- lapply(lr2dv, function(x) {
  x[(log(min(x$V11)) - log(x$V11)) > -10,]
})

lr2dv <- do.call(rbind, lr2dv)


dv2lr <- dv2lr[dv2lr$hit %in% lr2dv$hit,]
lr2dv <- lr2dv[lr2dv$hit %in% dv2lr$hit,]

#remove psuedovalue
dv2lr$V11 <- dv2lr$V11 - 1e-181

dv2lr$hit <- NULL

dv2lr <- dv2lr[,c(1:4,11)]

#load in LR and Dv SP annotations
lr.sp <- read.csv("LRv2_SP.csv")

dv.sp <- read.csv("../resources/expanded_dovetail_SP_annot.csv")

dv.sp$ID <- gsub("[.]","",dv.sp$ID)

dv.sp <- dv.sp[!duplicated(dv.sp$ID),]

dv2lr$dvSP <- mapvalues(dv2lr$V1, from  = dv.sp$ID, to = dv.sp$SP, warn_missing = F)

dv2lr$dvSP <- gsub(".*[|]","",dv2lr$dvSP)

dv2lr$dvSP[grepl("^g",dv2lr$dvSP)] <- NA

dv2lr$lrSP <- mapvalues(dv2lr$V2, from  = lr.sp$ID, to = lr.sp$sp_ID, warn_missing = F)

dv2lr$lrSP[grepl("^t",dv2lr$lrSP)] <- NA

colnames(dv2lr) <- c("dvID","lrID","percentIdentity","alignmentLength","eVal","dvSP","lrSP")

write.csv(dv2lr, file = "crossReferenceIDs.csv", row.names = F)

