suppressMessages(library(DiffBind, lib.loc = "/group/julianolab/jacazet/R/library"))
options(scipen=999)

args <- commandArgs(trailingOnly=TRUE)

prefix <- as.character(args[1])
print(prefix)

#get file list
file.list <- list.files()

file.list <- file.list[grepl(paste0("^",prefix),file.list)]

#isolate bio rep idr peakfiles
bio.idr <- file.list[grepl("IDR.narrowPeak$",file.list, ignore.case = T)]

#drop pseudoreps and pooled reps
bio.idr <- bio.idr[!grepl("_PR",bio.idr, ignore.case = T)]
bio.idr <- bio.idr[!grepl("_MG_",bio.idr, ignore.case = T)]

#get bam files
bam.files <- file.list[grepl("final_shift[.]bam$",file.list, ignore.case = T)]


bioR <- data.frame(SampleID = paste0("bioR_",1:length(bio.idr)),
                         Condition = rep("Peaks", length(bio.idr)),
                         Treatment = rep("regen", length(bio.idr)),
                         Replicate = 1:length(bio.idr),
                         bamReads = rep(bam.files[1], length(bio.idr)),
                         Peaks = bio.idr,
                         PeakCaller = rep("narrow",length(bio.idr))
)

bioPeaks <- dba(sampleSheet = bioR, minOverlap = 3)

bioPeaks <- dba.peakset(bioPeaks, bRetrieve = T)

df <- data.frame(seqnames=seqnames(bioPeaks),
                 starts=start(bioPeaks),
                 ends=end(bioPeaks),
                 names=1:length(bioPeaks),
                 scores=c(rep(".", length(bioPeaks))),
                 strands=c(rep(".", length(bioPeaks))))


write.table(df, file=paste0(prefix,"_","bioReps.consensus.bed"), quote=F, sep="\t", row.names=F, col.names=F)

#now do pseudorep idr from pooled reps
ps.idr <- list.files()
ps.idr <- ps.idr[grepl(paste0("^",prefix),ps.idr)]
ps.idr <- ps.idr[grepl(".*MG_PR.*MG_PR.*Peak$",ps.idr, ignore.case = T)]


psR <- data.frame(SampleID = paste0("psR_",1:length(ps.idr)),
                   Condition = rep("Peaks", length(ps.idr)),
                   Treatment = rep("regen", length(ps.idr)),
                   Replicate = 1:length(ps.idr),
                   bamReads = rep(bam.files[1], length(ps.idr)),
                   Peaks = ps.idr,
                   PeakCaller = rep("narrow",length(ps.idr))
)

psPeaks <- dba(sampleSheet = psR, minOverlap = 3)

psPeaks <- dba.peakset(psPeaks, bRetrieve = T)

df <- data.frame(seqnames=seqnames(psPeaks),
                 starts=start(psPeaks),
                 ends=end(psPeaks),
                 names=1:length(psPeaks),
                 scores=c(rep(".", length(psPeaks))),
                 strands=c(rep(".", length(psPeaks))))


write.table(df, file=paste0(prefix,"_","psReps.consensus.bed"), quote=F, sep="\t", row.names=F, col.names=F)
