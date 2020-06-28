# this script will set up the object for subsequent chromVar analysis by generating fragment counts per peak
# across the untreated ATAC-seq samples

####Setup####
library(chromVAR)
library(BiocParallel)
library(BSgenome.Hymag.NIH.Hymag2)
library(plyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

register(MulticoreParam(workers = 6))

#need this to prevent bed files from being written with coordinates in scientific notation
options(scipen=999)

####Create bed files with regularly sized peaks####
peakfile <- "resources/regen_peakset.bed"
peaks <- getPeaks(peakfile, sort_peaks = T)

#load seq lengths to granges object
contigLengths <- read.table("resources/Dovetail.genome")

peaks@seqinfo@seqlengths <- as.integer(mapvalues(peaks@seqinfo@seqnames, 
                                                 from = contigLengths$V1, 
                                                 to = contigLengths$V2, 
                                                 warn_missing = F))

peaks <- resize(peaks, width = 250, fix = "center")

peaks <- trim(peaks)

####Read in read counts####
bamfiles <- c("resources/Bams/0H1_final_shift.bam",
              "resources/Bams/0H2_final_shift.bam",
              "resources/Bams/0H3_final_shift.bam",
              "resources/Bams/0H4_final_shift.bam",
              "resources/Bams/0H5_final_shift.bam",
              "resources/Bams/0F1_final_shift.bam",
              "resources/Bams/0F2_final_shift.bam",
              "resources/Bams/0F3_final_shift.bam",
              "resources/Bams/0F4_final_shift.bam",
              "resources/Bams/0F5_final_shift.bam",
              "resources/Bams/3H1_final_shift.bam",
              "resources/Bams/3H2_final_shift.bam",
              "resources/Bams/3H3_final_shift.bam",
              "resources/Bams/3F1_final_shift.bam",
              "resources/Bams/3F2_final_shift.bam",
              "resources/Bams/3F3_final_shift.bam",
              "resources/Bams/8H1_final_shift.bam",
              "resources/Bams/8H2_final_shift.bam",
              "resources/Bams/8H3_final_shift.bam",
              "resources/Bams/8H4_final_shift.bam",
              "resources/Bams/8H5_final_shift.bam",
              "resources/Bams/8F1_final_shift.bam",
              "resources/Bams/8F2_final_shift.bam",
              "resources/Bams/8F3_final_shift.bam",
              "resources/Bams/8F4_final_shift.bam",
              "resources/Bams/8F5_final_shift.bam",
              "resources/Bams/12H1_final_shift.bam",
              "resources/Bams/12H2_final_shift.bam",
              "resources/Bams/12H3_final_shift.bam",
              "resources/Bams/12H4_final_shift.bam",
              "resources/Bams/12H5_final_shift.bam",
              "resources/Bams/12F1_final_shift.bam",
              "resources/Bams/12F2_final_shift.bam",
              "resources/Bams/12F3_final_shift.bam",
              "resources/Bams/12F4_final_shift.bam",
              "resources/Bams/12F5_final_shift.bam")

#configure chromVar count object
fragment_counts <- chromVAR::getCounts(bamfiles, peaks, 
                                       paired =  TRUE, 
                                       by_rg = FALSE, 
                                       format = "bam", 
                                       colData = DataFrame(sampleType = c(rep("Head_Uninjured", 5),
                                                                          rep("Foot_Uninjured", 5),
                                                                          rep("Head_Regeneration_3hpa", 3),
                                                                          rep("Foot_Regeneration_3hpa", 3),
                                                                          rep("Head_Regeneration_8hpa", 5),
                                                                          rep("Foot_Regeneration_8hpa", 5),
                                                                          rep("Head_Regeneration_12hpa", 5),
                                                                          rep("Foot_Regeneration_12hpa", 5)
                                                                          ),
                                                           timePoint = c(rep("0h",10),
                                                                         rep("3h",6),
                                                                         rep("8h",10),
                                                                         rep("12h",10)
                                                                         ),
                                                           
                                                           structure = c(rep("head",5),
                                                                         rep("foot",5),
                                                                         rep("head",3),
                                                                         rep("foot",3),
                                                                         rep("head",5),
                                                                         rep("foot",5),
                                                                         rep("Head", 5),
                                                                         rep("Foot", 5)
                                                                         )
                                                           )
                                       )

#filter overlapping peaks
fragment_counts <- filterPeaks(fragment_counts, non_overlapping = TRUE)

#calculate GC bias
fragment_counts <- addGCBias(fragment_counts, genome = BSgenome.Hymag.NIH.Hymag2)

#save chromvar object
save(list = "fragment_counts", file = "resources/untreated_chromVar_Counts.RData")

####repeat for all samples (including iCRT14 treatment)####

bamfiles <- c("resources/Bams/0H1_final_shift.bam",
              "resources/Bams/0H2_final_shift.bam",
              "resources/Bams/0H3_final_shift.bam",
              "resources/Bams/0H4_final_shift.bam",
              "resources/Bams/0H5_final_shift.bam",
              "resources/Bams/0F1_final_shift.bam",
              "resources/Bams/0F2_final_shift.bam",
              "resources/Bams/0F3_final_shift.bam",
              "resources/Bams/0F4_final_shift.bam",
              "resources/Bams/0F5_final_shift.bam",
              "resources/Bams/3H1_final_shift.bam",
              "resources/Bams/3H2_final_shift.bam",
              "resources/Bams/3H3_final_shift.bam",
              "resources/Bams/3F1_final_shift.bam",
              "resources/Bams/3F2_final_shift.bam",
              "resources/Bams/3F3_final_shift.bam",
              "resources/Bams/8H1_final_shift.bam",
              "resources/Bams/8H2_final_shift.bam",
              "resources/Bams/8H3_final_shift.bam",
              "resources/Bams/8H4_final_shift.bam",
              "resources/Bams/8H5_final_shift.bam",
              "resources/Bams/8F1_final_shift.bam",
              "resources/Bams/8F2_final_shift.bam",
              "resources/Bams/8F3_final_shift.bam",
              "resources/Bams/8F4_final_shift.bam",
              "resources/Bams/8F5_final_shift.bam",
              "resources/Bams/12H1_final_shift.bam",
              "resources/Bams/12H2_final_shift.bam",
              "resources/Bams/12H3_final_shift.bam",
              "resources/Bams/12H4_final_shift.bam",
              "resources/Bams/12H5_final_shift.bam",
              "resources/Bams/12F1_final_shift.bam",
              "resources/Bams/12F2_final_shift.bam",
              "resources/Bams/12F3_final_shift.bam",
              "resources/Bams/12F4_final_shift.bam",
              "resources/Bams/12F5_final_shift.bam",
              "resources/Bams/0iH1_final_shift.bam",
              "resources/Bams/0iH2_final_shift.bam",
              "resources/Bams/0iH3_final_shift.bam",
              "resources/Bams/0iH4_final_shift.bam",
              "resources/Bams/0iF1_final_shift.bam",
              "resources/Bams/0iF2_final_shift.bam",
              "resources/Bams/0iF3_final_shift.bam",
              "resources/Bams/0iF4_final_shift.bam",
              "resources/Bams/3iH1_final_shift.bam",
              "resources/Bams/3iH2_final_shift.bam",
              "resources/Bams/3iH3_final_shift.bam",
              "resources/Bams/3iH4_final_shift.bam",
              "resources/Bams/3iF1_final_shift.bam",
              "resources/Bams/3iF2_final_shift.bam",
              "resources/Bams/3iF3_final_shift.bam",
              "resources/Bams/3iF4_final_shift.bam",
              "resources/Bams/3iF5_final_shift.bam",
              "resources/Bams/8iH1_final_shift.bam",
              "resources/Bams/8iH2_final_shift.bam",
              "resources/Bams/8iH3_final_shift.bam",
              "resources/Bams/8iH4_final_shift.bam",
              "resources/Bams/8iH5_final_shift.bam",
              "resources/Bams/8iF1_final_shift.bam",
              "resources/Bams/8iF2_final_shift.bam",
              "resources/Bams/8iF3_final_shift.bam",
              "resources/Bams/8iF4_final_shift.bam",
              "resources/Bams/12iH1_final_shift.bam",
              "resources/Bams/12iH2_final_shift.bam",
              "resources/Bams/12iH3_final_shift.bam",
              "resources/Bams/12iH4_final_shift.bam",
              "resources/Bams/12iF1_final_shift.bam",
              "resources/Bams/12iF2_final_shift.bam",
              "resources/Bams/12iF3_final_shift.bam",
              "resources/Bams/12iF4_final_shift.bam",
              "resources/Bams/12iF5_final_shift.bam")

#configure chromVar count object
fragment_counts <- chromVAR::getCounts(bamfiles, peaks, 
                                       paired =  TRUE, 
                                       by_rg = FALSE, 
                                       format = "bam", 
                                       colData = DataFrame(sampleType = c(rep("Head_Uninjured", 5),
                                                                          rep("Foot_Uninjured", 5),
                                                                          rep("Head_Regeneration_3hpa", 3),
                                                                          rep("Foot_Regeneration_3hpa", 3),
                                                                          rep("Head_Regeneration_8hpa", 5),
                                                                          rep("Foot_Regeneration_8hpa", 5),
                                                                          rep("Head_Regeneration_12hpa", 5),
                                                                          rep("Foot_Regeneration_12hpa", 5),
                                                                          rep("Head_Uninjured", 4),
                                                                          rep("Foot_Uninjured", 4),
                                                                          rep("Head_Regeneration_3hpa", 4),
                                                                          rep("Head_Regeneration_3hpa", 5),
                                                                          rep("Head_Regeneration_8hpa", 5),
                                                                          rep("Foot_Regeneration_8hpa", 4),
                                                                          rep("Head_Regeneration_12hpa", 4),
                                                                          rep("Foot_Regeneration_12hpa", 5)
                                       ),
                                       timePoint = c(rep("0h",10),
                                                     rep("3h",6),
                                                     rep("8h",10),
                                                     rep("12h",10),
                                                     rep("0h",8),
                                                     rep("3h",9),
                                                     rep("8h",9),
                                                     rep("12h",9)
                                       ),
                                       
                                       structure = c(rep("head",5),
                                                     rep("foot",5),
                                                     rep("head",3),
                                                     rep("foot",3),
                                                     rep("head",5),
                                                     rep("foot",5),
                                                     rep("Head", 5),
                                                     rep("Foot", 5),
                                                     rep("head",4),
                                                     rep("foot",4),
                                                     rep("Head",4),
                                                     rep("Foot",5),
                                                     rep("head",5),
                                                     rep("foot",4),
                                                     rep("Head",4),
                                                     rep("Foot",5)
                                       ),
                                       treatment = c(rep("untreated",36),
                                                     rep("treated", 35))
                                       )
)
#filter overlapping peaks
fragment_counts <- filterPeaks(fragment_counts, non_overlapping = TRUE)

#calculate GC bias
fragment_counts <- addGCBias(fragment_counts, genome = BSgenome.Hymag.NIH.Hymag2)

#save chromvar object
save(list = "fragment_counts", file = "resources/chromVar_iCounts.RData")
