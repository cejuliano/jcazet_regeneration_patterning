#This script uses DiffBind to count the reads for each peak in the consensus peakset and export it as RData

library(DiffBind)

#Set working directory to be the location of this R script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

writeBed <- function(x,y) {
  z <- data.frame(chrom = x[,1], start = x[,2], end = x[,3], name = rownames(x), score = x[,4], strand = ".")
  write.table(z, file = y, quote = F, sep = "\t", row.names = F, col.names = F)
}

#build sampleSheet
ATACSample <- data.frame(SampleID = c("0hpa_H1", "0hpa_H2", "0hpa_H3", "0hpa_H4", "0hpa_H5",
                                      "0hpa_F1", "0hpa_F2", "0hpa_F3", "0hpa_F4", "0hpa_F5",
                                      "3hpa_H1", "3hpa_H2", "3hpa_H3",
                                      "3hpa_F1", "3hpa_F2","3hpa_F3",
                                      "8hpa_H1", "8hpa_H2", "8hpa_H3", "8hpa_H4", "8hpa_H5",
                                      "8hpa_F1", "8hpa_F2", "8hpa_F3", "8hpa_F4", "8hpa_F5",
                                      "12hpa_H1", "12hpa_H2", "12hpa_H3", "12hpa_H4", "12hpa_H5",
                                      "12hpa_F1", "12hpa_F2", "12hpa_F3", "12hpa_F4", "12hpa_F5",
                                      "0hpa_iH1", "0hpa_iH2", "0hpa_iH3", "0hpa_iH4",
                                      "0hpa_iF1", "0hpa_iF2", "0hpa_iF3", "0hpa_iF4",
                                      "3hpa_iH1", "3hpa_iH2", "3hpa_iH3", "3hpa_iH4",
                                      "3hpa_iF1", "3hpa_iF2","3hpa_iF3", "3hpa_iF4", "3hpa_iF5",
                                      "8hpa_iH1", "8hpa_iH2", "8hpa_iH3", "8hpa_iH4", "8hpa_iH5",
                                      "8hpa_iF1", "8hpa_iF2", "8hpa_iF3", "8hpa_iF4",
                                      "12hpa_iH1", "12hpa_iH2", "12hpa_iH3", "12hpa_iH4",
                                      "12hpa_iF1", "12hpa_iF2", "12hpa_iF3", "12hpa_iF4", "12hpa_iF5"),
                         Condition = c(rep("Head_Uninjured", 5),
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
                                       rep("Foot_Regeneration_3hpa", 5),
                                       rep("Head_Regeneration_8hpa", 5),
                                       rep("Foot_Regeneration_8hpa", 4),
                                       rep("Head_Regeneration_12hpa", 4),
                                       rep("Foot_Regeneration_12hpa", 5)),
                         Treatment = c(rep("Head_Uninjured", 5),
                                       rep("Foot_Uninjured", 5),
                                       rep("Head_Regeneration_3hpa", 3),
                                       rep("Foot_Regeneration_3hpa", 3),
                                       rep("Head_Regeneration_8hpa", 5),
                                       rep("Foot_Regeneration_8hpa", 5),
                                       rep("Head_Regeneration_12hpa", 5),
                                       rep("Foot_Regeneration_12hpa", 5),
                                       rep("Head_Uninjured_inhib", 4),
                                       rep("Foot_Uninjured_inhib", 4),
                                       rep("Head_Regeneration_3hpa_inhib", 4),
                                       rep("Foot_Regeneration_3hpa_inhib", 5),
                                       rep("Head_Regeneration_8hpa_inhib", 5),
                                       rep("Foot_Regeneration_8hpa_inhib", 4),
                                       rep("Head_Regeneration_12hpa_inhib", 4),
                                       rep("Foot_Regeneration_12hpa_inhib", 5)),
                         Replicate = c(1:5,1:5,1:3,1:3,1:5,1:5,1:5,1:5,1:4,1:4,1:4,1:5,1:5,1:4,1:4,1:5),
                         bamReads = c("resources/Bams/0H1_final_shift.bam",
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
                                      "resources/Bams/12iF5_final_shift.bam"),
                         Peaks = rep("resources/full_consensus.bed",71),
                         PeakCaller = rep("bed",71)
)


#create object for analysis
regenATAC<- dba(sampleSheet = ATACSample)

#count reads within peaks
regenATAC<- dba.count(regenATAC, score = "DBA_SCORE_TMM_READS_FULL", bParallel = T)

fullPeaks <- regenATAC$peaks[[1]]
fullPeaks$peak_id <- rownames(fullPeaks)

writeBed(fullPeaks, "resources/full_consensus_diffbind_labels.bed")

save(regenATAC, file = "resources/full_ATAC_Counts.rds")

#build sampleSheet
ATACSample <- data.frame(SampleID = c("0hpa_H1", "0hpa_H2", "0hpa_H3", "0hpa_H4", "0hpa_H5",
                                      "0hpa_F1", "0hpa_F2", "0hpa_F3", "0hpa_F4", "0hpa_F5",
                                      "3hpa_H1", "3hpa_H2", "3hpa_H3",
                                      "3hpa_F1", "3hpa_F2","3hpa_F3",
                                      "8hpa_H1", "8hpa_H2", "8hpa_H3", "8hpa_H4", "8hpa_H5",
                                      "8hpa_F1", "8hpa_F2", "8hpa_F3", "8hpa_F4", "8hpa_F5",
                                      "12hpa_H1", "12hpa_H2", "12hpa_H3", "12hpa_H4", "12hpa_H5",
                                      "12hpa_F1", "12hpa_F2", "12hpa_F3", "12hpa_F4", "12hpa_F5"),
                         Condition = c(rep("Head_Uninjured", 5),
                                       rep("Foot_Uninjured", 5),
                                       rep("Head_Regeneration_3hpa", 3),
                                       rep("Foot_Regeneration_3hpa", 3),
                                       rep("Head_Regeneration_8hpa", 5),
                                       rep("Foot_Regeneration_8hpa", 5),
                                       rep("Head_Regeneration_12hpa", 5),
                                       rep("Foot_Regeneration_12hpa", 5)),
                         Treatment = c(rep("Head_Uninjured", 5),
                                       rep("Foot_Uninjured", 5),
                                       rep("Head_Regeneration_3hpa", 3),
                                       rep("Foot_Regeneration_3hpa", 3),
                                       rep("Head_Regeneration_8hpa", 5),
                                       rep("Foot_Regeneration_8hpa", 5),
                                       rep("Head_Regeneration_12hpa", 5),
                                       rep("Foot_Regeneration_12hpa", 5)),
                         Replicate = c(1:5,1:5,1:3,1:3,1:5,1:5,1:5,1:5),
                         bamReads = c("resources/Bams/0H1_final_shift.bam",
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
                                      "resources/Bams/12F5_final_shift.bam"),
                         Peaks = rep("resources/untreated_consensus.bed",36),
                         PeakCaller = rep("bed",36)
)


#create object for analysis
regenATAC<- dba(sampleSheet = ATACSample)

#count reads within peaks
regenATAC<- dba.count(regenATAC, score = "DBA_SCORE_TMM_READS_FULL", bParallel = T)

untreatedPeaks <- regenATAC$peaks[[1]]
untreatedPeaks$peak_id <- rownames(untreatedPeaks)

writeBed(untreatedPeaks, "resources/untreated_consensus_diffbind_labels.bed")

#annotate peaks
system("source resources/venv/bin/activate;./annotate_peaks.sh")

save(regenATAC, file = "resources/untreated_ATAC_Counts.rds")


