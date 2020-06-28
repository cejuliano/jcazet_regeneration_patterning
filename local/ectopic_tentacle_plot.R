library(ggplot2)
library(agricolae)

#Set working directory to be the location of this R script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#import ectopic tentacle quantification data
tents <- read.csv("resources/tentacle_quant.csv")

#merge time and treatment
tents$timeTreat <- factor(paste(tents$Organizer,tents$Time, sep = "_"))

#sub for strain comparison
tents.strain <- tents[!grepl("WO",tents$Organizer) & tents$Time != 0,]

#drop strains from main tent df
tents <- tents[!grepl("AEP|105", tents$timeTreat),]

#reorder factors
tents$timeTreat <- factor(tents$timeTreat)

print(levels(tents$timeTreat))

tents$timeTreat <- factor(tents$timeTreat,levels(tents$timeTreat)[c(3,4,1,2)])

gg <- ggplot(tents, aes(timeTreat,Tentacles)) + geom_boxplot(outlier.shape = NA)
gg <- gg + geom_jitter(width = 0.35, height = 0.1)
gg <- gg + theme_classic()
gg

ggsave(file = "plots/ectopic_tentacles.pdf", width = 2.5, height = 4.5, useDingbats=FALSE)

#anova
ajuste <- aov(tents$Tentacles ~ tents$timeTreat)
HSD.test(ajuste, 'tents$timeTreat', alpha = 0.05, console = T, group = T)


#strain plot
#reorder factors
tents.strain$timeTreat <- factor(tents.strain$timeTreat)

print(levels(tents.strain$timeTreat))

tents.strain$timeTreat <- factor(tents.strain$timeTreat,levels(tents.strain$timeTreat)[c(2,1,3)])

gg <- ggplot(tents.strain, aes(timeTreat,Tentacles)) + geom_boxplot(outlier.shape = NA)
gg <- gg + geom_jitter(width = 0.35, height = 0.1)
gg <- gg + theme_classic()
gg

ggsave(file = "plots/ectopic_tentacles_by_Strain.pdf", width = 2.5, height = 4.5, useDingbats=FALSE)
