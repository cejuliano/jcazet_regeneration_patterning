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



####irt validation####

#tents
#recoreded at 60hpa
icrt <- read.csv("resources/icrt_tents.csv")

icrtTreat <- c(rep("DMSO", length(na.omit(icrt$DMSO))), rep("iCRT", length(na.omit(icrt$iCRT))))

icrt <- data.frame(treatment = icrtTreat, tents = c(na.omit(icrt$DMSO), icrt$iCRT))

gg <- ggplot(icrt, aes(treatment,tents)) + geom_boxplot(outlier.shape = NA)
gg <- gg + geom_jitter(width = 0.3, height = 0.2)
gg <- gg + theme_classic()
gg

ggsave(file = "plots/icrt_tents_plot.pdf", width = 2, height = 4.5, useDingbats=FALSE)

t.test(icrt$tents ~ icrt$treatment)

#feet
#recoreded at 36 hpa
icrt <- read.csv("resources/icrt_feet.csv")
icrt <- icrt[complete.cases(icrt),]

icrt$ave <- icrt$X..basal.disks/icrt$total.animals

icrt$treat <- c(rep("iCRT", 3), rep("DMSO", 3))

gg <- ggplot(icrt, aes(treat,ave)) + geom_boxplot(outlier.shape = NA)
gg <- gg + geom_jitter(width = 0.2)
gg <- gg + theme_classic()
gg

t.test(icrt$ave ~ icrt$treat)

ggsave(file = "plots/icrt_feet_plot.pdf", width = 2, height = 4.5, useDingbats=FALSE)

#icrt effect on skewering tents
icrt <- read.csv("resources/skewer_icrt.csv")

icrtTreat <- c(rep("DMSO", length(na.omit(icrt$DMSO))), rep("iCRT", length(na.omit(icrt$iCRT))))

icrt <- data.frame(treatment = icrtTreat, tents = c(na.omit(icrt$DMSO), icrt$iCRT))

gg <- ggplot(icrt, aes(treatment,tents)) + geom_boxplot(outlier.shape = NA)
gg <- gg + geom_jitter(width = 0.3, height = 0.2)
gg <- gg + theme_classic()
gg


t.test(icrt$tents ~ icrt$treatment)

ggsave(file = "plots/icrt_skewer_tents.pdf", width = 2.5, height = 4.5, useDingbats=FALSE)

#parallel skewer tents
#icrt effect on skewering tents
tents <- read.csv("resources/reamp_tents.csv")

tentsTreat <- c(rep("control", length(na.omit(tents$Control))), rep("reamp", length(na.omit(tents$Reamp))))

tents <- data.frame(treatment = tentsTreat, tents = c(tents$Control, na.omit(tents$Reamp)))

gg <- ggplot(tents, aes(treatment,tents)) + geom_boxplot(outlier.shape = NA)
gg <- gg + geom_jitter(width = 0.3, height = 0.2)
gg <- gg + theme_classic()
gg


t.test(tents$tents ~ tents$treatment)

ggsave(file = "plots/reamp_tents.pdf", width = 1.75, height = 4.5, useDingbats=FALSE)
