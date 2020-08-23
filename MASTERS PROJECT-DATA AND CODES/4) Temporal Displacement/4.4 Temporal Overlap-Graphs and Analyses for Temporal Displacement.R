#PLOT GRAPHS OF TEMPORAL OVERLAP BETWEEN HUMANS AND MAMMALS
#IN LOW AND HIGH HUMAN ABUNDANCE SITES

#CALCULATE WALD STATISTIC TO GET A FORMAL STATISTICS DESCRIBING IF TEMPORAL OVERLAPS SIGNIFICANTLY DIFFER BETWEEN SITES
#WALD STATISTIC: W = (a1-a2)^2 / (SE1^2+SE2^2)
#WITH a1 AND a2 BEING THE ESTIMATES, AND SE1 AND SE2 BEING THEIR STANDARD ERRORS

rm(list=ls())

#INSTALL NECESSSARY PACKAGES
install.packages("ggplot2") #for plotting graphs
library("ggplot2")

##### LOW AND HIGH HUMAN ABUNDANCE #####
setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/4) Temporal Displacement") #set working directory to get temporal displacement results
HighLowOverlap<-read.csv("Temporal Overlap Results-High and Low Abundance.csv", header=TRUE) #load temporal displacement results
str(HighLowOverlap)
View(HighLowOverlap)

#GRAPH-HUMAN ABUNDANCE
#PLOT BAR CHART OF TEMPORAL OVERLAP ESTIMATE BY HUMAN ABUNDANCE
AnimalHLO<-subset(HighLowOverlap, Species == "Fox" | Species == "Hedgehog") #subset by mammals (no humans)
ggplot(AnimalHLO, aes(Species, Overlap, fill=Human.Abundance)) + geom_bar(stat="identity", position = "dodge", color="black")+labs(title="Temporal Overlap: Foxes and Hedgehogs", subtitle="Human Abundances", y="Temporal Overlap")+geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI), width=0.2, position=position_dodge(.9))+theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))+scale_fill_discrete(type=c("deepskyblue3","gold1"),name = "Human Abundance")
ggsave("Temporal Overlap-High and Low Humans-Foxes and Hedgehog.jpeg") #save JPEG. Adjust size of graphics window to ensure a clear plot

#WALD STATISTIC-HIGH AND LOW ABUNDANCE
#FOXES
#Wald statistics formula at top of script
Wald_Fox_HL<-((HighLowOverlap[3,5]-HighLowOverlap[4,5])^2)/(((HighLowOverlap[3,7]-HighLowOverlap[3,5])/1.96)^2+((HighLowOverlap[4,7]-HighLowOverlap[4,5])/1.96)^2)
Wald_Fox_HL #w=3.229351
#CHECK WALD STATISTICS AGAINST CHI SQUARE TABLE WITH ONE DEGREE OF FREEDOM
#IF P IS LESS THAN 0.05, TEMPORAL OVERLAP SIGNIFICANTLY DIFFERS BY HUMAN ABUNDANCE
pchisq(Wald_Fox_HL, df=1, lower.tail=FALSE) #p value= 0.07232935 #ALMOST SIGNIFICANT

#HEDGEHOG
Wald_Hog_HL<-((HighLowOverlap[5,5]-HighLowOverlap[6,5])^2)/(((HighLowOverlap[5,7]-HighLowOverlap[5,5])/1.96)^2+((HighLowOverlap[6,7]-HighLowOverlap[6,5])/1.96)^2)
Wald_Hog_HL #w= 0.6228303
pchisq(Wald_Hog_HL, df=1, lower.tail=FALSE) #p value= 0.4299975
