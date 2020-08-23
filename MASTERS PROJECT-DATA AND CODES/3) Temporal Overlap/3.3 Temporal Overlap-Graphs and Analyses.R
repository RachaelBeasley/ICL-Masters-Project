#PLOT GRAPHS OF TEMPORAL OVERLAP BETWEEN HUMANS AND MAMMALS
#IN DIFFERENT HABITAT COVERS AND TYPES, AND LOW AND HIGH HUMAN ABUNDANCE SITES

#CALCULATE WALD STATISTIC TO GET A FORMAL STATISTICS DESCRIBING IF TEMPORAL OVERLAPS SIGNIFICANTLY DIFFER BETWEEN SITES
#WALD STATISTIC: W = (a1-a2)^2 / (SE1^2+SE2^2)
#WITH a1 AND a2 BEING THE ESTIMATES, AND SE1 AND SE2 BEING THEIR STANDARD ERRORS


rm(list=ls())

#INSTALL NECESSSARY PACKAGES
install.packages("ggplot2") #for plotting graphs
library("ggplot2")

##### HABITAT COVER #####
setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/3) Temporal Overlap/Habitat Cover") #set working directory to get results
HCoverlap<-read.csv("Temporal Overlap Results-Habitat Cover.csv", header=TRUE) #load results from estimating temporal overlap by habitat cover
str(HCoverlap) #check data structure
View(HCoverlap) #visualise examine data

#GRAPH-HABITAT COVER
#PLOT BAR CHART OF TEMPORAL OVERLAP ESTIMATE BY HABITAT COVER
AnimalHCO<-subset(HCoverlap, Species == "Fox" | Species == "Hedgehog") #Subset mammals only (no humans)
ggplot(AnimalHCO, aes(Species, Overlap, fill=HabitatCover)) + geom_bar(stat="identity", position = "dodge", color="black") +labs(title="Temporal Overlap: Foxes and Hedgehogs", subtitle="Habitat Cover",y="Temporal Overlap")+geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI), width=0.2, position=position_dodge(.9))+theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))+ scale_fill_discrete(name = "Habitat Cover") 
ggsave("Temporal Overlap-Habitat Cover-Foxes and Hedgehog.jpeg") #save image, after adjusting graphic window to get a clear plot

#WALD STATISTIC-HABITAT COVER
#FOXES
#Wald statistic formula at top of script
Wald_Fox_HC<-((HCoverlap[3,5]-HCoverlap[4,5])^2)/(((HCoverlap[3,7]-HCoverlap[3,5])/1.96)^2+((HCoverlap[4,7]-HCoverlap[4,5])/1.96)^2)
Wald_Fox_HC #w=2.71963
#CHECK WALD STATISTICS AGAINST CHI SQUARE TABLE WITH ONE DEGREE OF FREEDOM
#IF P IS LESS THAN 0.05, TEMPORAL OVERLAP SIGNIFICANTLY DIFFERS BY HABITAT COVER
pchisq(Wald_Fox_HC, df=1, lower.tail=FALSE) #p value=0.09912097 

#HEDGEHOGS
Wald_Hog_HC<-((HCoverlap[5,5]-HCoverlap[6,5])^2)/(((HCoverlap[5,7]-HCoverlap[5,5])/1.96)^2+((HCoverlap[6,7]-HCoverlap[6,5])/1.96)^2)
Wald_Hog_HC #W=0.07394966
pchisq(Wald_Hog_HC, df=1, lower.tail=FALSE) #p value=0.7856706


##### HABITAT TYPE ######
setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/3) Temporal Overlap/Habitat Type") #set working directory to get results
HToverlap<-read.csv("Temporal Overlap Results-Habitat Type.csv", header=TRUE) #load results from estimating temporal overlap by habitat type
str(HToverlap)
View(HToverlap)

#GRAPHS-HABITAT TYPE
#PLOT BAR CHARTS OF TEMPORAL OVERLAP ESTIMATE BY HABITAT COVER
#DONE SEPARATELY FOR FOXES AND HEDGEHOGS FOR CLARITY

#FOXES
HTFox<-subset(HToverlap, Species=="Fox") #subset by foxes
ggplot(HTFox, aes(Habitat, Overlap)) + geom_bar(stat="identity", position = "dodge", fill="orchid3", color="black") +labs(title="Temporal Overlap: Foxes and Humans", subtitle="Habitat Type", y="Temporal Overlap", x="Habitat Type")+theme(axis.text.x = element_text(vjust=-0.05,angle = 90))+geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI), width=0.2, position=position_dodge(.9))+theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#Saved by exporting with a height of 400 and graphic window stretched wide enough for plot to be clear

#HEDGEHOGS
HTHedgehog<-subset(HToverlap, Species=="Hedgehog") #subset by hedgehogs
ggplot(HTHedgehog, aes(Habitat, Overlap)) + geom_bar(stat="identity", position = "dodge", fill="peru", color="black") +labs(title="Temporal Overlap: Hedgehog and Humans", subtitle="Habitat Type", y="Temporal Overlap", x="Habitat Type")+theme(axis.text.x = element_text(vjust=-0.05,angle = 90))+geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI), width=0.2, position=position_dodge(.9))+theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#N.B. Error for missing value occurs because there is no temporal overlap estimate for Tall Herbs
#Saved by exporting with a height of 400 and graphic window stretched wide enough for plot to be clear


#WALD STATISTIC-HABITAT TYPE

#FOXES
#CREATE RESULTS TABLE FOR WALD STATISTICS (EACH PAIRWISE COMPARISON BETWEEN HABITAT TYPES)
HTwald_Fox<-matrix(0, nrow=6, ncol=6) #empty results matrix
colnames(HTwald_Fox)<-sort(unique(HToverlap$HabitatType)) #Habitat Types become column names 
rownames(HTwald_Fox)<-sort(unique(HToverlap$HabitatType)) #Habitat Types become row names
HTwald_Fox<-HTwald_Fox[,-6] #Remove last column to prevent unnecessary recalculation of wald statistics for Tall Herbs
View(HTwald_Fox) #check results matrix is formatted correctly

#CREATE ASSOCIATED RESULTS TABLE FOR P VALUES
HTpvalue_Fox<-matrix(0, nrow=6, ncol=6) #empty results matrix
colnames(HTpvalue_Fox)<-sort(unique(HToverlap$Habitat)) #Habitat Types become column names
rownames(HTpvalue_Fox)<-sort(unique(HToverlap$Habitat)) #Habitat Types become row names
HTpvalue_Fox<-HTpvalue_Fox[,-6] #Remove last column to prevent unnecessary recalculation of wald statistics for Tall Herbs
View(HTpvalue_Fox) #check results matrix is formatted correctly

#CALCULATE WALD TEST AND P VALUES FOR FOXES
  for(i in 1:5){ #for each habitat type (each habitat types takes a turn being the reference one) 
    for(j in (i+1):6){ #for each habitat type after the current reference habitat type
      #calculate the wald statistic between reference habitat type and each habitat type after it, and save in wald statistics results matrix
      #Wald statistic formula at top of script
      HTwald_Fox[j,i]<-((HToverlap[i+6,5]-HToverlap[j+6,5])^2)/(((HToverlap[i+6,7]-HToverlap[i+6,5])/1.96)^2+((HToverlap[j+6,7]-HToverlap[j+6,5])/1.96)^2) 
      #determine associated p values and save in p-value results matrix
      HTpvalue_Fox[j,i]<-pchisq(HTwald_Fox[j,i], df=1, lower.tail=FALSE) 
      }
  }

#HEDGEHOGS
#CREATE RESULTS TABLE FOR WALD STATISTICS (EACH PAIRWISE COMPARISON BETWEEN HABITAT TYPES)
HTwald_Hog<-matrix(0, nrow=6, ncol=6) #empty results matrix
colnames(HTwald_Hog)<-sort(unique(HToverlap$Habitat)) #Habitat Types become column names
rownames(HTwald_Hog)<-sort(unique(HToverlap$Habitat)) #Habitat Types become row names
HTwald_Hog<-HTwald_Hog[,-6] #Remove last column to prevent unnecessary recalculation of wald statistics for Tall Herbs
View(HTwald_Hog) #check results matrix is formatted correctly

#CREATE ASSOCIATED RESULTS TABLE FOR P VALUES
HTpvalue_Hog<-matrix(0, nrow=6, ncol=6) #empty results matrix
colnames(HTpvalue_Hog)<-sort(unique(HToverlap$Habitat)) #Habitat Types become column names
rownames(HTpvalue_Hog)<-sort(unique(HToverlap$Habitat)) #Habitat Types become row names
HTpvalue_Hog<-HTpvalue_Hog[,-6] #Remove last column to prevent unnecessary recalculation of wald statistics for Tall Herbs
View(HTpvalue_Hog) #check results matrix is formatted correctly

#CALCULATE WALD TEST AND P VALUES FOR HEDGEHOGS
for(i in 1:5){ #for each habitat type (each habitat types takes a turn being the reference one) 
  for(j in (i+1):6){ #for each habitat type after the current reference habitat type
    #calculate the wald statistic between reference habitat type and each habitat type after it, and save in wald statistics results matrix
    #Wald statistic formula at top of script
    HTwald_Hog[j,i]<-((HToverlap[i+12,5]-HToverlap[j+12,5])^2)/(((HToverlap[i+12,7]-HToverlap[i+12,5])/1.96)^2+((HToverlap[j+12,7]-HToverlap[j+12,5])/1.96)^2)
    #determine associated p values and save in p-value results matrix
    HTpvalue_Hog[j,i]<-pchisq(HTwald_Hog[j,i], df=1, lower.tail=FALSE)
  }
}

#SAVE ALL RESULTS MATRICES IN THE SAME WORKING DIRECTORY
write.csv(HTwald_Fox, "Habitat Type-Wald Statistic-FOX.csv")
write.csv(HTpvalue_Fox, "Habitat Type-P Value-FOX.csv")
write.csv(HTwald_Hog, "Habitat Type-Wald Statistic-HEDGEHOG.csv")
write.csv(HTpvalue_Hog, "Habitat Type-P Value-HEDGEHOG.csv")




