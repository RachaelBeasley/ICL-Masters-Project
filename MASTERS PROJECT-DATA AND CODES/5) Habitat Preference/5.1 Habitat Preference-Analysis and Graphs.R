#HABITAT PREFERENCE
#PLOT RELATIVE ABUNDANCE INDICES BY HABITAT WITH ERROR BARS
#CONDUCT STATISTICAL ANALYSIS TO SEE IF SIGNIFICANCE DIFFERENCE IN RAIS ACROSS SITES OF DIFFERENT HABITAT COVER AND TYPE

rm(list=ls())

setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/4) Temporal Displacement") #set working directory to retrieve results data
RAIResults<-read.csv("RAI Results MD EXIF-Two Minutes.csv" , header=TRUE) #load Relative Abundance Indices Data
View(RAIResults) #N.B. Place ID column name added in Excel
str(RAIResults)

#RELATIVE ABUNDANCES FROM EACH SITE OBTAINED 
#BUT NEED TO BE GROUPED BY HABITAT TYPE/COVER
#ADD HABITAT DATA
setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/Data") #set working directory to retrieve data
Habitats<-read.csv("Site Habitat.csv", header=TRUE) #load data with Habitat Types and Covers by camera site
View(Habitats)
str(Habitats)

#MERGE DATASETS
RAIResults<-merge(RAIResults, Habitats, by.x="PlaceID", by.y="TagPlaceID" )

####### GROUP BY HABITAT COVER ########
Closed<-subset(RAIResults, HabitatCover=="Closed")
length(Closed$PlaceID) #61 sites with Closed Habitat Cover
Open<-subset(RAIResults, HabitatCover=="Open")
length(Open$PlaceID) #53 sites with Open Habitat Cover

#EXAMINE DATA, CHECK STATISTICAL TEST ASSUMPTIONS AND DO ANALYSES
#DETERMINE WHETHER ASSUMPTIONS OF PARAMETRIC TESTS ARE MET (NORMALITY AND HOMOGENEITY OF VARIANCES) 

#INSTALL NECESSARY PACKAGES
install.packages("lme4") 
library("lme4")#To run mixed linear models
install.packages("predictmeans") 
library("predictmeans")#To see residual plots

###### HUMANS ######

#CLOSED HABITAT COVER
par(mar=c(4,4,1,1)) #change graphic window margins to allow figures to be plotted
hist(log(Closed$Human+1)) #Plot histogram of logged human RAI in Closed Habitats-positively Skewed
boxplot(log(Closed$Human+1)) #Plot boxplot of logged human RAI in Closed Habitats-No outliers
shapiro.test(log(Closed$Human+1)) #Statistically test if data is normally distributed using Shapiro Test
#Significant (p<0.05), therefore data is not normally distributed

#OPEN HABITAT COVER
hist(log(Open$Human+1))#Plot histogram of logged human RAI in Open Habitats-very slight positive Skew, but fairly well distributed
boxplot(log(Open$Human+1)) #Plot boxplot of logged human RAI in Open Habitats-No outliers
shapiro.test(log(Open$Human+1)) #Statistically test if data is normally distributed using Shapiro Test
#Not significant (p=0.02), therefore data is normally distributed

#CHECKING HOMOGENEITY OF VARIANCES
var(log(Open$Human+1))/var(log(Closed$Human+1)) #Ratio of variances is 1.63345. As this is less than four, variances are homogeneous.

#CONCLUSION: 
#DESPITE HOMOGENEOUS OF VARIANCES AND OPEN HABITAT COVER DATA BEING NORMALLY DISTRIBUTED,
#DATA FOR CLOSED HABITATS IS NOT NORMALLY DISTRIBUTED.
#THIS SUGGESTS A NON-PARAMETRIC TEST MUST BE USED TO TEST FOR SIGNIFICANT DIFFERENCES BETWEEN RAIS IN OPEN AND CLOSED HABITATS.
#I.E.USE A WILCOX TEST IN PLACE OF A T-TEST (NONPARAMETRIC ALTERNATIVE)

wilcox.test(log(Closed$Human+1), log(Open$Human+1)) #p.value<0.05. Significant difference
median(log(Closed$Human+1)) #Median in Closed Habitats: 0.7884574
median(log(Open$Human+1)) #Median in Open Habitats: 2.09404

#RESULTS: Significant difference; Humans are sighted more in open habitats than closed ones

#IMPORTANT TO ACCOUNT FOR HOW SURVEYS ACROSS OPEN AND CLOSED SITES OCCURRED AT DIFFERENT TIMES
#THEREFORE ALSO CHECK RESULT USING A MIXED LINEAR MODEL

#ADD MONTH TO DATASET
#LOAD SURVEY LENGTHS DATASET
#N.B. THE MAIN MONTH OF EACH CAMERA TRAP'S SURVEY HAS BEEN ADDED TO THE "SURVEY LENGTH.CSV" RESULTS FILE
#AS HAS THE DETECTION ZONE
setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/2) Miscellaneous") #set working directory to retrieve results
SL<-read.csv("Survey Lengths-Month and Detection Zone Added.csv", header=TRUE) #load Survey Length data with month and detection zone manually added to it
View(SL)
SL<-SL[,-c(1,8,9)] #Remove surplus columns; reduces the number of redundant columns in the dataset

#MERGE DATASETS INTO FINAL ONE FOR ANALYSIS
FINAL<-merge(RAIResults, SL, by.x="PlaceID", by.y="PlaceID")

mH<-lmer(log(Human+1)~HabitatCover+(1|Month), data=FINAL)
summary(mH) 
#As t-value is greater than two, significant difference between abundances in open and closed habitats is confirmed
residplot(mH)
#Examining residuals show that correcting for month results in greater normality (Normal QQ plot is fairly diagonal)
#Residuals are fairly homogeneous in variance, as shown by similar spread of residuals around y=0

###### FOXES ######
#CLOSED HABITAT
hist(log(Closed$Fox+1)) #Plot histogram of logged fox RAI in Closed Habitats-positively skewed
boxplot(log(Closed$Fox+1)) #Plot boxplot of logged fox RAI in Closed Habitats-four outliers
shapiro.test(log(Closed$Fox+1)) #Statistically test if data is normally distributed using Shapiro Test
#Significant (p<0.05), therefore data is not normally distributed

#OPEN HABITAT
hist(log(Open$Fox+1)) #Plot histogram of logged fox RAI in Open Habitats-slightly positively skewed
boxplot(log(Open$Fox+1)) #Plot boxplot of logged fox RAI in Open Habitats-One outlier
shapiro.test(log(Open$Fox+1)) #Statistically test if data is normally distributed using Shapiro Test
#Significant (p<0.05), therefore data is not normally distributed

#Check for homogeneity of variances:
var(log(Closed$Fox+1))/var(log(Open$Fox+1)) #Ratio of variances is 0.8672056. As it is less than four, variances are homogeneous.

#CONCLUSION: ASSUMPTION OF NORMALITY IS BROKEN FOR BOTH OPEN AND CLOSED HABITATS. 
#THEREFORE, DESPITE HOMOGENEOUS OF VARIANCES, A NON-PARAMETRIC TEST IS THE APPROPRIATE TEST TO USE
#I.E.USE A WILCOX TEST IN PLACE OF A T-TEST (NONPARAMETRIC ALTERNATIVE)

wilcox.test(log(Closed$Fox+1), log(Open$Fox+1)) #p-value=0.7547. Not-significant
median(log(Closed$Fox+1)) #0.5389965
median(log(Open$Fox+1)) #0.5260931
#RESULTS: Not significantly different. Fox RAIs do not vary between open or closed habitats

#IMPORTANT TO ACCOUNT FOR HOW SURVEYS ACROSS OPEN AND CLOSED SITES OCCURRED AT DIFFERENT TIMES
#THEREFORE ALSO CHECK RESULT USING A MIXED LINEAR MODEL
mF<-lmer(log(Fox+1)~HabitatCover+(1|Month), data=FINAL)
summary(mF) 
#As t-value is less than two, the result of no significant difference between abundances in open and closed habitats is confirmed
residplot(mF)
#Examining residuals show that correcting for month still results in some deviation from normality
#(Normal QQ plot is fairly diagonal, but curves at the start and end of the plot, reducing the amount of normality)
#Variance are fairly homogeneous with central line having smaller spread than the others on the residual-fitted plots
#Suggests wilcox.test was appropriate to use as slight deviations from the assumptions
#but these deviations are not too extreme to invalidate the results of the linear model


##### HEDGEHOGS #####
#CLOSED HABITAT
hist(log(Closed$Hedgehog+1)) #Plot histogram of logged hedgehog RAI in Closed Habitats-positively skewed
boxplot(log(Closed$Hedgehog+1)) #Plot boxplot of logged hedgehog RAI in Closed Habitats-four outliers
shapiro.test(log(Closed$Hedgehog+1)) #Statistically test if data is normally distributed using Shapiro Test
#Significant (p<0.05), therefore data is not normally distributed

#OPEN HABITAT
hist(log(Open$Hedgehog+1)) #Plot histogram of logged hedgehog RAI in open Habitats-positively skewed
boxplot(log(Open$Hedgehog+1)) #Plot boxplot of logged hedgehog RAI in Open Habitats-eight outliers
shapiro.test(log(Open$Hedgehog+1)) #Statistically test if data is normally distributed using Shapiro Test
#Significant (p<0.05), therefore data is not normally distributed

#Check for homogeneity of variances:
var(log(Closed$Hedgehog+1))/var(log(Open$Hedgehog+1)) #Ratio of variances is 1.119226. As it is less than four, variances are homogeneous.

#CONCLUSION: ASSUMPTION OF NORMALITY IS BROKEN FOR BOTH OPEN AND CLOSED HABITATS. 
#THEREFORE, DESPITE HOMOGENEOUS OF VARIANCES, A NON-PARAMETRIC TEST IS THE APPROPRIATE TEST TO USE
#I.E.USE A WILCOX TEST IN PLACE OF A T-TEST (NONPARAMETRIC ALTERNATIVE)

wilcox.test(log(Closed$Hedgehog+1), log(Open$Hedgehog+1)) #p-value=0.478
median(log(Closed$Hedgehog+1)) #0.06899287
median(log(Open$Hedgehog+1)) #0.04255961
#RESULTS: No significant difference found between hedgehog RAIs in open or close habitats

#IMPORTANT TO ACCOUNT FOR HOW SURVEYS ACROSS OPEN AND CLOSED SITES OCCURRED AT DIFFERENT TIMES
#THEREFORE ALSO CHECK RESULT USING A MIXED LINEAR MODEL
mHH<-lmer(log(Hedgehog+1)~HabitatCover+(1|Month), data=FINAL)
summary(mHH) 
#Month explains no variance, so do not need to correct for it.


##### GROUP BY HABITAT TYPE #####
#SIX HABITAT TYPE CATEGORIES
#SEE DISTRIBTUTION OF HABITAT TYPES ACROSS SITES
table(RAIResults$HabitatType) #range is 13 to 26 sites per habitat type. Unbalanced sample size.

#FIRST DECIDE IF A PARAMETRIC TEST OR NOT IS REQUIRED
#BY TESTING ASSUMPTIONS OF MODEL
install.packages("ggpubr") #to plot boxplots in order to examine variance and see outliers
library("ggpubr") 
install.packages("car") #to test the homogeniety of variances using the levene test
library(car) 

##### HUMANS ####

#CHECK FOR OUTLIERS
ggboxplot(RAIResults, x = "HabitatType", y = "Human") #Plot boxplots of human RAI by habitat type. Lots of outliers, and variances appear quite unequal
#CHECK STATISTICALLY FOR HOMOGENEITY OF VARIANCES
leveneTest(log(Human+1) ~ HabitatType, data = RAIResults) #Variances significantly different from eachother (p<0.05), therefore variances are statistically not homogeneous
#CHECK NORMALITY
res.aovH <- aov(log(Human+1) ~ HabitatType, data = RAIResults) #generate residuals from comparing human RAI by habitat types
aov_residuals <- residuals(object = res.aovH ) #extract residual values from above model
par(mar=c(4,4,1,1)) #change graphics window to allow figure to fit
qqnorm(aov_residuals) #Trend is sigmoidal, though many points appear to lie on the diagonal. Deviations at lower and upper corners. Not ideally normal, but likely acceptable.
#CONCLUSION: ASSUMPTION OF HOMOGENEITY OF VARIANCES (AND POSSIBLY) NORMALITY ARE NOT MET, AND THERE IS A VARIETY OF SAMPLE SIZES
#THIS SUGGESTS A NON-PARAMETRIC TEST SHOULD BE USED
#I.E. A KRUSKAL WALLIS TEST INSTEAD OF AN ANOVA.

kruskal.test(log(Human+1) ~ HabitatType, data = RAIResults) 
#Relative sightings by habitat type are significantly different (p<0.05), therefore do a posthoc test to determine which habitat types are significantly different from each other

#Post Hoc test: Dunn Test
install.packages("dunn.test") #package needed to conduct a dunn test
library("dunn.test")

dunn.test(log(RAIResults$Human+1), RAIResults$HabitatType)
#p value has been adjusted for number of comparisons, so is 0.025
#See results summary for abbreviations' full names
#Significant difference between AG and AS (p=0.0040), BLW (p=0.0019), and DS (p<0.0001)
#Significant difference between AS and DS (p=0.0076) and  HR (0.0202)
#Significant difference between BLW and DS (p=0.0035), and ***HR (p=0.0107)***
#Significant difference between DS and HR (p<0.001) and TH (p=0.0003)

#IMPORTANT TO ACCOUNT FOR HOW SURVEYS ACROSS OPEN AND CLOSED SITES OCCURRED AT DIFFERENT TIMES
#THEREFORE ALSO CHECK RESULT USING A MIXED LINEAR MODEL
HT_Human<-lmer(log(Human+1)~HabitatType+(1|Month), data=FINAL) #construct simple linear model with logged RAIs as the response varibles, habitat type as a fixed effect and month as a random effect
summary(HT_Human) #As some t-values are greater than two, need to do posthoc pairwise comparisons (see below)
residplot(HT_Human) #Normality plot is fairly straight with some deviations, and variances are fairly homogeneous.
#assumptions are acceptably met to validate the results of the linear model

#NEXT DETERMINE WHICH PAIRWISE HABITAT TYPES ARE SIGNIFICANTLY DIFFERENT FROM EACH OTHER
#PAIRWISE COMPARISONS
install.packages("emmeans") #package to compare estimates and gradients between categories pairwise
library("emmeans")

emmeans(HT_Human, pairwise ~ HabitatType) #pairwise comparison of RAIs between habitat type
#Significance difference declared if <0.05
#Significant difference between AG and BLW, DS (ALMOST AS)
#Significant difference between DS and HR, TH
#Correcting for month leads to some significant differences between habitats being lost,relative to the dunn.test. 
#The significant differences that remain are validated by correcting for month
#However, correcting for month makes a different, and assumptions of linear model is acceptably met
#Therefore overall results from linear model correcting for month have greater weight than non-parametric test


##### FOXES #####
#CHECK FOR OUTLIERS
ggboxplot(RAIResults, x = "HabitatType", y = "Fox") #Plot boxplots of fox RAIs by habitat type. Seven outliers, and variances appear unequal between some groups
#CHECK STATISTICALLY FOR HOMOGENEITY OF VARIANCES
leveneTest(log(Fox+1) ~ HabitatType, data = RAIResults) #Variances not significantly different from eachother, therefore variances are statistically homogeneous
#CHECK NORMALITY
res.aovF <- aov(log(Fox+1) ~ HabitatType, data = RAIResults) #generate residuals from comparing fox RAI by habitat types
aov_residualsF <- residuals(object = res.aovF ) #extract residual values from above model
qqnorm(aov_residualsF) #Trend is slightly curved, across whole plot. Appears that at least half of the points would not lie on the diagonal. Not normal.

#CONCLUSION: ASSUMPTION OF NORMALITY IS NOT MET, AND THERE IS A VARIETY OF SAMPLE SIZES
#THIS SUGGESTS A NON-PARAMETRIC TEST SHOULD BE USED
#I.E. A KRUSKAL WALLIS TEST INSTEAD OF AN ANOVA.

kruskal.test(log(Fox+1) ~ HabitatType, data = RAIResults) 
#Relative Abundance by habitat type are not significantly different.
#As not significant, post hoc tests not required

#IMPORTANT TO ACCOUNT FOR HOW SURVEYS ACROSS OPEN AND CLOSED SITES OCCURRED AT DIFFERENT TIMES
#THEREFORE ALSO CHECK RESULT USING A MIXED LINEAR MODEL
HT_Fox<-lmer(log(Fox+1)~HabitatType+(1|Month), data=FINAL) #construct simple linear model with logged RAIs as the response variables, habitat type as a fixed effect and month as a random effect
summary(HT_Fox) 

emmeans(HT_Fox, pairwise ~ HabitatType) #pairwise comparison of RAIs between habitat type
#No t-values are greater than two. Validates kruskal wallis test results.
residplot(HT_Fox)
#normal QQ plot: Deviates quite a lot from diagonal in upper right corner. Fairly straight elsewhere. Not completely normal, but acceptable for ecological data
#Residual fitted has good amount of scatter for categorical variables.
#Generally, diagnostics validate results of linear model, which complement the results of the kruskal wallis test

##### HEDGEHOGS #####

#CHECK FOR OUTLIERS
ggboxplot(RAIResults, x = "HabitatType", y = "Hedgehog") #Plot boxplots of hedgehog RAIs by habitat type. Nine outliers, and variances appear quite fairly equal except in Tall Herbs
#CHECK STATISTICALLY FOR HOMOGENEITY OF VARIANCES
leveneTest(log(Hedgehog+1) ~ HabitatType, data = RAIResults) #Variances not significantly different from eachother, therefore variances are statistically homogeneous
#CHECK NORMALITY
res.aovHH <- aov(log(Hedgehog+1)~ HabitatType, data = RAIResults) #generate residuals from comparing hedgehog RAI by habitat types
aov_residualsHH <- residuals(object = res.aovHH ) #extract residual values from above model
qqnorm(aov_residualsHH) #Trend is very curved, especially in upper right hand corner Not normal.

#CONCLUSION: ASSUMPTION OF NORMALITY IS NOT MET, AND THERE IS A VARIETY OF SAMPLE SIZES
#THIS SUGGESTS A NON-PARAMETRIC TEST SHOULD BE USED
#I.E. A KRUSKAL WALLIS TEST INSTEAD OF AN ANOVA.

kruskal.test(Hedgehog ~ HabitatType, data = RAIResults) 
#Relative Abundance by habitat type are not significantly different.
#As not significant, post hoc tests not required

#IMPORTANT TO ACCOUNT FOR HOW SURVEYS ACROSS OPEN AND CLOSED SITES OCCURRED AT DIFFERENT TIMES
#THEREFORE ALSO CHECK RESULT USING A MIXED LINEAR MODEL
HT_Hog<-lmer(log(Hedgehog+1)~HabitatType+(1|Month), data=FINAL) #construct simple linear model with logged RAIs as the response variables, habitat type as a fixed effect and month as a random effect
summary(HT_Hog)
#Month explains no variance, so do not need to correct for it.

emmeans(HT_Hog, pairwise ~ HabitatType) #pairwise comparison of RAIs between habitat type
#No t-values are greater than two. Validates kruskal wallis test results.
residplot(HT_Hog)
#normal QQ plot: Deviates quite a lot from diagonal especially in upper right corner. Fairly straight elsewhere. Residuals are not very normal.
#Residual fitted has poor amount of scatter compressed to right hand side
#Results of linear model, complement kruskal wallis results, but their weight is reduced by their poor diagnostics



#### GENERATE SUMMARY OF MEDIAN RAI VALUES BY HABITAT #####
### WITH GRAPHS-ALTERNATIVE WAY OF VISUALING THE DATA

#### HABITAT COVER ####
#CONSTRUCT TABLE TO PLOT GRAPHS
#MEDIANS
ResultsToPlotMED<-as.data.frame(matrix(0, ncol=5, nrow=6))
colnames(ResultsToPlotMED)<-c("Species", "MedianLog", "LowerCI", "UpperCI", "HabitatCover")
ResultsToPlotMED[,1]<-rep(c("Human", "Fox", "Hedgehog"), rep=2)
ResultsToPlotMED[,5]<-rep(c("Closed", "Open"), each=3)

#Add median values
Closed<-subset(RAIResults, HabitatCover=="Closed")
Open<-subset(RAIResults, HabitatCover=="Open")
ResultsToPlotMED[1:3,2]<-apply(Closed[,8:10], 2, median)
ResultsToPlotMED[4:6,2]<-apply(Open[,8:10], 2, median)

#ADD CONFIDENCE INTERVALS
#from: https://www.ucl.ac.uk/child-health/short-courses-events/about-statistical-courses/research-methods-and-statistics/chapter-8-content-8

LowerCIrankedvalueCLOSED<-(length(Closed$PlaceID)/2)-(1.96*sqrt(length(Closed$PlaceID))/2) #22.84596 rounds to 23
UpperCIrankedvalueCLOSED<-1+(length(Closed$PlaceID)/2)+(1.96*sqrt(length(Closed$PlaceID)/2)) #42.32445 rounds to 42

LowerCIrankedvalueOPEN<-(length(Open$PlaceID)/2)-(1.96*sqrt(length(Open$PlaceID))/2) #19.36549 rounds to 19
UpperCIrankedvalueOPEN<-1+(length(Open$PlaceID)/2)+(1.96*sqrt(length(Open$PlaceID)/2)) #37.58972 rounds to 38

#FIND CONFIDENCE INTERVALS AND ADD TO THE TABLE
for(i in 8:10){
  c<-sort(Closed[,i])
  o<-sort(Open[,i])
  ResultsToPlotMED[i-7, 3]<-c[23]
  ResultsToPlotMED[i-7, 4]<-c[42]
  ResultsToPlotMED[i-4, 3]<-o[19]
  ResultsToPlotMED[i-4, 4]<-o[38]
}

View(ResultsToPlotMED)
setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/5) Habitat Preference") #set working directory to save results
write.csv(ResultsToPlotMED, "Habitat Cover-Median Logged RAIs.csv") #save completed results matrix

#PLOT GRAPHS FOR HABITAT COVER
#WITH MEDIAN VALUES

#WITHOUT HUMANS
ggplot(ResultsToPlotMED[-c(1,4),], aes(Species, MedianLog, fill=HabitatCover)) + geom_bar(stat="identity", position = "dodge", color="black") +labs(title="Species Sightings by Habitat Cover", y="Median Log of RAI")+theme(axis.text.x = element_text(vjust=-0.05,angle = 90))+ scale_fill_discrete(name = "Habitat Cover")+geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI), width=.2,position=position_dodge(.9))
ggsave("Median Logged RAI by Habitat Cover-Mammal.jpeg")

#WITH HUMANS
ggplot(ResultsToPlotMED[c(1,4),], aes(Species, MedianLog, fill=HabitatCover)) + geom_bar(stat="identity", position = "dodge", color="black") +labs(title="Human Sightings by Habitat Cover", y="Median Log of RAI")+theme(axis.text.x = element_text(vjust=-0.05,angle = 90))+ scale_fill_discrete(name = "Habitat Cover")+geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI), width=.2,position=position_dodge(.9))
ggsave("Median Logged RAI by Habitat Cover-Human.jpeg")

#### HABITAT TYPE ####
#WITH MEDIANS

#CONSTRUCT TABLE FOR VALUES TO PLOT
HabType<-sort(unique(Habitats$HabitatType))
ResultsToPlotMEDHT<-as.data.frame(matrix(0, ncol=5, nrow=18))
colnames(ResultsToPlotMEDHT)<-c("Species", "MedianLog", "LowerCI", "UpperCI", "HabitatType")
ResultsToPlotMEDHT[,1]<-rep(c("Human", "Fox", "Hedgehog"), rep=6)
ResultsToPlotMEDHT[,5]<-rep(HabType, each=3)

#ADD MEDIAN VALUES
for(i in 1:6){
  HT<-subset(RAIResults, HabitatType==HabType[i])
  rows<-c(1,2,3)+3*(i-1)
  ResultsToPlotMEDHT[rows,2]<-apply(HT[,8:10], 2, median) #add median to results
}

#ADD CONFIDENCE INTERVALS 
for(i in 1:6){
  HT<-subset(RAIResults, HabitatType==HabType[i])
  for(j in 8:10){
    d<-sort(HT[,j]) #select that species
    LowerCIrankedvalue<-(length(d)/2)-((1.96*sqrt(length(d)))/2)
    UpperCIrankedvalue<-1+(length(d)/2)+((1.96*sqrt(length(d)))/2)
    if(LowerCIrankedvalue<0.5){LowerCIrankedvalue<-1}
    if(UpperCIrankedvalue>length(d)+0.5){UpperCIrankedvalue<-length(d)}
    rows<-(j-7)+3*(i-1)
    ResultsToPlotMEDHT[rows, 3]<-d[round(LowerCIrankedvalue)]
    ResultsToPlotMEDHT[rows, 4]<-d[round(UpperCIrankedvalue)]
  }
}

setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/5) Habitat Preference") #set working directory to save results
write.csv(ResultsToPlotMEDHT, "Habitat Type-Median Logged RAIs.csv")

#PLOT GRAPHS
#WITH MEDIANS

#WITH FOXES
Fox<-subset(ResultsToPlotMEDHT, Species =="Fox")
ggplot(Fox, aes(HabitatType, MedianLog)) + geom_bar(stat="identity", position = "dodge", color="black", fill="orchid3") +labs(title="Fox Sightings by Habitat Type", y="Median Log of RAI")+theme(axis.text.x = element_text(vjust=-0.05,angle = 90), plot.title=element_text(hjust = 0.5))+geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI), width=.2,position=position_dodge(.9))
ggsave("Median Logged RAI by Habitat Type-Fox.jpeg")

#WITH HEDGEHOGS
Hedgehog<-subset(ResultsToPlotMEDHT, Species =="Hedgehog")
ggplot(Hedgehog, aes(HabitatType, MedianLog)) + geom_bar(stat="identity", position = "dodge", color="black", fill="peru") +labs(title="Hedgehog Sightings by Habitat Type", y="Median Log of RAI")+theme(axis.text.x = element_text(vjust=-0.05,angle = 90), plot.title=element_text(hjust = 0.5))+geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI), width=.2,position=position_dodge(.9))
ggsave("Median Logged RAI by Habitat Type-Hedgehog.jpeg")

#WITH HUMANS
Human<-subset(ResultsToPlotMEDHT, Species =="Human")
ggplot(Human, aes(HabitatType, MedianLog)) + geom_bar(stat="identity", position = "dodge", color="black", fill="lightgoldenrod") +labs(title="Human Sightings by Habitat Type", y="Median Log of RAI")+theme(axis.text.x = element_text(vjust=-0.05,angle = 90), plot.title=element_text(hjust = 0.5))+geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI), width=.2,position=position_dodge(.9))
ggsave("Median Logged RAI by Habitat Type-Human.jpeg")



#### PLOT GRAPHS ####
install.packages("ggplot2")
library(ggplot2)

#LOG VALUES WITHIN DATASET TO SIMPLIFY GRAPHING FUNCTION
RAIResults$HumanLog<-log(RAIResults$Human+1)
RAIResults$FoxLog<-log(RAIResults$Fox+1)
RAIResults$HedgehogLog<-log(RAIResults$Hedgehog+1)

####### GRAPHS WITH HABITAT COVER #######

##### PLOT BOXPLOTS #####
install.packages("tidyr") #package so that dataset can be converted from a wide to a narrow dataset. required for ggplot to interpret the data
library("tidyr")

setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/5) Habitat Preference") #set working directory to save graphs

RAIResultsF<-RAIResults[,6:10] #filter only the data that is required (logged RAIs and habitat types and covers)
gatheredDATA<-gather(RAIResultsF, key="Species", value="LogRAI", HumanLog, FoxLog, HedgehogLog) #convert dataset into a narrow dataframe so that ggplot can read it correctly
gatheredDATA_species<-subset(gatheredDATA, Species !="HumanLog") #remove humans from the dataset for the mammal only plots

ticklabels<-c("Fox", "Hedgehog", "Human") #define tick labels in order to change x axis tick labels from "SpeciesLog" to "Species"
#Plot boxplot of logged mammal RAIs
ggplot(gatheredDATA_species, aes(x=Species, y=LogRAI,fill=HabitatCover)) + geom_boxplot()+labs(title="Species RAIs by Habitat Cover", y="Logged Relative Abundance Index")+theme(plot.title = element_text(hjust = 0.5))+scale_fill_discrete(name = "Habitat Cover")+ scale_x_discrete(labels= ticklabels)
ggsave("Boxplot-Mammal Logged RAIs-Habitat Cover.jpeg") #save graph as JPEG. ensure graphics windows is adjusted to make figure clear

gatheredDATA_human<-subset(gatheredDATA, Species =="HumanLog") #filter data to have logged human RAIs only
#Plot boxplot of logged mammal RAIs
ggplot(gatheredDATA_human, aes(x=Species, y=LogRAI,  fill=HabitatCover)) + geom_boxplot()+labs(title="Human RAIs by Habitat Cover", y="Logged Relative Abundance Index", x="Humans")+theme(plot.title = element_text(hjust = 0.5))+scale_fill_discrete(name = "Habitat Cover")
ggsave("Boxplot-Human RAIs Log-Habitat Cover.jpeg") #save graph as JPEG. ensure graphics windows is adjusted to make figure clear


#### PLOT BARCHARTS WITH MEANS ####
#Not presented in write up, but an alternative way to view data
#CONSTRUCT TABLE FOR MEAN VALUES AND CONFIDENCE INTERVALS TO PLOT
ResultsToPlotMEAN<-matrix(0, nrow=6, ncol=5) #empty results matrix
colnames(ResultsToPlotMEAN)<-c("Species", "MeanOfLog", "LowerCI", "UpperCI", "HabitatCover") #add column names to results matrix
ResultsToPlotMEAN[,1]<-rep(c("Human", "Fox", "Hedgehog"), rep=2) #add species names to results matrix
ResultsToPlotMEAN[,5]<-rep(c("Closed", "Open"), each=3) #add habitit cover names to results matrix

#CALCULATE MEAN LOGGED RAIS
Closed<-subset(RAIResults, HabitatCover=="Closed") #subset data by Closed habitat cover
Open<-subset(RAIResults, HabitatCover=="Open") #subset data by Open habitat cover
ResultsToPlotMEAN[1:3,2]<-apply(Closed[,8:10], 2, mean) #calculate mean logged RAIs of each species in Closed habitats
ResultsToPlotMEAN[4:6,2]<-apply(Open[,8:10], 2, mean) #calculate mean logged RAIs of each species in Open habitats

#CALCULATE CONFIDENCE INTERVALS
for(i in 8:10){ #for each species (taken from columns with logged RAIs)
  c<-Closed[,i] #extract logged species RAIs from Closed Habitat Sites only
  o<-Open[,i] #extract logged species RAIs from Open Habitat Sites only
  standarderrorC<-sd(c)/sqrt(length(c)) #calculate standard error (standard deviation/square root of sample size)
  standarderrorO<-sd(o)/sqrt(length(o)) #for both closed and open habitat sites
  #For each species in closed and open habitats, 
  #calculate confidence intervals by adding and subtracting standard errors from the mean values
  #Save results in results matrix
  ResultsToPlotMEAN[i-7, 3]<-as.numeric(ResultsToPlotMEAN[i-7,2])-standarderrorC*1.96 
  ResultsToPlotMEAN[i-7, 4]<-as.numeric(ResultsToPlotMEAN[i-7,2])+standarderrorC*1.96
  ResultsToPlotMEAN[i-4, 3]<-as.numeric(ResultsToPlotMEAN[i-4,2])-standarderrorO*1.96
  ResultsToPlotMEAN[i-4, 4]<-as.numeric(ResultsToPlotMEAN[i-4,2])+standarderrorO*1.96
}

setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/5) Habitat Preference") #set working directory to save data
write.csv(ResultsToPlotMEAN, "Habitat Cover-MEAN Logged RAIS.csv") #save completed results matrix as CSV file

#PLOT GRAPHS
#MEANS
#Convert data from data class/type to another, so graphing function can read them 
ResultsToPlotMEAN<-as.data.frame(ResultsToPlotMEAN) #convert results matrix to a dataframe
ResultsToPlotMEAN$MeanOfLog<-as.numeric(as.character(ResultsToPlotMEAN$MeanOfLog)) #convert Mean values from factors to characters
ResultsToPlotMEAN$LowerCI<-as.numeric(as.character(ResultsToPlotMEAN$LowerCI)) #convert Lower confidence interval values from factors to characters
ResultsToPlotMEAN$UpperCI<-as.numeric(as.character(ResultsToPlotMEAN$UpperCI)) #convert Upper confidence interval values from factors to characters

#BARCHART OF MEAN LOGGED RAIS OF FOXES AND HEDGEHOGS IN CLOSED AND OPEN HABITATS
ggplot(ResultsToPlotMEAN[-c(1,4),], aes(Species, MeanOfLog, fill=HabitatCover)) + geom_bar(stat="identity", position = "dodge", color="black") +labs(title="Species RAIs by Habitat Cover", y="Mean of Logged RAI")+theme(axis.text.x = element_text(vjust=-0.05,angle = 90))+ scale_fill_discrete(name = "Habitat Cover")+geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI), width=.2,position=position_dodge(.9))
ggsave("Mean RAI by Habitat Cover-Mammal Logged RAI.jpeg") #save graph as JPEG. ensure graphics windows is adjusted to make figure clear

#BARCHART OF MEAN LOGGED RAIS OF HUMANS IN CLOSED AND OPEN HABITATS
ggplot(ResultsToPlotMEAN[c(1,4),], aes(Species, MeanOfLog, fill=HabitatCover)) + geom_bar(stat="identity", position = "dodge", color="black") +labs(title="Human RAIs by Habitat Cover", y="Mean RAI")+theme(axis.text.x = element_text(vjust=-0.05,angle = 90))+ scale_fill_discrete(name = "Habitat Cover")+geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI), width=.2,position=position_dodge(.9))
ggsave("Mean RAI by Habitat Cover-Humans Logged RAI.jpeg") #save graph as JPEG. ensure graphics windows is adjusted to make figure clear


###### GRAPHS WITH HABITAT TYPES ########

#### PLOT BOXPLOTS ####
#Saved graphs by exporting with a height of 400 as ggsave JPEGs are a bit cramped
#Separate graphs for each species to make graphs clearer

gatheredDATA_Fox<-subset(gatheredDATA, Species =="FoxLog") #filter narrowed dataset, made with 'gather' function earlier, by foxes
#Plot boxplot of logged fox RAIs
ggplot(gatheredDATA_Fox, aes(x=HabitatType, y=LogRAI)) + geom_boxplot(fill="orchid3")+labs(title="Fox Logged RAIs by Habitat Type", y="Logged Relative Abundance Index")+theme(axis.text.x = element_text(vjust=-0.05,angle = 90), plot.title=element_text(hjust = 0.5))
ggsave("Boxplot-Fox RAIs-Habitat Type.jpeg") #save graph as JPEG. ensure graphics windows is adjusted to make figure clear

gatheredDATA_Hedgehog<-subset(gatheredDATA, Species =="HedgehogLog") #filter narrowed dataset, made with 'gather' function earlier, by hedgehogs
#Plot boxplot of logged hedgehog RAIs
ggplot(gatheredDATA_Hedgehog, aes(x=HabitatType, y=LogRAI)) + geom_boxplot(fill="peru")+labs(title="Hedgehog Logged RAIs by Habitat Type", y="Logged Relative Abundance Index")+theme(axis.text.x = element_text(vjust=-0.05,angle = 90), plot.title=element_text(hjust = 0.5))
ggsave("Boxplot-Hedgehog RAIs-Habitat Type.jpeg") #save graph as JPEG. ensure graphics windows is adjusted to make figure clear

#Plot boxplot of logged human RAIs
ggplot(gatheredDATA_human, aes(x=HabitatType, y=LogRAI)) + geom_boxplot(fill="lightgoldenrod")+labs(title="Human Logged RAIs by Habitat Type", y="Logged Relative Abundance Index")+theme(axis.text.x = element_text(vjust=-0.05,angle = 90), plot.title=element_text(hjust = 0.5))
ggsave("Boxplot-Human RAIs-Habitat Type.jpeg") #save graph as JPEG. ensure graphics windows is adjusted to make figure clear


#### PLOT BARCHARTS WITH MEANS ####
#Not presented in write up, but an alternative way to view data
#CONSTRUCT TABLE FOR MEAN VALUES AND CONFIDENCE INTERVALS TO PLOT
HabType<-sort(unique(Habitats$HabitatType)) #vector of habitat type in alphabetical order
ResultsToPlotMEANHT<-as.data.frame(matrix(0, ncol=5, nrow=18)) #empty results  matrix
colnames(ResultsToPlotMEANHT)<-c("Species", "MeanLog", "LowerCI", "UpperCI", "HabitatType") #add column names to results matrix
ResultsToPlotMEANHT[,1]<-rep(c("Human", "Fox", "Hedgehog"), rep=6) #add species names to results matrix
ResultsToPlotMEANHT[,5]<-rep(HabType, each=3) #add habitat type names to results matrix

##CALCULATE MEAN LOGGED RAIS
for(i in 1:6){ #for each habitat type
  HT<-subset(RAIResults, HabitatType==HabType[i]) #filter results by that habitat type
  rows<-c(1,2,3)+3*(i-1) #define the correct row in the results matrix to put the newly calculated mean 
  ResultsToPlotMEANHT[rows,2]<-apply(HT[,8:10], 2, mean) #calculate mean of logged RAI for each species and add to the results matrix 
}

#CALCULATE CONFIDENCE INTERVALS
for(i in 1:6){ #for each habitat type
  HT<-subset(RAIResults, HabitatType==HabType[i]) #filter results by that habitat type
  for(j in 8:10){ #for each species (take from column with logged RAIs)
    d<-HT[,j] #select that species from results dataset filtered by the habitat type
    standarderror<-sd(d)/sqrt(length(d)) #calculate standard error (standard deviation/square root of sample size)
    rows<-(j-7)+3*(i-1) #define the correct row in the results matrix to put the newly calculated confidence intervals
    ResultsToPlotMEANHT[rows, 3]<-ResultsToPlotMEANHT[rows,2]-(standarderror*1.96) #calculate the lower confidence interval by subtracting standard error from the mean, and add to the results matrix
    ResultsToPlotMEANHT[rows, 4]<-ResultsToPlotMEANHT[rows,2]+(standarderror*1.96) #calculate the upper confidence interval by adding standard error to the mean, and add to the results matrix
  }
}

setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/5) Habitat Preference") #set working directory to save data
write.csv(ResultsToPlotMEANHT, "Habitat Type-MEAN Logged RAIS.csv") #save completed results matrix as CSV file

#PLOT GRAPHS

#WITH FOXES
Fox<-subset(ResultsToPlotMEANHT, Species =="Fox") #filter mean results by foxes
#BARCHART OF MEAN LOGGED RAIS OF FOXES BY HABITAT TYPE
ggplot(Fox, aes(HabitatType, MeanLog)) + geom_bar(stat="identity", position = "dodge", color="black", fill="orchid3") +labs(title="Fox RAIs by Habitat Type", y="Mean of Logged RAI")+theme(axis.text.x = element_text(vjust=-0.05,angle = 90), plot.title=element_text(hjust = 0.5))+geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI), width=.2,position=position_dodge(.9))
ggsave("Mean Sightings by Habitat Type-Fox.jpeg") #save graph as JPEG. ensure graphics windows is adjusted to make figure clear

#WITH HEDGEHOGS
Hedgehog<-subset(ResultsToPlotMEANHT, Species =="Hedgehog") #filter mean results by hedgehogs
#BARCHART OF MEAN LOGGED RAIS OF HEDGEHOGS BY HABITAT TYPE
ggplot(Hedgehog, aes(HabitatType, MeanLog)) + geom_bar(stat="identity", position = "dodge", color="black", fill="peru") +labs(title="Hedgehog RAIs by Habitat Type", y="Mean of Logged RAI")+theme(axis.text.x = element_text(vjust=-0.05,angle = 90), plot.title=element_text(hjust = 0.5))+geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI), width=.2,position=position_dodge(.9))
ggsave("Mean Sightings by Habitat Type-Hedgehog.jpeg") #save graph as JPEG. ensure graphics windows is adjusted to make figure clear

#WITH HUMANS
Human<-subset(ResultsToPlotMEANHT, Species =="Human") #filter mean results by humans
#BARCHART OF MEAN LOGGED RAIS OF HUMANS BY HABITAT TYPE
ggplot(Human, aes(HabitatType, MeanLog)) + geom_bar(stat="identity", position = "dodge", color="black", fill="lightgoldenrod") +labs(title="Human RAIs by Habitat Type", y="Mean of Logged RAI")+theme(axis.text.x = element_text(vjust=-0.05,angle = 90), plot.title=element_text(hjust = 0.5))+geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI), width=.2,position=position_dodge(.9))
ggsave("Mean Sightings by Habitat Type-Humans.jpeg") #save graph as JPEG. ensure graphics windows is adjusted to make figure clear





