#NB THE ORIGINAL DATASETS DO NOT HAVE PERMISSION TO BE SHARED
#THESE ORIGINAL DATASETS ARE NEEDED FOR THIS CODE TO RUN

#DETERMINE WHICH CONFIDENCE VALUE TO FILTER THE MEGADETECTOR HUMAN TAGS BY
#AT LOWER CONFIDENCE VALUES IT IS LIKELY TREES AND ANIMALS MAY HAVE BEEN TAGGED AS
#HUMANS, REDUCING THE ACCURACY OF THE RESULTS. THEREFORE NEED TO FILTER OUT THESE 
#INACCURATE TAGS.

#TO ASSESS THE ACCURACY OF THE MEGADETECTOR, I COMPARED THE HUMAN SIGHTINGS DETECTED 
#MANUALLY (BY DR JEFF WAAGE) AND AUTOMATICALLY BY THE MEGADETECTOR 
#COMPARISONS OCCURRED ACROSS DIFFERENT CONFIDENCE VALUES OF THE MEGADETECTOR TAGS 
#TO DETERMINE WHICH IS THE BEST CONFIDENCE VALUE TO USE
#N.B. MegaDetector Dataset has already been filtered by a confidence value of 0.2
#as the initial output file would have been too large too save

rm(list=ls())

####### LOAD AND CHECK DATA #######
setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/Data") #set working directory to retrieve data
data<-read.csv("J.Waage Data-27 Sites.csv", header=TRUE) #Human Tags by Dr Jeff Waage DETECTION
View(data)
str(data)

MD<-read.csv("MegaDetector-FINAL DATASET.csv" , header=TRUE) #Human Tags by the MegaDetector-Full Original Dataset
str(MD)
View(MD)

####### DATA PREPARATION #####
#TRUNCATE SET UP AND TAKE DOWN TIMES FROM MEGADETECTOR DATA
#AS DR JEFF WAAGE IGNORED THESE IN HIS HUMAN TAGGING
#I.E. REMOVE FIRST AND LAST 15 MINUTES FROM THE SURVEYS OF EACH SITE
PlaceIDMD<-sort(unique(MD$TagPlaceID)) #vector of camera sites scanned by megadetector
MDTruncated<-as.data.frame(matrix(0, nrow=2, ncol=7)) #starting point for the new truncated dataset i.e. blank starting rows where truncated surveys by site can be added 
colnames(MDTruncated)<-colnames(MD) #Add column names to the new truncated dataset

for(i in 1:length(PlaceIDMD)){ #for each camera site
  d<-subset(MD, TagPlaceID==PlaceIDMD[i]) #select by that site
  for(j in 2:length(d$File)){ #for each successive image after the first image
    if(d$TimeDecimal[j]>=(d$TimeDecimal[1]+15/1440)){ #As data is in chronological order, and camera setups/takedowns occur during the day (so will not be affected by the return to TimeDecimal==0 at midnight), compare succesive image time stamps until one is greater than 15 minutes after the first image's timestamp. 
      newstart<-j #new starting image for the truncated dataset at that site is at jth row (jth image is the first one that is at least 15 minutes after the timestamp of the first image)
      break}} #once new starting image found, cease the search
  
  for(k in 1:length(d$File)){ #for each file before the final one
    if(identical(d$Date[length(d$File)],d$Date[length(d$File)-k])){ #if the date is the same for the images being compared
      if(d$TimeDecimal[length(d$File)-k]<=(d$TimeDecimal[length(d$File)]-15/1440)){ #compare successive images (from newest to oldest) to the final image, to determine the first image that was taken at least 15 minutes before the last image 
        newend<-length(d$File)-k #New final image for the truncated dataset at that site is at kth row (kth image is the first one that is at least 15 minutes before the timestamp of the last image)
        break}}else{ #once a new final image is found, cease the search
          newend<-length(d$File) #if on separate days, files will be greater than 15 minutes apart as set up or take downs never occur at night. So that image on a different day becomes the new final image.
        }
  }
  
  New<-d[newstart:newend,] #truncate the data set at that site by the new starting and final tags found
  MDTruncated<-as.data.frame(rbind(MDTruncated, New)) #bind the truncated dataset at this sites to ones at other sites
}
MDTruncated<-MDTruncated[-c(1,2),] #remove blank starting rows
rownames(MDTruncated)<-1:length(MDTruncated$File) #rename the rows based on the number that are now present
View(MDTruncated) #check final result has formed correctly

#SELECT ONLY THE SITES THAT DR JEFF WAAGE TAGGED
PlaceIDJ<-sort(unique(data$TagPlaceID)) #vector of sites tagged by Dr Jeff Waage
length(PlaceIDJ) #27 sites tagged by Dr Jeff Waage
PlaceIDMD<-unique(MDTruncated$TagPlaceID) #vector of sites scanned by MegaDetector
length(PlaceIDMD) #136 sites tagged by the MegaDetector

MDTJ<-subset(MDTruncated, TagPlaceID %in% PlaceIDJ) #filter MegaDetector data to include Dr Jeff Waage's sites only
length(sort(unique(MDTJ$TagPlaceID))) #26 sites shared between the MegaDetector and Dr Jeff Waage's dataset
PlaceIDJ #view Dr Jeff Waage's sites
sort(unique(MDTJ$TagPlaceID)) #view MegaDetector's sites
#From examining camera sites in Dr Waage's and the MegaDetector's dataset
#MegaDetector did not scan site 129, leaving 26 sites in total to compare
#In addition, sites 43 and 87 were incompletely tagged by Dr Jeff Waage due to the sheer number of humans present
#Sites 43 and 87 removed for the final analysis, but their effect can be seen in the first analysis below
#IN TOTAL 24 SITES USED IN THE ANALYSIS. 

###### EXAMINE RELATIONSHIPS BETWEEN MANUAL AND AUTOMATIC TAGS ######
#### USES MEGADETECTOR WITH A 0.5 CONFIDENCE VALUE (ARBITRARY START POINT) ####
MDTJ0.5<-subset(MDTJ, conf>=0.5) #Filter for a minimum confidence value fo 0.5
PlaceID<-sort(unique(MDTJ0.5$TagPlaceID)) #vector of shared sites 
SightingsJ<-rep(0, length(PlaceID)) #empty results vector for human tags by Dr Jeff Waage at each site
SightingsM<-rep(0, length(PlaceID)) #empty results vector for human tags by MegaDetector at each site

for(i in 1:length(PlaceID)){ #for each camera site
  j<-subset(data, TagPlaceID==PlaceID[i]) #filter by that site
  j1<-subset(j, TagSpecies1=="Human") #subset by Human Tags
  j2<-subset(j, TagSpecies2=="Human") 
  j3<-subset(j, TagSpecies3=="Human")
  SightingsJ[i]<-sum(length(j1$File), length(j2$File), length(j3$File)) #Sum number of number of human tags tagged manually
  m<-subset(MDTJ0.5, TagPlaceID==PlaceID[i] & category=="person") #filter MegaDetector by that site and for human tags
  SightingsM[i]<-length(m$File) #Count number of humans tags tagged by the MegaDetector
  }
sightingsRESULT0.5<-as.data.frame(cbind(PlaceID, SightingsJ, SightingsM)) #combine sightings into a results data frame
View(sightingsRESULT0.5) #view results data frame. 3 sites with a few humans detected when Dr Jeff Waage tagged none (likely an animal mistaken for a human)

#PLOT GRAPHS AND STATISTICAL ANALYSIS
#ANALYSIS: SIMPLE LINEAR REGRESSION
plot(sightingsRESULT0.5$SightingsM~sightingsRESULT0.5$SightingsJ, main="Jeff vs MD Human Sightings (0.5)", xlab="Jeff's Sightings", ylab="MD Sightings", pch=19, col="lightblue",cex=2) #plot manual human tags against automatic human tags
text(sightingsRESULT0.5$SightingsM~sightingsRESULT0.5$SightingsJ, labels=sightingsRESULT0.5$PlaceID,cex=0.9, font=2) #label each point by its camera site
lines(x = c(0,15000), y = c(0,15000), col="red") #add y=x line to more clearly see if there is a 1:1 relationship 

model<-lm(sightingsRESULT0.5$SightingsM~sightingsRESULT0.5$SightingsJ)
summary(model) #significant but slope is 2.2106 (i.e. the MegaDetector is detecting more than twice the humans there actually are)

#THEREFORE ESSENTIAL TO REMOVE SITES 43 AND 87 WHICH WERE INCOMPLETELY MANUALLY TAGGED
sightingsRESULT1<-sightingsRESULT0.5[-c(6,17),]

#GRAPHS AND ANALYSES REDONE WITHOUT THESE SITES
#FOR ANNOTATIONS OF PLOTTING SEE ABOVE (0.5 WITH ALL SITES)
plot(sightingsRESULT1$SightingsM~sightingsRESULT1$SightingsJ, main="Jeff vs MD Human Sightings (no 87/43 0.5)", xlab="Jeff's Sightings", ylab="MD Sightings", pch=19, col="lightblue",cex=2)
text(sightingsRESULT1$SightingsM~sightingsRESULT1$SightingsJ, labels=sightingsRESULT1$PlaceID,cex=0.9, font=2)
lines(x = c(0,15000), y = c(0,15000), col="red") 

model1<-lm(sightingsRESULT1$SightingsM~sightingsRESULT1$SightingsJ)
summary(model1) #slope is 1.04796, showing close to a 1 to 1 relationship. R squared=0.95

#REMOVE SITES 43 AND 87 FOR THE REMAINDER OF ANALYSES
j<-subset(j,TagPlaceID != 43 & TagPlaceID != 87) #remove sites from Dr Jeff Waage's dataset
MDTJ<-subset(MDTJ, TagPlaceID != 43 & TagPlaceID != 87) #remove sites from MegaDetector dataset

##### REPEAT PROCEDURE AGAIN WITH HIGHER CONFIDENCE VALUES #####
#### CONFIDENCE VALUE-0.6 ####
#N.B. For annotations see code above, as procedure is repeated
MDTJ0.6<-subset(MDTJ, conf>=0.6) #filter by confidence value >=0.6

#EXAMINE RELATIONSHIPS BETWEEN SIGHTINGS
#USES MEGADETECTOR FILTERD BY A 0.6 CONFIDENCE VALUE
PlaceID<-sort(unique(MDTJ0.6$TagPlaceID)) 
SightingsJ<-rep(0, length(PlaceID)) 
SightingsM<-rep(0, length(PlaceID))

for(i in 1:length(PlaceID)){ 
  j<-subset(data, TagPlaceID==PlaceID[i]) 
  j1<-subset(j, TagSpecies1=="Human") 
  j2<-subset(j, TagSpecies2=="Human")
  j3<-subset(j, TagSpecies3=="Human")
  SightingsJ[i]<-sum(length(j1$File), length(j2$File), length(j3$File)) 
  m<-subset(MDTJ0.6, TagPlaceID==PlaceID[i] & category=="person") 
  SightingsM[i]<-length(m$File) 
}
sightingsRESULT0.6<-as.data.frame(cbind(PlaceID, SightingsJ, SightingsM)) 
View(sightingsRESULT0.6) #3 sites with a few humans detected where Dr Jeff Waage tagged none (likely an animal mistaken for a human)

#GRAPHS AND STATISTICAL ANALYSIS
par(mar=c(4,4,1,1)) #adjust margins to ensure figure can be plotted in graphic windows
plot(sightingsRESULT0.6$SightingsM~sightingsRESULT0.6$SightingsJ, main="Jeff vs MD Human Sightings (0.6)", xlab="Jeff's Sightings", ylab="MD Sightings", pch=19, col="lightblue",cex=2)
text(sightingsRESULT0.6$SightingsM~sightingsRESULT0.6$SightingsJ, labels=sightingsRESULT0.6$PlaceID,cex=0.9, font=2)
lines(x = c(0,15000), y = c(0,15000), col="red") 

model2<-lm(sightingsRESULT0.6$SightingsM~sightingsRESULT0.6$SightingsJ)
summary(model2) #slope=1.01344 (almost 1:1). R Squared=0.96

#### CONFIDENCE VALUE-0.7 ####
MDTJ0.7<-subset(MDTJ, conf>=0.7) #filter by confidence value >=0.7

#EXAMINE RELATIONSHIPs BETWEEN SIGHTINGS
#USES MEGADETECTOR FILTERED BY A 0.7 CONFIDENCE VALUE
PlaceID<-sort(unique(MDTJ0.7$TagPlaceID)) 
SightingsJ<-rep(0, length(PlaceID)) 
SightingsM<-rep(0, length(PlaceID)) 

for(i in 1:length(PlaceID)){ 
  j<-subset(data, TagPlaceID==PlaceID[i]) 
  j1<-subset(j, TagSpecies1=="Human") 
  j2<-subset(j, TagSpecies2=="Human")
  j3<-subset(j, TagSpecies3=="Human")
  SightingsJ[i]<-sum(length(j1$File), length(j2$File), length(j3$File)) 
  m<-subset(MDTJ0.7, TagPlaceID==PlaceID[i] & category=="person") 
  SightingsM[i]<-length(m$File) 
}
sightingsRESULT0.7<-as.data.frame(cbind(PlaceID, SightingsJ, SightingsM)) 
View(sightingsRESULT0.7) #2 sites with a few humans detected where Dr Jeff Waage tagged none (likely an animal mistaken for a human)

#GRAPHS AND STATISTICAL ANALYSIS
plot(sightingsRESULT0.7$SightingsM~sightingsRESULT0.7$SightingsJ, main="Jeff vs MD Human Sightings (0.7)", xlab="Jeff's Sightings", ylab="MD Sightings", pch=19, col="lightblue",cex=2)
text(sightingsRESULT0.7$SightingsM~sightingsRESULT0.7$SightingsJ, labels=sightingsRESULT0.7$PlaceID,cex=0.9, font=2)
lines(x = c(0,15000), y = c(0,15000), col="red") 

model3<-lm(sightingsRESULT0.7$SightingsM~sightingsRESULT0.7$SightingsJ)
summary(model3) #slope=0.98553 (almost 1:1). R squared=0.96

#### CONFIDENCE VALUE-0.8 ####
MDTJ0.8<-subset(MDTJ, conf>=0.8) #filter by confidence value >=0.8

#EXAMINE RELATIONSHIPs BETWEEN SIGHTINGS
#USES MEGADETECTOR FILTERED BY A 0.8 CONFIDENCE VALUE
PlaceID<-sort(unique(MDTJ0.8$TagPlaceID)) 
SightingsJ<-rep(0, length(PlaceID)) 
SightingsM<-rep(0, length(PlaceID)) 

for(i in 1:length(PlaceID)){ 
  j<-subset(data, TagPlaceID==PlaceID[i]) 
  j1<-subset(j, TagSpecies1=="Human") 
  j2<-subset(j, TagSpecies2=="Human")
  j3<-subset(j, TagSpecies3=="Human")
  SightingsJ[i]<-sum(length(j1$File), length(j2$File), length(j3$File)) 
  m<-subset(MDTJ0.8, TagPlaceID==PlaceID[i] & category=="person") 
  SightingsM[i]<-length(m$File) 
}
sightingsRESULT0.8<-as.data.frame(cbind(PlaceID, SightingsJ, SightingsM))
View(sightingsRESULT0.8) #2 sites with humans detected when Dr Jeff Waage tagged none

#GRAPHS AND STATISTICAL ANALYSIS
plot(sightingsRESULT0.8$SightingsM~sightingsRESULT0.8$SightingsJ, main="Jeff vs MD Human Sightings (0.8)", xlab="Jeff's Sightings", ylab="MD Sightings", pch=19, col="lightblue",cex=2)
text(sightingsRESULT0.8$SightingsM~sightingsRESULT0.8$SightingsJ, labels=sightingsRESULT0.8$PlaceID,cex=0.9, font=2)
lines(x = c(0,15000), y = c(0,15000), col="red") 

model4<-lm(sightingsRESULT0.8$SightingsM~sightingsRESULT0.8$SightingsJ)
summary(model4) #slope=0.96, R squared 0.96

#### CONFIDENCE VALUE-0.9 ####
MDTJ0.9<-subset(MDTJ, conf>=0.9) #filter by confidence VALUE >=0.9

#EXAMINE RELATIONSHIPs BETWEEN SIGHTINGS
#USES MEGADETECTOR FILTERED BY A 0.9 CONFIDENCE VALUE
PlaceID<-sort(unique(MDTJ0.9$TagPlaceID)) 
SightingsJ<-rep(0, length(PlaceID)) 
SightingsM<-rep(0, length(PlaceID)) 

for(i in 1:length(PlaceID)){ 
  j<-subset(data, TagPlaceID==PlaceID[i]) 
  j1<-subset(j, TagSpecies1=="Human") 
  j2<-subset(j, TagSpecies2=="Human")
  j3<-subset(j, TagSpecies3=="Human")
  SightingsJ[i]<-sum(length(j1$File), length(j2$File), length(j3$File)) 
  m<-subset(MDTJ0.9, TagPlaceID==PlaceID[i] & category=="person") 
  SightingsM[i]<-length(m$File) 
}
sightingsRESULT0.9<-as.data.frame(cbind(PlaceID, SightingsJ, SightingsM)) 
View(sightingsRESULT0.9) #2 sites with humans detected where Dr Jeff Waage tagged none

#GRAPHS AND STATISTICAL ANALYSIS
plot(sightingsRESULT0.9$SightingsM~sightingsRESULT0.9$SightingsJ, main="Jeff vs MD Human Sightings (0.9)", xlab="Jeff's Sightings", ylab="MD Sightings", pch=19, col="lightblue",cex=2)
text(sightingsRESULT0.9$SightingsM~sightingsRESULT0.9$SightingsJ, labels=sightingsRESULT0.9$PlaceID,cex=0.9, font=2)
lines(x = c(0,15000), y = c(0,15000), col="red") 

model5<-lm(sightingsRESULT0.9$SightingsM~sightingsRESULT0.9$SightingsJ)
summary(model5) #slope=0.92, R squared 0.96

#### REFINING THE VALUES BY EXAMINING AT 0.01 VALUES ####
#FROM SLOPES, A CONFIDENCE VALUE OF BETWEEN 0.6 OR 0.7 WILL LEAD TO A CLOSE 
#TO A 1:1 RELATIONSHIP AS POSSIBLE

#### CONFIDENCE VALUE-0.63 ####
#REPEAT AGAIN WITH AN INTERMEDIATE CONFIDENCE VALUE
#ARBITRARY STARTING POINTS OF 0.63-TO REDUCE NUMBER OF MODELS RUN
MDTJ0.63<-subset(MDTJ, conf>=0.63) #filter by confidence value >=0.63

#EXAMINE RELATIONSHIPS BETWEEN SIGHTINGS
#USES MEGADETECTOR FILTERED BY A 0.63 CONFIDENCE VALUE
PlaceID<-sort(unique(MDTJ0.63$TagPlaceID)) 
SightingsJ<-rep(0, length(PlaceID)) 
SightingsM<-rep(0, length(PlaceID)) 

for(i in 1:length(PlaceID)){ 
  j<-subset(data, TagPlaceID==PlaceID[i]) 
  j1<-subset(j, TagSpecies1=="Human") 
  j2<-subset(j, TagSpecies2=="Human")
  j3<-subset(j, TagSpecies3=="Human")
  SightingsJ[i]<-sum(length(j1$File), length(j2$File), length(j3$File)) 
  m<-subset(MDTJ0.63, TagPlaceID==PlaceID[i] & category=="person") 
  SightingsM[i]<-length(m$File)
}
sightingsRESULT0.63<-as.data.frame(cbind(PlaceID, SightingsJ, SightingsM)) 
View(sightingsRESULT0.63) #3 sites with humans detected where Dr Jeff Waage tagged none

#GRAPHS AND STATISTICAL ANALYSIS
plot(sightingsRESULT0.63$SightingsM~sightingsRESULT0.63$SightingsJ, main="Jeff vs MD Human Sightings (0.63)", xlab="Jeff's Sightings", ylab="MD Sightings", pch=19, col="lightblue",cex=2)
text(sightingsRESULT0.63$SightingsM~sightingsRESULT0.63$SightingsJ, labels=sightingsRESULT0.63$PlaceID,cex=0.9, font=2)
lines(x = c(0,15000), y = c(0,15000), col="red") 

model6<-lm(sightingsRESULT0.63$SightingsM~sightingsRESULT0.63$SightingsJ)
summary(model6) #slope=1.00581, R squared 0.96

#### CONFIDENCE VALUE-0.64 ####
MDTJ0.64<-subset(MDTJ, conf>=0.64) #filter by confidence value >=0.64

#EXAMINE RELATIONSHIPS BETWEEN SIGHTINGS
#USES MEGADETECTOR FILTERED BY A 0.64 CONFIDENCE VALUE
PlaceID<-sort(unique(MDTJ0.64$TagPlaceID)) 
SightingsJ<-rep(0, length(PlaceID)) 
SightingsM<-rep(0, length(PlaceID)) 

for(i in 1:length(PlaceID)){ 
  j<-subset(data, TagPlaceID==PlaceID[i]) 
  j1<-subset(j, TagSpecies1=="Human") 
  j2<-subset(j, TagSpecies2=="Human")
  j3<-subset(j, TagSpecies3=="Human")
  SightingsJ[i]<-sum(length(j1$File), length(j2$File), length(j3$File)) 
  m<-subset(MDTJ0.64, TagPlaceID==PlaceID[i] & category=="person") 
  SightingsM[i]<-length(m$File) 
}
sightingsRESULT0.64<-as.data.frame(cbind(PlaceID, SightingsJ, SightingsM))
View(sightingsRESULT0.64) #3 sites with humans detected where Dr Jeff Waage tagged none

#GRAPHS AND STATISTICAL ANALYSIS
plot(sightingsRESULT0.64$SightingsM~sightingsRESULT0.64$SightingsJ, main="Jeff vs MD Human Sightings (0.64)", xlab="Jeff's Sightings", ylab="MD Sightings", pch=19, col="lightblue",cex=2)
text(sightingsRESULT0.64$SightingsM~sightingsRESULT0.64$SightingsJ, labels=sightingsRESULT0.64$PlaceID,cex=0.9, font=2)
lines(x = c(0,15000), y = c(0,15000), col="red") 

model7<-lm(sightingsRESULT0.64$SightingsM~sightingsRESULT0.64$SightingsJ)
summary(model7) #slope=1.0043 R squared 0.95

#### CONFIDENCE VALUE-0.65 ####
MDTJ0.65<-subset(MDTJ, conf>=0.65) #filter by confidence value >=0.65

#EXAMINE RELATIONSHIPS BETWEEN SIGHTINGS
#USES MEGADETECTOR FILTERED BY A 0.65 CONFIDENCE VALUE
PlaceID<-sort(unique(MDTJ0.65$TagPlaceID)) 
SightingsJ<-rep(0, length(PlaceID)) 
SightingsM<-rep(0, length(PlaceID)) 

for(i in 1:length(PlaceID)){ 
  j<-subset(data, TagPlaceID==PlaceID[i]) 
  j1<-subset(j, TagSpecies1=="Human") 
  j2<-subset(j, TagSpecies2=="Human")
  j3<-subset(j, TagSpecies3=="Human")
  SightingsJ[i]<-sum(length(j1$File), length(j2$File), length(j3$File)) 
  m<-subset(MDTJ0.65, TagPlaceID==PlaceID[i] & category=="person") 
  SightingsM[i]<-length(m$File) 
}
sightingsRESULT0.65<-as.data.frame(cbind(PlaceID, SightingsJ, SightingsM)) 
View(sightingsRESULT0.65) #2 sites with humans detected where Dr Jeff Waage tagged none

#GRAPHS AND STATISTICAL ANALYSIS
plot(sightingsRESULT0.65$SightingsM~sightingsRESULT0.65$SightingsJ, main="Jeff vs MD Human Sightings (0.65)", xlab="Jeff's Sightings", ylab="MD Sightings", pch=19, col="lightblue",cex=2)
text(sightingsRESULT0.65$SightingsM~sightingsRESULT0.65$SightingsJ, labels=sightingsRESULT0.65$PlaceID,cex=0.9, font=2)
lines(x = c(0,15000), y = c(0,15000), col="red") 

model8<-lm(sightingsRESULT0.65$SightingsM~sightingsRESULT0.65$SightingsJ)
summary(model8) #slope=0.99917, R squared 0.96

#### CONFIDENCE VALUE-0.66 ####
MDTJ0.66<-subset(MDTJ, conf>=0.66) #filter by confidence VALUE >=0.66

#EXAMINE RELATIONSHIPS BETWEEN SIGHTINGS
#USES MEGADETECTOR FILTERED BY A 0.66 CONFIDENCE VALUE
PlaceID<-sort(unique(MDTJ0.66$TagPlaceID)) 
SightingsJ<-rep(0, length(PlaceID)) 
SightingsM<-rep(0, length(PlaceID)) 

for(i in 1:length(PlaceID)){ 
  j<-subset(data, TagPlaceID==PlaceID[i]) 
  j1<-subset(j, TagSpecies1=="Human") 
  j2<-subset(j, TagSpecies2=="Human")
  j3<-subset(j, TagSpecies3=="Human")
  SightingsJ[i]<-sum(length(j1$File), length(j2$File), length(j3$File)) 
  m<-subset(MDTJ0.66, TagPlaceID==PlaceID[i] & category=="person") 
  SightingsM[i]<-length(m$File) 
}
sightingsRESULT0.66<-as.data.frame(cbind(PlaceID, SightingsJ, SightingsM)) 
View(sightingsRESULT0.66) #2 sites with humans detected when Dr Jeff Waage tagged none

#GRAPHS AND STATISTICAL ANALYSIS
plot(sightingsRESULT0.66$SightingsM~sightingsRESULT0.66$SightingsJ, main="Jeff vs MD Human Sightings (0.66)", xlab="Jeff's Sightings", ylab="MD Sightings", pch=19, col="lightblue",cex=2)
text(sightingsRESULT0.66$SightingsM~sightingsRESULT0.66$SightingsJ, labels=sightingsRESULT0.66$PlaceID,cex=0.9, font=2)
lines(x = c(0,15000), y = c(0,15000), col="red") 

model9<-lm(sightingsRESULT0.66$SightingsM~sightingsRESULT0.66$SightingsJ)
summary(model9) #slope=0.99665, R squared 0.95

#### CONFIDENCE VALUE-0.67 ####
MDTJ0.67<-subset(MDTJ, conf>=0.67) #filter by confidence VALUE >=0.67

#EXAMINE RELATIONSHIPS BETWEEN SIGHTINGS
#USES MEGADETECTOR FILTERED WITH A 0.67 CONFIDENCE VALUE
PlaceID<-sort(unique(MDTJ0.67$TagPlaceID)) 
SightingsJ<-rep(0, length(PlaceID)) 
SightingsM<-rep(0, length(PlaceID)) 

for(i in 1:length(PlaceID)){ 
  j<-subset(data, TagPlaceID==PlaceID[i]) 
  j1<-subset(j, TagSpecies1=="Human") 
  j2<-subset(j, TagSpecies2=="Human")
  j3<-subset(j, TagSpecies3=="Human")
  SightingsJ[i]<-sum(length(j1$File), length(j2$File), length(j3$File)) 
  m<-subset(MDTJ0.67, TagPlaceID==PlaceID[i] & category=="person") 
  SightingsM[i]<-length(m$File) 
}

sightingsRESULT0.67<-as.data.frame(cbind(PlaceID, SightingsJ, SightingsM)) 
View(sightingsRESULT0.67) #2 sites with humans detected where Dr Jeff Waage tagged none

#GRAPHS AND STATISTICAL ANALYSIS
plot(sightingsRESULT0.67$SightingsM~sightingsRESULT0.67$SightingsJ, main="Jeff vs MD Human Sightings (0.67)", xlab="Jeff's Sightings", ylab="MD Sightings", pch=19, col="lightblue",cex=2)
text(sightingsRESULT0.67$SightingsM~sightingsRESULT0.67$SightingsJ, labels=sightingsRESULT0.67$PlaceID,cex=0.9, font=2)
lines(x = c(0,15000), y = c(0,15000), col="red") 

model10<-lm(sightingsRESULT0.67$SightingsM~sightingsRESULT0.67$SightingsJ)
summary(model10) #slope=0.99299, R squared 0.96







