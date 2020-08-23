#NB THE ORIGINAL DATASETS DO NOT HAVE PERMISSION TO BE SHARED
#THESE ORIGINAL DATASETS ARE NEEDED FOR THIS CODE TO RUN

#COMBINED SCRIPT TO DETERMINE THE:
#1) SURVEY LENGTH AND DATES OF EACH CAMERA SITE (LINES 8 to 41)
#2) TOTAL NUMBER OF SPECIES IMAGES DETECTED ACROSS ALL SURVEYS AND THE NUMBER OF CAMERA SITES OBSERVED AT(LINE 42 ONWARDS)
#CALCULATED FOR FOXES AND HEDGEHOGS BETWEEN 6PM AND 8AM, AND FOR FULL DAY FOR HUMANS

rm(list=ls())

##### DETERMINE SURVEY LENGTH AT EACH CAMERA SITE #####
#START AND END DATE FOR EACH CAMERA SITE
setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/Data") #set working directory
data<-read.csv("ExifPro-FINAL DATASET.csv", header=TRUE) #FULL, ORIGINAL EXIFPRO DATASET AT 136 SITES

#REMOVE SITES THAT HAVE POOR PLACEMENTS AND INCONSISTENT TAGGING AS DONE FOR THE ANALYSES
SitesToRemove<-c(15,26,54,61,62,63,66,67,69,70,71,72,75,76,77,78,86,108,117,118,146,147) #Vector of camera sites to remove
length(SitesToRemove) #22 sites to remove in total
PlaceIDFULL<-sort(unique(data$TagPlaceID)) #vector of all the camera sites with images
length(PlaceIDFULL) #136 sites with images in total
SitesToKeepAreFalse<-PlaceIDFULL %in% SitesToRemove #create a vector with sites to keep marked as FALSE and marked to remove are TRUE
SitesToKeep<-PlaceIDFULL[SitesToKeepAreFalse==FALSE] #Make a vector of sites to keep in for the analyses
length(SitesToKeep) #114 sites remaining once the 22 sites are removed
data<-subset(data, TagPlaceID %in% SitesToKeep) #subset the ExifPro data by only the sites to be kept in the analyses

#DETERMINE LENGTH OF SURVEY AT EACH SIGHT
PlaceID<-sort(unique(data$TagPlaceID)) #Vector of sites used in the analyses
data$Date<-as.character(data$Date) #convert Date from a factor to a character
Days<-rep(0, length(PlaceID)) #create an empty results vector for "Survey Lenths" to be entered into
StartDate<-rep(0, length(PlaceID)) #create an empty results vector for "Start Date of Survey" to be entered into
EndDate<-rep(0, length(PlaceID)) #create an empty results vector for "End Date of Survey" to be entered into
for(i in 1:length(PlaceID)){ #for each camera site
  d<-subset(data, TagPlaceID==PlaceID[i]) #subset data by each camera site
  Days[i]<-length(unique(d$Date)) #determine survey length at that site and save in results vector
  Dates<-unique(d$Date) #vector of unique dates kept in chronological order
  StartDate[i]<-Dates[1] #determine start date of survey at that camera site and save in results vector
  EndDate[i]<-Dates[length(Dates)]} #determine end date of survey at that camera site and save in results vector
range(Days) #Survey lengths range between 4 to 43 days
mean(Days) #Mean survey length is 15.08772 days
SURVEY<-as.data.frame(cbind(PlaceID, Days, StartDate, EndDate)) #combine survey lengths and dates into one dataframe

setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/2) Miscellaneous") #set working directory to save results
write.csv(SURVEY, "Survey Lengths.csv") #save csv file containing survey lengths, and start and end dates of each survey


##### DETERMINE THE NUMBER OF SIGHTINGS OF EACH SPECIES AND THE NUMBER OF CAMERA SITES SEEN AT #####
#LOAD AND CHECK data
setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/Data") #set working directory
exifNight<-read.csv("exifNight.csv", header=TRUE) #ExifPro data for mammal sightings
View(exifNight)
str(exifNight)

mdNIGHT<-read.csv("mdNIGHT.csv", header=TRUE) #MegaDetector data between 6pm and 8am-to determine number of human sightings at night
View(mdNIGHT)
str(mdNIGHT)

mdFULL<-read.csv("mdFULL.csv", header=TRUE) #MegaDetector data between 6pm and 8am-to determine number of human sightings during the day
View(mdFULL)
str(mdFULL)

#FUNCTION TO DETERMINE THE NUMBER OF IMAGES OF EACH SPECIES ACROSS ALL SITES
SIGHTINGS<-function(exifNight, mdNIGHT, mdFULL){
  
  #CONVERT RELEVANT COLUMNS FROM FACTORS TO CHARACTERS
  exifNight$TagSpecies1<-as.character(exifNight$TagSpecies1)
  exifNight$TagSpecies2<-as.character(exifNight$TagSpecies2)
  exifNight$TagSpecies3<-as.character(exifNight$TagSpecies3)
  exifNight$Date<-as.character(exifNight$Date)
  mdNIGHT$category<-as.character(mdNIGHT$category)
  mdNIGHT$Date<-as.character(mdNIGHT$Date)
  
  #MAMMALS TO EXAMINE
  Species<-c("Fox", "Hedgehog")
  
  #CREATE RESULTS MATRIX
  resultsN<-matrix(0, ncol=2, nrow=4) #empty results matrix
  rownames(resultsN)<-c( Species, "Human (Night)", "Human (Day)") #row names for results matrix
  colnames(resultsN)<-c("TotalImages", "Sites Seen At") #column names for results matrix
  
  #COUNT NUMBER OF IMAGES OF EACH SPECIES AND NUMBER OF CAMERA SITES SEEN AT
  for(s in 1:2){ #for each mammal
    d<-subset(exifNight, subset = TagSpecies1== Species[s] | TagSpecies2 == Species[s]  | TagSpecies3 == Species[s]) #filter by specified mammal
    
    #APPLY ONE MINUTE INTERVAL TO CORRECT FOR CAMERA DIFFERENCES
    if(length(d$File)>1){ #if there is at least two images to compare
      
      record<-1 #for the initial comparison, use the first image as the reference image
      repeat{
        i<-record[length(record)] #for the last image marked as being greater than one minute from the last reference image, mark it as the new reference image
        for(j in (i+1):length(d$File)){ #compare each successive image's timestamp to the reference image's timestamp
          if(d$Date[i]==d$Date[j]){ #if images were taken on the same date (images in chronological order to enable these comparisons)
            if(d$TimeDecimal[j]>=(d$TimeDecimal[i]+1/1440)){#and if time taken was one minute or more after the reference image
              record<-c(record,j) #record image as the new reference image and stop comparisons with the current reference image
              break}}else{ #if images were taken on different days
                
                if(d$TimeDecimal[i]>=1-1/1440 & d$TimeDecimal[j]<=1/1440){ #if image on the first day was taken later than 11:59pm and if the image taken on the second day was taken before 0:01am (assuming no day without at least one image)
                  if(((1-d$TimeDecimal[i])+d$TimeDecimal[j])>=1/1440){ #if the time difference between the images being taken is greater than one minute 
                    record<-c(record,j) #make the next day image the new reference image, and stop comparisons with the current reference image
                    break}}else{
                      record<-c(record,j) #else time difference is larger than one minute between images taken on different days, so next day image becomes the new reference image
                      break}
              }
        }
        if(j==length(d$File)){break} #once comparisons have made until the very last image, stop the entire procedure
      }
    }else{ if(length(d$File)==1){record<-1}} #if only one image of that species, only this image is counted for mammal number
    
    OM<-d[record,] #new dataset contains only the images of a particular mammal with at least one minute between them
    
    #TOTAL NUMBER OF IMAGES
    resultsN[s,1]<-length(OM$File) #number of images remaining is the number of images of that specified mammal
    resultsN[s,2]<-length(unique(OM$TagPlaceID)) #count the number of unique camera sites the specified mammal has been sighted at
  }
  
  #TOTAL NUMBER OF HUMAN IMAGES DURING THE NIGHT (6PM TO 8AM)
  d<-subset(mdNIGHT, category=="person") #filter the MegaDetector, with images between 6pm and 8am, for human tags
  
  #APPLY ONE MINUTE INTERVAL TO CORRECT FOR CAMERA DIFFERENCES
  #FOR ANNOTATIONS SEE ABOVE WITH MAMMAL TAGS
  if(length(d$File)>1){
    record<-1
    repeat{
      i<-record[length(record)]
      for(j in (i+1):length(d$File)){
        if(d$Date[i]==d$Date[j]){
          if(d$TimeDecimal[j]>=(d$TimeDecimal[i]+1/1440)){ 
            record<-c(record,j)
            break}}else{ 
              if(d$TimeDecimal[i]>=1-1/1440 & d$TimeDecimal[j]<=1/1440){ 
                if(((1-d$TimeDecimal[i])+d$TimeDecimal[j])>=1/1440){ 
                  record<-c(record,j)
                  break}}else{
                    record<-c(record,j)
                    break}
            }
      }
      if(j==length(d$File)){break}
    }
  }else{ if(length(d$File)==1){record<-1}}
  
  OM<-d[record,]
  
  #TOTAL NUMBER OF IMAGES
  resultsN[3,1]<-length(OM$File)
  resultsN[3,2]<-length(unique(OM$TagPlaceID))
  
  #NUMBER OF HUMAN IMAGES DURING THE DAY (8AM TO 6PM)
  d<-subset(mdFULL, category=="person") #filter the MegaDetector for human tags
  d<-subset(d, TimeDecimal>8/24 & TimeDecimal <18/24) #filter the MegaDetector for hours during the day
  
  #APPLY ONE MINUTE INTERVAL TO CORRECT FOR CAMERA DIFFERENCES
  #FOR ANNOTATIONS SEE ABOVE WITH MAMMAL TAGS
  if(length(d$File)>1){
    record<-1
    repeat{
      i<-record[length(record)]
      for(j in (i+1):length(d$File)){
        if(d$Date[i]==d$Date[j]){
          if(d$TimeDecimal[j]>=(d$TimeDecimal[i]+1/1440)){ 
            record<-c(record,j)
            break}}else{ 
              if(d$TimeDecimal[i]>=1-1/1440 & d$TimeDecimal[j]<=1/1440){ 
                if(((1-d$TimeDecimal[i])+d$TimeDecimal[j])>=1/1440){ 
                  record<-c(record,j)
                  break}}else{
                    record<-c(record,j)
                    break}
            }
      }
      if(j==length(d$File)){break}
    }
  }else{ if(length(d$File)==1){record<-1}}
  
  OM<-d[record,]
  
  #TOTAL NUMBER OF IMAGES
  resultsN[4,1]<-length(OM$File)
  resultsN[4,2]<-length(unique(OM$TagPlaceID))
  
  setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/2) Miscellaneous") #set working directory to save results matrix
  write.csv(resultsN, "Species Sightings.csv") #save CSV of completed results matrix
}

SIGHTINGS(exifNight, mdNIGHT, mdFULL) #run function
#HUMANS PRESENT AT 110 SITES ACROSS THE FULL 24 HOURS. THEY WERE NOT DETECTED AT SITES 44, 49, 90 AND 116

##### PLOT GRAPH OF NUMBER OF IMAGES #####
#Not used in the write up, but useful in case graphs are preferred
#INSTALL AND LOAD GGPLOT FOR PLOTTING GRAPHS
install.packages("ggplot2") 
library(ggplot2)

setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/2) Miscellaneous") #set working directory to retrieve results file
sightings<-read.csv("Species Sightings.csv", header=TRUE) #load results file of species sightings

#PLOT GRAPHS
#PLOT NUMBER OF IMAGES FOR FOXES AND HEDGEHOGS AS A BAR CHART
ggplot(sightings[1:2,], aes(X,TotalImages)) + geom_bar(stat="identity", position = "dodge", fill="blue")+labs(title="Total Species Sightings (Night)", x="Species", y="Total Images")+theme(axis.text.x = element_text(vjust=-0.05,angle = 90), plot.title=element_text(hjust = 0.5))
#PLOT NUMBER OF IMAGES FOR ALL SPECIES, INCLUDING DAY AND NIGHT IMAGES OF HUMANS
ggplot(sightings[1:4,], aes(X,TotalImages)) + geom_bar(stat="identity", position = "dodge", fill="coral")+labs(title="Total Species Sightings", x=" Species")+theme(axis.text.x = element_text(vjust=-0.05,angle = 90), plot.title=element_text(hjust = 0.5))
#PLOT NUMBER OF IMAGES OF HUMANS DURING THE DAY AND NIGHT
ggplot(sightings[3:4,], aes(X,TotalImages)) + geom_bar(stat="identity", position = "dodge", fill="lightgoldenrod")+labs(title="Total Human Sightings (Night/Day)", y="Total Images", x=" ")+theme(plot.title=element_text(hjust = 0.5))

#PLOTS NOT SAVED AS NOT NEEDED FOR WRITE UP.


