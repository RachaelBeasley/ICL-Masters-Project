#DATA PREPARATION FOR ANALYSES
#DATASETS: EXIF PRO AND MEGADETECTOR

#ACTIONS COMPLETED IN THIS R SCRIPT:
#REMOVE SITES WITH POOR PLACEMENT, IMAGE CAPTURE AND INCONSISTENT TAGGING
#REMOVES FIRST AND LAST 15 MINUTES OF EACH CAMERA SITE'S SURVEY-REMOVES SET UP AND TAKE DOWN TIMES
#MEGADETECTOR FILTERED BY A CONFIDENCE VALUE OF 0.65 (SEE 'Confidence Value Threshold' R SCRIPT)
#MEGADETECTOR WITH FULL 24 HOURS SAVED (mdFULL)
#MEGADETECTOR AND EXIFPRO FILTERED BY NIGHT TIME HOURS-6PM TO 8AM (exifNIGHT, mdNIGHT)
#N.B. data was still truncated despite using 'night hours' (6pm to 8am) as some cameras were taken down later than 6pm

rm(list=ls())

####### LOAD AND CHECK DATA  #######
setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/Data") #set working directory to retrieve original datasets
exif<-read.csv("ExifPro-FINAL DATASET.csv", header=TRUE) #ExifPro Data-Species Tags
View(exif) #visualise examine data
str(exif) #see structure of data

MD<-read.csv("MegaDetector-FINAL DATASET.csv" , header=TRUE) #MegaDetector Data-Human Tags
View(MD)
str(MD)

####### EXIF PRO PREPARATION #######

#REMOVE SITES WITH POOR PLACEMENT, IMAGE CAPTURE AND INCONSISTENT TAGGING
SitesToRemove<-c(15,26,54,61,62,63,66,67,69,70,71,72,75,76,77,78,86,108,117,118,146,147) #Vector of camera sites to remove
length(SitesToRemove) #22 Camera sites to be removed 
PlaceIDFULL<-sort(unique(exif$TagPlaceID)) #Get a list of all the sites where images have been taken
length(PlaceIDFULL) #136 Sites where images have been taken
SitesToKeepAreFalse<-PlaceIDFULL %in% SitesToRemove #Creates logical of sites to remove (TRUE) and sites to keep (false)
SitesToKeep<-PlaceIDFULL[SitesToKeepAreFalse==FALSE] #Creates vector of sites to keep for analyses
length(SitesToKeep) #114 Sites remaining for the analyses

exif<-subset(exif, TagPlaceID %in% SitesToKeep) #Filter ExifPro data so it contains only sites to be kept              
length(unique(exif$TagPlaceID)) #114 Sites to keep (Checks filtering has worked)

#TRUNCATE SETUP/TAKE DOWN TIME FROM EXIFPRO DATA 
#I.E. REMOVE FIRST AND LAST 15 MINUTES FROM THE SURVEYS OF EACH SITE
PlaceID<-sort(unique(exif$TagPlaceID)) #Vector of  sites from ExifPro data
exifTruncated<-as.data.frame(matrix(0, nrow=2, ncol=14)) #Starting point for the new truncated dataset i.e. blank starting rows where truncated surveys by site can be added
colnames(exifTruncated)<-colnames(exif) #Add column names to the new truncated dataset

for(i in 1:length(PlaceID)){ #for each camera site
  d<-subset(exif, TagPlaceID==PlaceID[i]) #select by that site
  for(j in 2:length(d$File)){ #for each successive image after the first image
    if(d$TimeDecimal[j]>=(d$TimeDecimal[1]+15/1440)){ #As data is in chronological order, and camera setups/takedowns occur during the day (so will not be affected by the return to TimeDecimal==0 at midnight), compare succesive image time stamps until one is greater than 15 minutes after the first image's timestamp. 
      newstart<-j #New starting image for the truncated dataset at that site is at jth row (jth image is the first one that is at least 15 minutes after the timestamp of the first image)
      break}} #once new starting image found, cease the search
  
  for(k in 1:length(d$File)){ #for each image before the final image
    if(identical(d$Date[length(d$File)],d$Date[length(d$File)-k])){ #if the date is the same for the images being compared
      if(d$TimeDecimal[length(d$File)-k]<=(d$TimeDecimal[length(d$File)]-15/1440)){   #compare successive images (from newest to oldest) to the final image, to determine the first image that was taken at least 15 minutes before the last image 
        newend<-length(d$File)-k #New final image for the truncated dataset at that site is at kth row (kth image is the first one that is at least 15 minutes before the timestamp of the last image)
        break}}else{ #once a new final image is found, cease the search
          newend<-length(d$File)#if on separate days, images will be greater than 15 minutes apart as camera set ups and take downs never occurred around midnight
        }
  }
  NewEnd<-d[newstart:newend,] #truncate the data set at that site by the new starting and final images found
  exifTruncated<-as.data.frame(rbind(exifTruncated, NewEnd)) #bind the truncated dataset at this site to ones at other sites
}
exifTruncated<-exifTruncated[-c(1,2),] #remove blank starting rows
rownames(exifTruncated)<-1:length(exifTruncated$File) #Rename the rows based on the number that are now present
View(exifTruncated) #check final result has formed correctly

#FILTER BY 'NIGHT HOURS' ONLY (6PM TO 8AM)
exifNIGHT<-subset(exifTruncated,TimeDecimal<=8/24| TimeDecimal>=18/24) #time expressed as a proportion of 24 hours

write.csv(exifNIGHT, "exifNIGHT.csv") #save new file

######## MEGADETECTOR PREPARATION ########
#FILTER BY 0.65 CONFIDENCE VALUE (See 'Confidence Value Threshold' R script)
#DO NOT FILTER BY 'PERSON' HERE AS NEED SITE WITH ZERO SIGHTINGS OF HUMANS TO BE RECORDED TOO
MD<-subset(MD, conf>=0.65)

#REMOVE SITES WITH POOR PLACEMENT, IMAGE CAPTURE AND INCONSISTENT TAGGING
MD<-subset(MD, TagPlaceID %in% SitesToKeep) 

#TRUNCATE SETUP/TAKE DOWN TIME FROM EXIFPRO DATA 
#I.E. REMOVE FIRST AND LAST 15 MINUTES FROM THE SURVEYS OF EACH SITE

#FOR ANNOTATIONS OF BELOW CODE, SEE EXIFPRO ABOVE.
PlaceIDMD<-sort(unique(MD$TagPlaceID))
MD$Date<-as.character(MD$Date)
MDTruncated<-as.data.frame(matrix(0, nrow=2, ncol=7)) 
colnames(MDTruncated)<-colnames(MD)

for(i in 1:length(PlaceIDMD)){
  d<-subset(MD, TagPlaceID==PlaceIDMD[i])
  for(j in 2:length(d$File)){
    if(d$TimeDecimal[j]>=(d$TimeDecimal[1]+15/1440)){
      newstart<-j 
      break}}
  
  for(k in 1:length(d$File)){
    if(identical(d$Date[length(d$File)],d$Date[length(d$File)-k])){
      if(d$TimeDecimal[length(d$File)-k]<=(d$TimeDecimal[length(d$File)]-15/1440)){ 
        newend<-length(d$File)-k 
        break}}else{
          newend<-length(d$File)
        }
  }
  New<-d[newstart:newend,] 
  MDTruncated<-as.data.frame(rbind(MDTruncated, New)) 
}
MDTruncated<-MDTruncated[-c(1,2),] 
rownames(MDTruncated)<-1:length(MDTruncated$File)
View(MDTruncated)

#DOUBLE CHECK SITES ARE THE SAME IN EXIF PRO AND MEGADETECTOR DATASETS
PlaceIDE<-sort(unique(exif$TagPlaceID)) #Vector of ExifPro camera sites
PlaceIDM<-sort(unique(MDTruncated$TagPlaceID)) #Vector of MegaDetector camera sites
PlaceIDE-PlaceIDM #Result is a vector of zeroes, so camera sites are the same in each dataset

#SAVE TIME DATA
write.csv(MDTruncated, "mdFULL.csv") #save MegaDetector file across the full 24 hours

#FILTER BY 'NIGHT HOURS' ONLY (6PM TO 8AM)
mdNIGHT<-subset(MDTruncated,TimeDecimal<=8/24 | TimeDecimal>=18/24) #time expressed as a proportion of 24 hours

write.csv(mdNIGHT, "mdNIGHT.csv") #save MegaDetector file across night hours

