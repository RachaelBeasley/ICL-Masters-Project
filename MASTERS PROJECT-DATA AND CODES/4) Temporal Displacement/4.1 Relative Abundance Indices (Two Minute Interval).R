#NB THE ORIGINAL DATASETS DO NOT HAVE PERMISSION TO BE SHARED
#THESE ORIGINAL DATASETS ARE NEEDED FOR THIS CODE TO RUN

#RELATIVE ABUNDANCE INDEX
#=NUMBER OF SPECIES TAGS/NUMBER OF DAYS RECORDED
#TWO MINUTE TIME INTERVAL APPLIED TO DELIMIT INDEPENDENT SIGHTINGS
#MAMMAL TAGS BETWEEN 6PM TO 8AM
#HUMAN TAGS FOR FULL 24 HOURS

rm(list=ls())

######## LOAD AND PREPARATION OF DATA ########
#LOAD EXIF PRO DATA (I.E. MAMMAL SIGHTINGS)
setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/Data")
exifNIGHT<-read.csv("exifNIGHT.csv", header=TRUE) 
str(exifNIGHT) #examine data structure
View(exifNIGHT) #visualise examine data

#LOAD MEGADETECTOR DATA (I.E. HUMAN SIGHTINGS)
mdFULL<-read.csv("mdFULL.csv" , header=TRUE) 
str(mdFULL)
View(mdFULL)

#DETERMINE RELATIVE ABUNDANCE INDEX FOR EACH SPECIES AT EACH CAMERA SITE. 
#TWO MINUTE TIME INTERVAL TO DELIMIT INDEPENDENT SIGHTINGS
findRAITM<-function(exifNIGHT, mdFULL){
  
  #DETERMINE NUMBER OF DAYS RECORDED AT EACH CAMERA SITES
  PlaceID<-sort(unique(mdFULL$TagPlaceID)) #Get vector of camera sites in the datasets
  Days<-rep(0, length(PlaceID)) #Empty results vector to record number of days recorded at each camera sites
  for(i in 1:length(PlaceID)){ #For each camera site
    d<-subset(exifNIGHT, TagPlaceID==PlaceID[i]) #Filter by that camera site
    Days[i]<-length(unique(d$Date))} #Record the number of days recorded at that site
  
  #CONVERT RELEVANT COLUMN VALUES TO CHARACTERS FROM FACTORS
  exifNIGHT$TagSpecies1<-as.character(exifNIGHT$TagSpecies1)
  exifNIGHT$TagSpecies2<-as.character(exifNIGHT$TagSpecies2)
  exifNIGHT$TagSpecies3<-as.character(exifNIGHT$TagSpecies3)
  exifNIGHT$Date<-as.character(exifNIGHT$Date)
  mdFULL$category<-as.character(mdFULL$category)
  mdFULL$Date<-as.character(mdFULL$Date)
  
  #VECTOR OF MAMMAL SPECIES TO SUBSET BY
  #N.B.DOGS ARE ALSO INCLUDED TO DETERMINE CORRELATIONS WITH HUMANS IN LATER ANALYSES
  Species<-c("Dog","Fox", "Hedgehog")
  
  #CREATE RESULTS MATRIX FOR RELATIVE ABUNDANCE INDICES 
  RAIResultsTM<-matrix(0, nrow=length(PlaceID), ncol=length(Species)+1) #empty results matrix
  colnames(RAIResultsTM)<-c("Human", Species) #name columns with species names
  rownames(RAIResultsTM)<-PlaceID #name rows with camera site numbers
  
  #FIRST, DETERMINE RELATIVE ABUNDANCE INDICES FOR HUMANS
  d<-subset(mdFULL, category=="person") #filter MegaDetector dataset by human tags
  
  #APPLY TWO MINUTE TIME INTERVAL
  for(k in 1:length(PlaceID)){ #for each camera site
    p<-subset(d, TagPlaceID==PlaceID[k]) #filter data by that camera site
    
    if(length(p$File)>1){ #if there are at least two images to compare
      
      record<-1 #for the initial comparison, use the first image as the reference image
      repeat{
        i<-record[length(record)] #for the last image marked as being greater than two minutes from the last reference image, mark it as the new reference image
        for(j in (i+1):length(p$File)){ #compare each successive image's timestamp to the reference image's timestamp
          if(p$Date[i]==p$Date[j]){ #if images were taken on the same date (images in chronological order to enable these comparisons)
            if(p$TimeDecimal[j]>=(p$TimeDecimal[i]+2/1440)){ #and if time taken was two minutes or more after the reference image
              record<-c(record,j) #record image as the new reference image and stop comparisons with the current reference image
              break}}else{  #else if images were taken on different days
                
                if(p$TimeDecimal[i]>=1-2/1440 & p$TimeDecimal[j]<=2/1440){ #if image on the first day was taken later than 11:58pm and if the image taken on the second day was taken before 0:02am (assuming no day without at least one image)
                  if(((1-p$TimeDecimal[i])+p$TimeDecimal[j])>=2/1440){ #and if the time difference between the images being taken is greater than two minutes 
                    record<-c(record,j) #make the next day image the new reference image, and stop comparisons with the current reference image
                    break}}else{
                      record<-c(record,j) #else time difference is larger than two minutes between images taken on different days, so next day image becomes the new reference image
                      break} 
              }
        }
        if(j==length(p$File)){break} #once comparisons have been made until the very last image, stop the entire procedure
      }
      RAIResultsTM[k,1]<-length(record)/Days[k] #Calculate Relative Abundance Index of humans at the camera site, and save in results matrix (number of human tags/number of days recorded at each site).
    }else{if(length(p$File)==1){RAIResultsTM[k,1]<-1/Days[k]}} #else if only one image of a human at that camera site, calculate Relative Abundance Index by dividing 1 by number of days recorded at that site
  }
  
  
  #SECOND, DETERMINE RELATIVE ABUNDANCE INDICES FOR NON-HUMAN SPECIES
  for(s in 1:length(Species)){ #for each non-human species
    d<-subset(exifNIGHT, subset = TagSpecies1== Species[s] | TagSpecies2 == Species[s]  | TagSpecies3 == Species[s]) #filter ExifPro data by specified species
    
    #APPLY TWO MINUTE TIME INTERVAL
    #FOR ANNOTATIONS, SEE HUMANS ABOVE AS PROCEDURE REPEATED
    for(k in 1:length(PlaceID)){ 
      p<-subset(d, TagPlaceID==PlaceID[k])
      
      if(length(p$File)>1){
        
        record<-1
        repeat{
          i<-record[length(record)]
          for(j in (i+1):length(p$File)){
            if(p$Date[i]==p$Date[j]){ 
              if(p$TimeDecimal[j]>=(p$TimeDecimal[i]+2/1440)){ 
                record<-c(record,j)
                break}}else{ 
                  if(p$TimeDecimal[i]>=1-2/1440 & d$TimeDecimal[j]<=2/1440){ 
                    if(((1-p$TimeDecimal[i])+p$TimeDecimal[j])>=2/1440){  
                      record<-c(record,j)
                      break}}else{
                        record<-c(record,j)
                        break} 
                }
          }
          if(j==length(p$File)){break}
        }
        
        RAIResultsTM[k,s+1]<-length(record)/Days[k] #Calculate Relative Abundance Index of specified mammal at the camera site, and save in results matrix (number of mammal tags/number of days recorded at each site).
      }else{if(length(p$File)==1){RAIResultsTM[k,s+1]<-1/Days[k]}} #else if only one image of the specified mammal at that camera site, calculate Relative Abundance Index by dividing 1 by number of days recorded at that site
    }
  }
  
  setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/4) Temporal Displacement")
  write.csv(RAIResultsTM, "RAI Results MD EXIF-Two Minutes.csv") #Save completed results matrix as a CSV file
}

findRAITM(exifNIGHT, mdFULL) #run function

#N.B. RESULTS FILE NEEDS 'PlaceID' ADDING THE NAME OF THE FIRST COLUMN