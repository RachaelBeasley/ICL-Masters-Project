#DETERMINE THE PROPORTION OF NOCTURNAL MAMMALS SIGHTED
#ACROSS THE NIGHT (6PM TO 8AM) AND DAY (8AM TO 6PM)
#USES DR JEFFS WAAGE FULL DATASET FOR 27 CAMERA SITES

rm(list=ls())

#LOAD AND CHECK data
setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/Data") #set working directory
data<-read.csv("J.Waage Data-27 Sites.csv", header=TRUE) #Dr Jeff Waage's dataset with 24 hour tagging of all species
View(data)
str(data)

DAYNIGHT<-function(data){
  
  #CONVERT RELEVANT COLUMNS FROM FACTORS TO CHARACTERS
  data$TagSpecies1<-as.character(data$TagSpecies1)
  data$TagSpecies2<-as.character(data$TagSpecies2)
  data$TagSpecies3<-as.character(data$TagSpecies3)
  data$Date<-as.character(data$Date)
  
  #SPECIES TO EXAMINE
  Species<-c("Fox", "Hedgehog", "Human")
  
  #CREATE RESULTS MATRIX
  resultsDN<-matrix(0, ncol=3, nrow=3)#empty results matrix
  colnames(resultsDN)<-Species #row names for results matrix
  rownames(resultsDN)<-c("TotalImages", "DayProp", "NightProp") #column names for results matrix

 for(s in 1:3){ #for each species
   d<-subset(data, subset = TagSpecies1== Species[s] | TagSpecies2 == Species[s]  | TagSpecies3 == Species[s]) #filter data by tags for specified species
   
  #APPLY ONE MINUTE INTERVAL TO CORRECT FOR CAMERA DIFFERENCES
   if(length(d$File)>1){ #if there is at least two images to compare
     record<-1 #for the initial comparison, use the first image as the reference image
     
     repeat{ 
       i<-record[length(record)] #for the last image marked as being greater than one minute from the last reference image, mark it as the new reference image
       for(j in (i+1):length(d$File)){ #compare each successive image's timestamp to the reference image's timestamp
         if(d$Date[i]==d$Date[j]){ #if images were taken on the same date (images in chronological order to enable these comparisons)
           if(d$TimeDecimal[j]>=(d$TimeDecimal[i]+1/1440)){ #and if time taken was one minute or more after the reference image
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
     }else{ 
       if(length(d$File)==1){OM<-d}} #if only one image of that species, only this image is counted for species number
     
     OM<-d[record,] #new dataset contains only the images of specified species with at least one minute between them

   #TOTAL NUMBER OF IMAGES ADDED TO RESULTS MATRIX
   resultsDN[1,s]<-length(OM$File)
   
   #PROPROTION OF THESE IMAGES TAKEN DURING THE DAY (8AM TO 6PM)
   resultsDN[2,s]<-length(which(OM$TimeDecimal>8/24 & OM$TimeDecimal <18/24))/resultsDN[1,s]
   
   #PROPROTION OF THESE IMAGES TAKEN DURING THE NIGHT (6PM AND 8AM)
   resultsDN[3,s]<-1-resultsDN[2,s]
  setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/Miscellaneous") #set working directory to save results file
  write.csv(resultsDN, "Day and Night Images-By Species (Jeffs).csv") #save results matrix as a CSV files
 }
}

DAYNIGHT(data) #run function

