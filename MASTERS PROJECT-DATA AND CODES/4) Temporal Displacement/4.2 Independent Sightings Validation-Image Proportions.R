#TO JUSTIFY THE TIME INTERVAL FOR INDEPENDENT SIGHTINGS
#DETERMINE THE NUMBER OF SPECIES TAGS REMAINING RELATIVE TO TOTAL NUMBER OF IMAGES OF THAT SPECIES
#FOR EACH TIME INTERVAL USED AND FOR EACH SPECIES
#BASED ON METHODOLOGY FROM YASUDA (2004)

#AS PART OF YASUDA'S (2004) METHOD, VERY SMALL DIFFERENCES BETWEEN THE 
#PROPORTION OF IMAGES REMAINING ACROSS A RANGE OF TIME INTERVALS SHOW
#THAT ANY OF THE TIME INTERVALS, ACROSS THAT RANGE, CAN BE CONSIDERED AS 
#THE THRESHOLD FOR INDEPENDENT SIGHTINGS
#REFERENCE: Yasuda, M. (2004). Monitoring diversity and abundance of mammals with camera traps: A case study on Mount Tsukuba, central Japan. Mammal Study. 29(1): 37-46. doi: 10.3106/mammalstudy.29.37

rm(list=ls())

######## LOAD DATA ########
#LOAD EXIF PRO DATA (I.E. MAMMAL SIGHTINGS)
setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/Data") #set working directory to retrieve data
exifNIGHT<-read.csv("exifNIGHT.csv", header=TRUE) #load data with mammal tags
View(exifNIGHT)
str(exifNIGHT)

#LOAD MEGADETECTOR DATA (I.E. HUMAN SIGHTINGS)
mdFULL<-read.csv("mdFULL.csv" , header=TRUE) #load data with human tags
str(mdFULL)
View(mdFULL)

#DETERMINE PROPORTION OF TOTAL IMAGES REMAINING AT EACH TIME INTERVAL
ImageProp<-function(exifNIGHT, mdFULL){
  
  #CONVERT COLUMN VALUES FROM FACTORS TO CHARACTERS 
  exifNIGHT$TagSpecies1<-as.character(exifNIGHT$TagSpecies1)
  exifNIGHT$TagSpecies2<-as.character(exifNIGHT$TagSpecies2)
  exifNIGHT$TagSpecies3<-as.character(exifNIGHT$TagSpecies3)
  exifNIGHT$Date<-as.character(exifNIGHT$Date)
  mdFULL$category<-as.character(mdFULL$category)
  mdFULL$Date<-as.character(mdFULL$Date)
  
  #SPECIES TO EXAMINE, PUT INTO A VECTOR
  Species<-c("Fox", "Hedgehog", "Human")
  
  #CREATE RESULTS MATRIX
  RESULTS<-matrix(0, nrow=9, ncol=3) #empty results matrix
  rownames(RESULTS)<-c("Total Images", "Full","One Minute", "Two Minutes", "Five Minutes", "Ten Minutes", "Twenty Minutes", "Thirty Minutes", "Sixty Minutes") #name rows with the time intervals to be used
  colnames(RESULTS)<-c(Species) #name columns with species name
  
  #DETERMINE THE TOTAL NUMBER OF IMAGES TAKEN OF EACH SPECIES AND SAVE IN RESULTS VECTOR
  RESULTS[1,1]<-length(which(exifNIGHT$TagSpecies1 == "Fox"| exifNIGHT$TagSpecies2 == "Fox" | exifNIGHT$TagSpecies3 == "Fox"))
  RESULTS[1,2]<-length(which(exifNIGHT$TagSpecies1 == "Hedgehog"| exifNIGHT$TagSpecies2 == "Hedgehog"| exifNIGHT$TagSpecies3 == "Hedgehog"))
  RESULTS[1,3]<-length(which(mdFULL$category=="person"))
  
  RESULTS[2,]<-1 #As time intervals have not been applied yet, the proportion of images remaining is 100% (i.e. 1). Add this value to results matrix for plotting later
  
  #DETERMINE PROPORTION OF IMAGES REMAINING IN EACH TIME INTERVAL
  TimeInterval<-c(1,2,5,10,20,30,60) #vector of numerical time intervals to filter by
  
  for(k in 1:3){ #for each species
    if(k<3){d<-subset(exifNIGHT, subset = TagSpecies1== Species[k] | TagSpecies2 == Species[k]  | TagSpecies3 == Species[k])} #filter by the specified mammal
    if(k==3){d<-subset(mdFULL, category=="person")} #filter by human tags once mammal tags have been processed
    
    #APPLY TIME INTERVAL
    for(t in 1:7){ #for each time interval  
      
      if(length(d$File)>1){ #if there are at least two images to compare
        record<-1 #for the initial comparison, use the first image as the reference image
        repeat{
          i<-record[length(record)] #for the last image marked as being greater than 'TimeInterval[t]' minutes from the last reference image, mark it as the new reference image
          for(j in (i+1):length(d$File)){ #compare each successive image's timestamp to the reference image's timestamp
            if(d$Date[i]==d$Date[j]){ #if images were taken on the same date (images in chronological order to enable these comparisons)
              if(d$TimeDecimal[j]>=(d$TimeDecimal[i]+TimeInterval[t]/1440)){ #and if time taken was 'TimeInterval[t]' minutes or more after the reference image
                record<-c(record,j) #record image as the new reference image and stop comparisons with the current reference image
                break}}else{ #else if images were taken on different days
                  if(d$TimeDecimal[i]>=1-TimeInterval[t]/1440 & d$TimeDecimal[j]<=TimeInterval[t]/1440){ #if reference image and comparison image were taken within 'TimeInterval[t]' minutes of midnight (assuming no day without at least one image)
                    if(((1-d$TimeDecimal[i])+d$TimeDecimal[j])>=TimeInterval[t]/1440){  #and if the time difference between the images being taken is greater than 'TimeInterval[t]' minutes 
                      record<-c(record,j) #make the next day image the new reference image, and stop comparisons with the current reference image
                      break}}else{
                        record<-c(record,j) #else time difference is larger than 'TimeInterval[t]' minutes between images taken on different days, so next day image becomes the new reference image
                        break}
                }
          }
          if(j==length(d$File)){break} #once comparisons have been made until the very last image, stop the entire procedure
          }
        p<-d[record,] #new dataset for all camera sites contains only the images with at least TimeInterval[t] minute's between them
        RESULTS[t+2,k]<-length(p$File)/RESULTS[1,k] #calculate the number of images remaining after time interval is applied/total images. Save this in the results matrix
      }else{if(length(d$File)==1){RESULTS[t+2,k]<-1/RESULTS[1,k]}} #if only one image of that species at a specific time interval, divide 1 by the total number of images
    }
  }
  RESULTS<-as.data.frame(RESULTS) #convert results matrix to a dataframe to allow combination of letters and numbers
  RESULTS$TimeInterval.no<-c("NA",0,1,2,5,10,20,30,60) #Add Time Interval in Numbers for plotting purposes
  RESULTS$MeanImages<-apply(RESULTS[,1:3],1, mean)#determine mean proportion of images across all species remaining at each time interval
  RESULTS$MinImages<-apply(RESULTS[,1:3],1, min) #determine minimum proportion of images across all species remaining at each time interval
  RESULTS$MaxImages<-apply(RESULTS[,1:3],1, max) #determine maximum proportion of images across all species remaining at each time interval
  View(RESULTS) #view final results dataframe
  
  setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/4) Temporal Displacement")
  write.csv(RESULTS, "Independent Sightings Validation.csv") #save results matrix as a CSV file
}

ImageProp(exifNIGHT, mdFULL) #run function

##### PLOTTING DATA AND DETERMINING GRADIENTS ####
#Not presented in write up, but useful for visualisation and further confirmation of the time 
#interval used to assume independent sightings

#LOAD DATA
#N.B 'Time Interval' added as column label for the first column in Excel
setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/4) Temporal Displacement") #set working directory to retrieve results file
RESULTS<-read.csv("Independent Sightings Validation.csv", header=TRUE) #load results file
View(RESULTS)

#PLOT DATA: PROPORTION OF IMAGES OF EACH SPECIES REMAINING AT EACH TIME INTERVAL
jpeg("Image Proportions at different time intervals.jpeg") #prime R to save a JPEG image
par(mar=c(4,4,1,1)) #change margins of graphics window to fit the figure

#PLOT IMAGE PROPORTIONS OF THE TOTAL NUMBER OF IMAGES FOR EACH SPECIES AT EACH TIME INTERVAL ('plot')
#ADD LINES BETWEEN THE DATA POINTS IN A DIFFERENT COLOUR FOR EACH SPECIES
#FOR FOXES
plot(RESULTS$Fox[2:9]~as.numeric(RESULTS$TimeInterval.no[2:9]), col="orchid3",pch=19,cex=0.5, main="Image Proportions vs Time Intervals", xlab="Time Interval (minutes)", ylab="Proportion of Images Remaining")
lines(RESULTS$Fox[2:9]~as.numeric(RESULTS$TimeInterval.no[2:9]), col="orchid3")
#FOR HEDGEHOGS
points(RESULTS$Hedgehog[2:9]~as.numeric(RESULTS$TimeInterval.no[2:9]), col="peru",pch=19,cex=0.5, main="Image Proportions vs Time Intervals", xlab="Time Interval (minutes)", ylab="Proportion of Images Remaining")
lines(RESULTS$Hedgehog[2:9]~as.numeric(RESULTS$TimeInterval.no[2:9]), col="peru")
#FOR HUMANS
points(RESULTS$Human[2:9]~as.numeric(RESULTS$TimeInterval.no[2:9]), col="black",pch=19,cex=0.5, main="Image Proportions vs Time Intervals", xlab="Time Interval (minutes)", ylab="Proportion of Images Remaining")
lines(RESULTS$Human[2:9]~as.numeric(RESULTS$TimeInterval.no[2:9]), col="black")
#ADD LEGEND TO SHOW WHICH COLOUR REPRESENTS WHICH SPECIES
legend("topright", legend=c("Fox", "Hedgehog", "Human"), col=c("orchid3", "peru", "black"), lty=1, lwd=2, cex=1)
dev.off() #FINISH SAVING THE JPEG

#CALCULATE GRADIENTS BETWEEN EACH CONSECUTIVE DATA POINT 
#AS PART OF YASUDA'S (2004) METHOD, VERY SMALL GRADIENTS BETWEEN DATA POINTS SHOW THERE IS 
#LITTLE DIFFERENCE BETWEEN THE PROPORTION OF IMAGES REMAINING BETWEEN THOSE TIME INTERVALS.
#IF GRADIENT REMAINS CONSISTENTLY SMALL, ANY OF THE TIME INTERVALS ACROSS THE GRADIENT 
#CAN BE CONSIDERED AS THE THRESHOLD FOR INDEPENDENT SIGHTINGS
gradients<-matrix(0, nrow=7, ncol=3) #empty results matrix for gradient values
Species<-c( "Fox", "Hedgehog", "Human") #vector of species examined
colnames(gradients)<-Species #column names of results matrix put as species names
rownames(gradients)<-c("0 to 1", "1 to 2", "2 to 5", "5 to 10", "10 to 20", "20 to 30", "30 to 60") #add row names to results matrix showing the time intervals the gradient has been calculated between
View(gradients) #view empty results matrix

for(k in 1:7){ #for each time interval
  for(i in 1:3){ #for each species
    gradients[k,i]<-(RESULTS[k+2, i+1]-RESULTS[k+1,i+1])/(as.numeric(RESULTS[k+2,4])-as.numeric(RESULTS[k+1,4])) #calculate the gradient as (change in proportion of images remaining)/(change in time interval)
  }
  }

gradients<-as.data.frame(gradients) #convert gradients results matrix to a data frame to allow apply function to be used
gradients$MeanGdt<-apply(gradients[,1:3],1, mean)#Mean gradient across each time interval
gradients$MinGdt<-apply(gradients,1, min) #minimum gradient across each time interval
gradients$MaxGdt<-apply(gradients,1, max) #maximum gradient across each time interval
View(gradients) #view completed gradients results matrix 
write.csv(gradients, "Image Proportions-Gradients.csv") #save completed results matrix as CSV file

#PLOT DATA WITH ONE MINUTE REMOVED
#THIS IS WHERE THE MOST DRAMATIC DROP OCCURS, SO REMOVING IT WILL ALLOW THE SMALLER GRADIENTS
#BETWEEN DATA POINTS TO BE SEEN
jpeg("Image Proportions with 1 Minute Removed.jpeg") #prime R to save a JPEG
par(mar=c(4,4,1,1)) #change margins of graphics window to fit the figure
#PLOT IMAGE PROPORTIONS OF THE TOTAL NUMBER OF IMAGES FOR EACH SPECIES AT EACH TIME INTERVAL ('plot')
#ADD LINES BETWEEN THE DATA POINTS IN A DIFFERENT COLOUR FOR EACH SPECIES
#FOR FOXES
plot(RESULTS$Fox[3:9]~as.numeric(RESULTS$TimeInterval.no[3:9]), col="orchid3",pch=19,cex=0.5, main="Image Proportions vs Time Intervals", xlab="Time Interval (minutes)", ylab="Proportion of Images Remaining")
lines(RESULTS$Fox[3:9]~as.numeric(RESULTS$TimeInterval.no[3:9]), col="orchid3")
#FOR HEDGEHOGS
points(RESULTS$Hedgehog[3:9]~as.numeric(RESULTS$TimeInterval.no[3:9]), col="peru",pch=19,cex=0.5, main="Image Proportions vs Time Intervals", xlab="Time Interval (minutes)", ylab="Proportion of Images Remaining")
lines(RESULTS$Hedgehog[3:9]~as.numeric(RESULTS$TimeInterval.no[3:9]), col="peru")
#FOR HUMANS
points(RESULTS$Human[3:9]~as.numeric(RESULTS$TimeInterval.no[3:9]), col="black",pch=19,cex=0.5, main="Image Proportions vs Time Intervals", xlab="Time Interval (minutes)", ylab="Proportion of Images Remaining")
lines(RESULTS$Human[3:9]~as.numeric(RESULTS$TimeInterval.no[3:9]), col="black")
#ADD LEGEND TO SHOW WHICH COLOUR REPRESENTS WHICH SPECIES
legend("topright", legend=c("Fox", "Hedgehog","Human"), col=c("orchid3","peru", "black"), lty=1, lwd=2, cex=1)
dev.off() #FINISH SAVING THE JPEG
