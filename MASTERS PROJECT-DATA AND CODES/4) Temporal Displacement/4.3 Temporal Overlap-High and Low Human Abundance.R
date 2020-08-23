#TEMPORAL OVERLAP CALCULATED BETWEEN HUMAN AND EACH MAMMAL IN SITES OF LOW AND HIGH HUMAN ABUNDANCE
#USING 24 HOUR HUMANS AND NOCTURNAL SPECIES (FOXES AND HEDGEHOGS)
#CORRECTED FOR SUNRISE AND SUNSET TIMES
#With bootstrapping to get confidence intervals, code takes approximately 10 minutes to run

rm(list=ls())

######## LOAD AND PREPARATION OF DATA ########
#LOAD EXIF PRO DATA (I.E. MAMMAL SIGHTINGS)
setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/Data") #set working directory
exifNIGHT<-read.csv("exifNIGHT.csv", header=TRUE) #load mammal tag data
str(exifNIGHT) #examine data structure
View(exifNIGHT) #visualise examine data

#LOAD MEGADETECTOR DATA (I.E. HUMAN SIGHTINGS)
mdFULL<-read.csv("mdFULL.csv" , header=TRUE) 
str(mdFULL)
View(mdFULL)

#LOAD RAI VALUES TO SEPERATE OUT HIGH AND LOW HUMAN ABUNDANCE
setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/4) Temporal Displacement")
RAI<-read.csv("RAI Results MD EXIF-Two Minutes.csv", header=TRUE)
str(RAI)
View(RAI)

#DETERMINE THRESHOLD RELATIVE ABUNDANCE VALUE TO DISTINGUISH SITES OF LOW AND HIGH HUMAN ABUNDANCE
#CHECK FOR NATURAL BREAKPOINTS
range(RAI$Human) #Range of Human RAI: 0 to 176.8. Used to determine number of bins
#PLOT HISTOGRAM WITH BIN INTERVALS OF ONE TO SEE DISTRIBUTION OF RAIs
hist(RAI$Human, breaks=177, col="red", xlab="Human RAI", main="Distribution of Human RAI") 

#MAJORITY OF RAIs ARE BELOW 25
#EXAMINE RAIs ON A FINER SCALE BY SEEING THOSE ONLY BELOW 25
RAI25<-subset(RAI, Human<=25) #Subset RAIs by those less than 25
range(RAI25$Human) #New range of Human RAIs: 0 to 24.35714
#PLOT HISTOGRAM WITH BIN INTERVALS OF ONE TO SEE DISTRIBUTION OF RAIs
hist(RAI25$Human, breaks=25,col="red", xlab="Human RAI", main="Distribution of Human RAI") 

#MAJORITY OF RAIS ARE BELOW 5
#EXAMINE ON A FINER SCALE BY SEEING THOSE ONLY BELOW 5
RAI5<-subset(RAI, Human<=5) #Subset RAIs by those less than 5
range(RAI5$Human) #New range of Human RAIs: 0 to 5
#PLOT HISTOGRAM WITH BIN INTERVALS OF ONE TO SEE DISTRIBUTION OF RAIs
hist(RAI5$Human, breaks=5,col="red", xlab="Human RAI", main="Distribution of Human RAI") 

#OVERALL, THERE IS A CLEAR PEAK WITH RAIs LESS THAN 1. THIS PEAK IS FOLLOWED BY A SHARP DROP IN FREQUENCIES.
#THIS SUGGESTS THAT AN RAI OF ONE IS AN APPROPRIATE THRESHOLD TO DISTINGUISH BETWEEN SITES 
#OF LOW AND HIGH HUMAN ABUNDANCE


##### TEMPORAL OVERLAP FUNCTION ########
#POOLS HUMAN AND SPECIES SIGHTINGS BY EACH SITE OF EITHER HIGH OR LOW HUMAN ABUNDANC
#COUNTS TOTAL NUMBER OF SIGHTINGS OF EACH SPECIES IN SITES OF EITHER HIGH OR LOW HUMAN ABUNDANCE
#CALCULATES TEMPORAL OVERLAP BETWEEN HUMANS AND EACH SPECIES IN SITES OF EITHER HIGH OR LOW HUMAN ABUNDANCE, WITH LOWER AND UPPER CIs
#PLOTS OVERLAP GRAPHS FOR EACH SPECIES IN IN SITES OF EITHER HIGH OR LOW HUMAN ABUNDANCE

#INSTALL AND LOAD "OVERLAP" PACKAGE
install.packages("overlap") #for estimating temporal overlap
library("overlap")

install.packages("maptools") #for conversion of time-in-radians to SunTime (i.e. times are standardised to account for the changing sunrise and sunset time across 2.5 months)
library("maptools")
install.packages("sp")
library("sp")

HLtemporalOVERLAP<-function(exifNIGHT, mdFULL, RAI){
  
  #CONVERT COLUMNS TO FORMATS NEEDED LATER IN THE CODE
  exifNIGHT$Radians<-exifNIGHT$TimeDecimal*2*pi #Converts time-in-decimals to radians. Necessary for OVERLAP package
  exifNIGHT$TagSpecies1<-as.character(exifNIGHT$TagSpecies1) #Convert species tags to characters from factors
  exifNIGHT$TagSpecies2<-as.character(exifNIGHT$TagSpecies2)
  exifNIGHT$TagSpecies3<-as.character(exifNIGHT$TagSpecies3)
  mdFULL$Radians<-mdFULL$TimeDecimal*2*pi #Converts time-in-decimals to radians. Necessary for OVERLAP package
  mdFULL$category<-as.character(mdFULL$category) #Convert human tags to characters from factors
  
  #CONVERT TIMES TO SUNTIME (ACCOUNTING FOR CHANGING PHOTOPERIOD ACROSS THE TOTAL SURVEY PERIOD)
  exifNIGHT$DatesP <- as.POSIXct(exifNIGHT$Date, tz="GMT") ##Convert dates to a POSIXct object
  mdFULL$DatesP <- as.POSIXct(mdFULL$Date, tz="GMT") ##Convert dates to a POSIXct object
  
  #DEFINE COORDINATES OF HAMPSTEAD HEATH AS A SPATIAL POINT, SO PACKAGE KNOWS WHICH SUNRISE AND SUNSET TIMES TO CORRECT FOR
  coords<-matrix(c(51.5608,-0.1629), nrow=1) #Coordinates of the centre of Hampstead Heath (Site is only 275 ha, so sunrise and sunset will vary negligibly across sites)
  Coords <- SpatialPoints(coords, proj4string=CRS("+proj=longlat +datum=WGS84")) #convert coordinates to a spatial point object
  
  #FINAL CONVERSION TO SUNTIME (SUNRISE MAPPED TO PI/2 AND SUNSET TO 3*PI/2)
  exifNIGHT$SunTime<-sunTime(exifNIGHT$Radians, exifNIGHT$DatesP, Coords) #Convert time-in-radians to sunTime for ExifPro data
  mdFULL$SunTime<-sunTime(mdFULL$Radians, mdFULL$DatesP, Coords) #Convert time-in-radians to sunTime for MegaDetector data
  
  
  #CONSTRUCTION OF RESULTS MATRIX
  RESULTS<-as.data.frame(matrix(0, ncol=6, nrow=6)) #Create empty results matrix
  colnames(RESULTS)<-c("Human Abundance", "Species", "Total Sightings", "Overlap", "LowerCI", "UpperCI") #Name columns of new results matrix
  HIGHLOW<-c("High", "Low") #Create reference vector defining Human Abundances to group by
  SpeciesNames<-c("Human", "Fox", "Hedgehog") #Vector of species for row names and later filtering
  RESULTS[,1]<-rep(HIGHLOW, 3) #Column 1 ("Human Abundance") filled systematically with Human Abundance names
  RESULTS[,2]<-rep(SpeciesNames, each=2) #Column 2 ("Species") filled with systematically with Species names
  
  #FILTER OUT HIGH AND LOW ABUNDANCE SITES
  HighRAI<-subset(RAI, Human>=1) #An RAI of at least one human a day denotes a site as having High Human Abundance (see above before temporal overlap function)
  HighSites<-sort(unique(HighRAI$PlaceID))#Vector of the 71 sites of High Human Abundance
  
  LowRAI<-subset(RAI, Human<1) #An RAI of less than one human a day denotes a site as having Low Human Abundance (see above before temporal overlap function)
  LowSites<-sort(unique(LowRAI$PlaceID))##Vector of the 43 sites of Low Human Abundance 
  
  #HUMAN-APPLY ONE MINUTE TIME INTERVAL TO EACH SITE AND RECOMBINE TO FORM NEW DATASET 
  #ONE MINUTE TIME INTERVAL IS APPLIED TO CORRECT FOR CAMERA MODEL DIFFERENCES IN PHOTOGRAPHIC RATE
  # I.E. ONLY RECOGNISED TAGS THAT ARE AT LEAST ONE MINUTE BETWEEN EACHOTHER
  #TO CORRECT FOR CAMERA DIFFERENCES
  Human<-subset(mdFULL, category=="person") #Filter MegaDetector data for human tags
  
  PlaceID<-sort(unique(mdFULL$TagPlaceID)) #Vector of camera sites to examine
  HumanOM<-as.data.frame(matrix(0, nrow=2, ncol=11)) #starting point for resulting dataset with one minute time interval applied
  colnames(HumanOM)<-colnames(Human) #Name columns of new dataset with column names from the full dataset
  
  for(k in 1:length(PlaceID)){ #for each camera site
    p<-subset(Human, TagPlaceID==PlaceID[k]) #select each camera site
    
    if(length(p$File)>1){ #if there is at least two images to compare
      
      record<-1 #for the initial comparison, use the first image as the reference image
      repeat{
        i<-record[length(record)] #for the last image marked as being greater than one minute from the last reference image, mark it as the new reference image
        for(j in (i+1):length(p$File)){ #compare each successive image's timestamp to the reference image's timestamp
          if(p$Date[i]==p$Date[j]){ #if images were taken on the same date (images in chronological order to enable these comparisons)
            if(p$TimeDecimal[j]>=(p$TimeDecimal[i]+1/1440)){ #and if time taken was one minute or more after the reference image
              record<-c(record,j) #record image as the new reference image and stop comparisons with the current reference image
              break}}else{ #if images were taken on different days
                
                if(p$TimeDecimal[i]>=1-1/1440 & p$TimeDecimal[j]<=1/1440){ #if image on the first day was taken later than 11:59pm and if the image taken on the second day was taken before 0:01am (assuming no day without at least one image)
                  if(((1-p$TimeDecimal[i])+p$TimeDecimal[j])>=1/1440){ #if the time difference between the images being taken is greater than one minute 
                    record<-c(record,j) #make the next day image the new reference image, and stop comparisons with the current reference image
                    break}}else{
                      record<-c(record,j) #else time difference is larger than one minute between images taken on different days, so next day image becomes the new reference image
                      break} 
              }
        }
        if(j==length(p$File)){break} #once comparisons have made until the very last image, stop the entire procedure
      }
      new<-p[record,] #new dataset for that camera site contains only the images with at least one minute between them
      HumanOM<-as.data.frame(rbind(HumanOM, new))}else{ #bind filtered dataset for that camera site to those from other camera sites
        if(length(p$File)==1){HumanOM<-as.data.frame(rbind(HumanOM, p))}} #if only one image at that particular camera site, add that image to the new combined dataset
  }
  HumanOM<-HumanOM[-c(1,2),] #remove blank starting rows
  rownames(HumanOM)<-1:length(HumanOM$File) #Rename the rows based on the number that are now present
  
  RESULTS[1:2,4:6]<-1 #fill in results with human vs human temporal overlap values i.e. 1 (meaning complete temporal overlap)
  
  #MAMMALS-APPLY ONE MINUTE TIME INTERVAL TO EACH SITE AND RECOMBINE TO FORM NEW DATASET 
  #FOR ANNOTATIONS SEE ABOVE FOR WHEN HUMAN TAGS WERE FILTERED BY ONE MINUTE INTERVALS
  for(s in 2:3){
    animal<-subset(exifNIGHT, TagSpecies1 == SpeciesNames[s] | TagSpecies2 == SpeciesNames[s] | TagSpecies3 == SpeciesNames[s]) #filter ExifPro dataset for the tags of either foxes or hedgehogs
    
    #APPLY ONE MINUTE INTERVAL TO TAGS OF SPECIFIED MAMMAL
    PlaceID<-sort(unique(exifNIGHT$TagPlaceID))
    OM<-as.data.frame(matrix(0, nrow=2, ncol=18)) 
    colnames(OM)<-colnames(exifNIGHT)
    
    for(k in 1:length(PlaceID)){ 
      p<-subset(animal, TagPlaceID==PlaceID[k])
      
      if(length(p$File)>1){
        
        record<-1
        repeat{
          i<-record[length(record)]
          for(j in (i+1):length(p$File)){
            if(p$Date[i]==p$Date[j]){ 
              if(p$TimeDecimal[j]>=(p$TimeDecimal[i]+1/1440)){ 
                record<-c(record,j)
                break}}else{ 
                  if(p$TimeDecimal[i]>=1-1/1440 & p$TimeDecimal[j]<=1/1440){ 
                    if(((1-p$TimeDecimal[i])+p$TimeDecimal[j])>=1/1440){ 
                      record<-c(record,j)
                      break}}else{
                        record<-c(record,j)
                        break} 
                }
          }
          if(j==length(p$File)){break}
        }
        new<-p[record,]
        OM<-as.data.frame(rbind(OM, new))
      }else{if(length(p$File)==1){OM<-as.data.frame(rbind(OM, p))}}
    }
    
    OM<-OM[-c(1,2),] 
    rownames(OM)<-1:length(OM$File)
   
    #WITH ONE MINUTE TIME INTERVAL APPLIED, SUBSET BY HIGH OR LOW ABUNDANCE SITES FOR HUMAN AND EACH MAMMAL
    #DETERMINE NUMBER OF TAGS, AND OVERLAP INDEX WITH CONFIDENCE INTERVALS AND PLOT GRAPHS
    for(j in 1:2){ #for either High or Low Human Abundance 
      if(j==1){h<-subset(HumanOM, TagPlaceID %in% HighSites)} #if j==1, filter human tags by those from High Human Abundance Sites
      if(j==1){a<-subset(OM, TagPlaceID %in% HighSites)} #if j==1, filter specified mammal tags by those from High Human Abundance Sites
      if(j==2){h<-subset(HumanOM, TagPlaceID %in% LowSites)} #if j==2, filter human tags by those from Low Human Abundance Sites
      if(j==2){a<-subset(OM, TagPlaceID %in% LowSites)} #if j==2, filter specified mammal tags by those from Low Human Abundance Sites
      RESULTS[j,3]<-length(h$File) #record number of human tags in either High or Low Human Abundance Sites
      RESULTS[j+(s-1)*2,3]<-length(a$File) #record number of specified mammal tags in either High or Low Human Abundance Sites
      
      if(RESULTS[j,3]>=10 & RESULTS[j+(s-1)*2,3]>=10){ #only generate temporal overlap estimates if at least ten tags of human and foxes/hedgehogs
        
        if(RESULTS[j,3]>=75 & RESULTS[j+(s-1)*2,3]>=75){  #if there are at least 75 tags, for both human and mammal tags, 
          t="Dhat4"; ad=1}else{ #use temporal overlap estimator Dhat4
            t="Dhat1"; ad=0.8} #if either species has less than 75 tags, use temporal overlap estimator Dhat1
        
        overlapINDEX<-overlapEst(h$SunTime, a$SunTime, type=t) #calculate temporal overlap estimate using sunTimes
        RESULTS[j+(s-1)*2,4]<-overlapINDEX #add temporal overlap estimate to results matrix
        
        Bootstrap<-bootstrap(h$SunTime, a$SunTime, nb=1000, smooth=TRUE, adjust=ad, type=t) #calculate bootstrap values for temporal overlap estimates
        CI<-bootCI(overlapINDEX, Bootstrap, conf=0.95) #calculate confidence intervals using the bootstrap values
        RESULTS[j+(s-1)*2,5]<-CI[4,1] #add lower confidence interval to results matrix
        RESULTS[j+(s-1)*2,6]<-CI[4,2] #add upper confidence interval to results matrix
     
        #PLOTTING OVERLAP GRAPHs
        #PLOT GRAPHS ACROSS THE FULL 24 HOURS
        setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/4) Temporal Displacement") #change working directory to Results folder
        if(j==1){t<-"High Human Abundance"}; if(j==2){t<-"Low Human Abundance"} #select the Human Abundance presently being processed, in order to add the correct titles 
        jpeg(paste(SpeciesNames[s], "vs Human Temporal Overlap-", t, ".jpeg", sep=" "), width = 175, height = 150, units='mm', res = 300) #prime R to save temporal graph as a JPEG
        title= paste(SpeciesNames[s], "vs Human Temporal Overlap", sep= " ") #generate title for the temporal overlap plot
        #Graphs not changed to suntime as time-in-radians is easier to view
        overlapPlot(h$Radians, a$Radians, xcenter = c("noon"), linewidth=c(2,2), linecol=c("darkblue", "red"), olapcol="lightgray", main=title, ylab="Proportion of Activity", xlab="Time of Day", ylim=c(0,0.165)) #plot temporal overlap graph, with graph centred on noon, different coloured and textured line by species, and appropriate labels
        grid(ny=NULL,nx=NA, lty="dotted",col="lightgray") #add gridlines to plot
        legend('top', c("Human", SpeciesNames[s]), lty=c(1,2), lwd=c(2,2), col=c(4,2), horiz=TRUE, box.col="black", bg="white")  #add legend at the top of the plot
        mysubtitle= t; mtext(side=3,line=0.3, mysubtitle) #add subtitle defining Human Abundance used for the graph
        dev.off()
        
        #FOR BETTER VISUALS, ALSO PLOT GRAPHS FOCUSED BETWEEN 6PM AND 8AM WHEN NOCTURNAL MAMMALS ARE ACTIVE
        jpeg(paste(SpeciesNames[s], "vs Human Temporal Overlap-", t, "night.jpeg", sep=" "), width = 175, height = 150, units='mm', res = 300) #prime R to save temporal graph as a JPEG
        axisvalues<-c("17:00", "18:00","19:00", "20:00", "21:00", "22:00", "23:00", "00:00", "01:00", "02:00", "03:00", "04:00", "05:00", "06:00", "07:00", "08:00", "09:00") #manually define values for x axis, as overlap defaults plotting the full 24 hours and the 14 hours of interest
        tickpos<-c(-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9) #define the tick positions of the new x axis values
        title= paste(SpeciesNames[s], "vs Human Temporal Overlap", sep= " ") #generate title for the temporal overlap plot
        #Graphs not changed to suntime as time-in-radians is easier to view
        overlapPlot(h$Radians, a$Radians, xcenter = c("midnight"), linewidth=c(2,2), linecol=c("darkblue", "red"), olapcol="lightgray", main=title, ylab="Proportion of Activity", xlab="Time of Day", xlim=c(-7,9), ylim=c(0,0.2), xaxt = "n") #plot temporal overlap graph, with graph centred on midday and with x limit between the hours of 6pm and 8am, different coloured and textured line by species, and appropriate labels
        axis(side=1, at=tickpos, labels=axisvalues) #add x axis values as defined above
        grid(ny=NULL,nx=NA); abline(v=tickpos[seq(1,length(tickpos),2)],lty="dotted",col="lightgray") #add gridlines that line up with the new tick positions on the x axis
        legend('top', c("Human", SpeciesNames[s]), lty=c(1,2), lwd=c(2,2), col=c(4,2), horiz=TRUE, box.col="black", bg="white") #add legend at the top of the plot
        mysubtitle= t; mtext(side=3,line=0.3, mysubtitle) #add subtitle defining Human Abundance used for the graph
        dev.off()}else{ #else if less than ten tags for any species
          RESULTS[j+(s-1)*2,4]<-"NA" #set all values in the results matrix to NA
          RESULTS[j+(s-1)*2,5]<-"NA"
          RESULTS[j+(s-1)*2,6]<-"NA"}
    }
    write.csv(RESULTS, "Temporal Overlap Results-High and Low Abundance.csv") #save completed results matrix to directory
  }
}

HLtemporalOVERLAP(exifNIGHT,mdFULL, RAI) #run function




