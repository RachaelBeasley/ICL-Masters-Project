#SPATIAL DISPLACEMENT
#LINEAR MIXED MODELS
#DETERMINES THE RELATIONSHIP BETWEEN HUMAN RELATIVE ABUNDANCE INDICES AND MAMMAL RELATIVE ABUNDANCE INDICES 
#AND HOW THIS RELATIONSHIP VARIES BY HABITAT COVER AND TYPE
rm(list=ls())

##### LOAD DATA #####
setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/4) Temporal Displacement") #set working directory to retrieve data
RAIResults<-read.csv("RAI Results MD EXIF-Two Minutes.csv" , header=TRUE) #load Relative Abundance Indices Results
View(RAIResults) #visualise examine data
str(RAIResults)#see data structure

#HABITATS
setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/Data") #set working directory to retrieve data
Habitats<-read.csv("Site Habitat.csv", header=TRUE) #load data with Habitat Types and Covers by camera site
View(Habitats)
str(Habitats)

#SURVEY LENGTHS
#N.B. THE MAIN MONTH OF EACH CAMERA TRAP'S SURVEY HAS BEEN ADDED TO THE "SURVEY LENGTH.CSV" RESULTS FILE
#AS HAS THE DETECTION ZONE
setwd("C:/Users/beasl/Documents/Project/Project Data/FINAL CODES AND RESULTS/Final Codes/DATA AND CODES/2) Miscellaneous") #set working directory to retrieve results
SL<-read.csv("Survey Lengths-Month and Detection Zone Added.csv", header=TRUE) #load Survey Length data with month and detection zone manually added to it
View(SL)
SL<-SL[,-c(1,8,9)] #Remove surplus columns; reduces the number of redundant columns in the dataset

#MERGE DATASETS INTO A FINAL ONE FOR ANALYSIS
FINAL<-merge(RAIResults, Habitats, by.x="PlaceID", by.y="TagPlaceID") #Merge RAIs and Habitat Type/Cover
FINAL<-merge(FINAL, SL, by.x="PlaceID", by.y="PlaceID") #Merge Survey Lengths data to final dataset
View(FINAL)

#INSTALL AND LOAD PACKAGES
install.packages("lme4") #To run mixed linear models
library("lme4")
install.packages("predictmeans") #To see residual plots
library("predictmeans")

#### EXAMINE DATA ####
#PLOT DATA
#ORIGINAL UNTRANSFORMED DATA
plot(Fox~Human, data=FINAL) #no linear trend and clustered around small RAI values
plot(Hedgehog~Human, data=FINAL) #no linear trend and clustered around small RAI values

#PLOT LOGGED DATA
#TO SEE IF TRANSFORMING MAKES IT EASIER TO VISUALISE A RELATIONSHIP, ESPECIALLY GIVEN THE NUMBER OF SMALL RAIs VALUES
#1 ADDED TO RAIs TO ALLOW ZERO VALUE TO STILL BE PLOTTED
plot(log(Fox + 1) ~ log(Human + 1), data=FINAL) #less clustered and perhaps slightly positive
plot(log(Hedgehog + 1) ~ log(Human + 1), data=FINAL) #Lot of zeroes, with no clear pattern

#ALTHOUGH THERE ARE NO CLEAR TRENDS, LINEAR MODELS OFFER AN INITIAL AND SIMPLE PRELIMINARY ANALYSIS
#TO IDENTIFY TENTATIVE TRENDS IN THE DATA. FURTHER STUDIES WITH MORE TIME CAN DETERMINE
#IF MORE COMPLEX MODELS CAN FIT THIS DATA BETTER

#NUMBER OF ZEROES IN THE MODEL
length(which(FINAL$Human > 0)) #110/114 camera sites have human sightings
length(which(FINAL$Hedgehog > 0)) #64/114 camera sites have hedgehog sightings
length(which(FINAL$Fox > 0)) #109/114 camera sites have fox sightings

#HEDGEHOGS HAVE A LOT OF ZEROS, SUGGESTING THAT A HURDLE MODEL MAY BE USEFUL 
#I.E. SPLITTING THE MODEL INTO A BINOMIAL ABSENCE/PRESENCE OF HEDGEHOGS
#AND A SECOND LINEAR MODEL FOR JUST HEDGEHOG PRESENCE DATA.

#SEE DISTRIBUTION OF RANDOM EFFECTS
#CHECK IF THEY ARE WELL DISTRUBUTED BETWEEN THE CATEGORIES
table(FINAL$Month) #April: 47, May: 55, June: 12
table(FINAL$Detection.Zone) #Far: 16, Moderate: 38, Close: 36, Limited: 24
table(FINAL$HabitatCover) #Closed: 61, Open: 53
table(FINAL$HabitatType) #Amenity Grassland: 13, Artifical Structure: 17; Broadleaved Woodland: 26; Dense Scrub: 24; Hedgerows: 21; Tall Herbs: 13
#ALL SEEM REASONABLE AND MATCH WHAT WOULD BE EXPECTED FROM THE THE DATASET

####### BUILD FULL MODELS WITH INTERACTIONS ######
#AIM: TO SEE IF SPECIES ABUNDANCE AND HUMAN ABUNDANCES RELATIONSHIPS CHANGE BETWEEN HABITAT TYPES AND COVERS

#### HABITAT COVER AS INTERACTION #####

#### FOXES ####
#FULL MODEL
lm_foxHCI = lmer(log(Fox + 1) ~ log(Human + 1)*HabitatCover+Detection.Zone+(1|Month), data=FINAL)
summary(lm_foxHCI)
(0.0131/(0.0131+0.1764))*100 #Variance explained by month is 6.91%
residplot(lm_foxHCI) 
#Check Diagnostics from residplot
#Normal QQ plot: residuals lie on the diagonal, with some deviation in the upper right hand and lower left hand corners 
#Residual vs fitted plot: Scatter is fairly uniform like a starry night, though not many residuals towards the far right of plot
#Red line is quite sigmoidal, though the curve is only especially pronounced at the higher fitted values. 
#Rest of the plot is considered linear enough to assume homogeneity in the variances and that linear modelling is acceptable to identify a tentative trend
#CONCLUSION: ASSUMPTION OF NORMALITY IS MET AND SO IS HOMOGENEITY OF VARIANCES. ECOLOGICAL DATA CAN NOT BE
#EXPECTED TO MEET THE ASSUMPTIONS PERFECTLY, AND THE RESIDUALS ARE DISTRUBUTED REASONABLY CLOSE ENOUGH TO 
#THE IDEAL EXPECTATIONS TO GIVE CONFIDENCE IN THE RESULTS OF THIS MODEL. 

#### HEDGEHOGS ####
lm_hogHCI = lmer(log(Hedgehog + 1) ~ log(Human + 1)*HabitatCover+Detection.Zone+(1|Month), data=FINAL, )
summary(lm_hogHCI)
(0.0007034/(0.0007034+0.0709854))*100 #Variance explained by month is 0.98%
residplot(lm_hogHCI) 
#Check Diagnostics from residplot
#Normal QQ plot is very curved at its right hand side with a distinct outlier (point 14), though fairly straight before this section
#Residual vs Fitted plot is mostly concentrated in the bottom right hand corner, with poor scatter and evidence of funnelling
#Red line is fairly flat suggesting an overall linear relationship
#Model not deemed to have sufficient homogeneity of variances, and there is slight uncertainty about its normality

#TO TRY AND IMPROVE MODEL'S RESIDUALS
#IDENTIFY SITE 14 AS AN OUTLIER?
lm_hog = lm(log(Hedgehog + 1) ~ log(Human + 1)*HabitatCover+Detection.Zone, data=FINAL) #to see cook's distance remove random effect. as month explains only 0.98% of variance, unlikely to affect results dramatically.
plot(lm_hog)
#RESIDUALS VS LEVERAGE SHOWS SITE 14 IS NOT OUTSIDE THE 0.5 COOK'S DISTANCE. 
#DECIDED TO KEEP SITE 14 IN

#ALTERNATIVE ATTEMPT TO MEET ASSUMPTIONS OF THE MODEL
#The bias in the residuals may be due to the excessive zeroes, so model will be split into two (Hurdle Model):
#Binomial model with absence and presence of hedgehogs, 
#and a linear model showing the relationships between abundances with only those present 

#MODEL WITH ABSENCE AND PRESENCE
#GENERATE DATA AS A BINARY INPUT (I.E. 0 IF ABSENT, 1 IF PRESENT)
hog_binary<-rep(0, length(FINAL$Hedgehog)) #empty results vector for binary presence (1) and absence (0) data
for(i in 1:length(FINAL$Hedgehog)){ #for every value in the vector hedgehog RAIs
  if(FINAL$Hedgehog[i]>0){hog_binary[i]=1} #save a 1 in the results vector if it is greater than zero
  }

#BINOMIAL MODEL WITH ABSENCE AND PRESENCE OF HEDGEHOGS AS RESPONSE VARIABLE
#MODEL OF THE STRUCTURE IS SAME AS MODELS ABOVE
glm_hog_HC<-glmer(hog_binary~log(FINAL$Human + 1)*FINAL$HabitatCover+FINAL$Detection.Zone+(1|FINAL$Month), family="binomial")
summary(glm_hog_HC)
#MONTH REMOVED AS IT EXPLAINS NO VARIANCE, SO CAUSING SINGULARITY
glm_hog_HC1<-glm(hog_binary~log(FINAL$Human + 1)*FINAL$HabitatCover+FINAL$Detection.Zone, family="binomial")
summary(glm_hog_HC1) #Now significant-Limited Zones and (almost) moderate zones are less likely to detect hedgehogs
residplot(glm_hog_HC1) #Diagnostic plots are not too informative due to using a two level category variable

#MODEL WITH ONLY PRESENCE (ZERO TRUNCATED MODEL)
FINAL_with_hogs<-subset(FINAL, Hedgehog >0) #subset hedgehog RAI vector to remove all zero values
lm_hog_presentHCI = lmer(log(Hedgehog+1) ~ log(Human+1)*HabitatCover+Detection.Zone+(1|Month), data=FINAL_with_hogs)
residplot(lm_hog_presentHCI) 
#Check diagnostics with residplot
#Normal QQ plot does not differ much from using original model with all datapoints
#Residuals vs Fitted plot has less obvious funneling than original model with all datapoints, though residuals are now skewed towards the left and not the right
#Red line is less flat than original interaction model

#OVERALL CONCLUSION OF WHICH MODEL TO USE:
#USE FULL MODEL, AS LITTLE DIFFERENCE IN THE NORMALITY OF THEIR RESIDUALS AND RESIDUALS ARE MORE LINEARLY DISTRIBUTED
#SCATTER IS SLIGHTLY IMPROVED IN THE HURDLE MODEL VS THE ORIGINAL MODEL, BUT NOT GREATLY SO
#TO KEEP CONSISTENCY BETWEEN FOXES AND HEDGEHOGS MODELS, AND AS THE HURDLE OFFERS ONLY VERY SLIGHT IMPROVEMENT OVER
#THE ORIGINAL INTERACTION MODEL, THE ORIGINAL MODEL WITH ALL THE DATAPOINTS WILL BE USED
#RESULTS WILL BE INTERPRETED CAUTIOUSLY AS THE ASSUMPTION OF HOMOGENEOUS VARIANCES HAS NOT BEEN MET, TO REDUCE
#THE PRECISION OF THE RESULTS. THIS WILL ESPECIALLY AFFECT ESTIMATES STRADDLING SIGNIFICANCE (PVALUE=0.05)

##### HABITAT TYPE AS INTERACTION #####

#### FOXES ####
lm_foxHTI = lmer(log(Fox + 1) ~ log(Human + 1)*HabitatType+Detection.Zone+(1|Month), data=FINAL)
summary(lm_foxHTI)
(0.01374 /(0.01374 +0.17059))*100 #Variance explained by month is 7.45%
residplot(lm_foxHTI) 
#Check diagnostics with residplot
#Normal QQ Plot is fairly linear with some deviation at the ends of the plot.
#Fitted Residual Plot is fairly scattered around the centre of the plot, though residuals are not fully towards the far right and left of plot
#The red line is quite S-shaped, though for the purposes of examining a general trend, linear modelling deemed sufficient

#CONCLUSION: MODEL DEEMED TO HAVE SUFFICIENT NORMALITY AND FAIRL HOMOGENEOUS VARIANCES FOR LINEAR MODELLING
#ECOLOGICAL DATA CANNOT BE EXPECTED TO MEET THE ASSUMPTIONS PERFECTLY, AND THE RESIDUALS ARE DISTRUBUTED REASONABLY CLOSE ENOUGH TO 
#THE IDEAL EXPECTATIONS TO GIVE CONFIDENCE IN THE RESULTS OF THIS MODEL. 

#NEXT DETERMINE IF THE SLOPES BETWEEN HABITAT TYPES ARE SIGNIFICANTLY DIFFERENT FROM EACHOTHER, AND NOT JUST THE REFERENCE, AMENITY GRASSLANDS
#PAIRWISE COMPARISONS
install.packages("emmeans") #package that enables specific comparisons between estimate means and slopes
library("emmeans")

#examine gradients between logged relative abundance indices of humans and foxes in different habitat types
#and compare them pairwise to see if any are significantly different from eachother
emtrends(lm_foxHTI, pairwise ~ HabitatType, var = "log(Human+1)") 
#results typed up and placed in thesis

#### HEDGEHOGS ####
lm_hogHTI = lmer(log(Hedgehog + 1) ~ log(Human + 1)*HabitatType+Detection.Zone+(1|Month), data=FINAL)
summary(lm_hogHTI)
#MONTH EXPLAINS NO VARIANCE IN HEDGEHOGS RAI,SO CAUSING SINGULARITY
#MONTH REMOVED FROM THE MODEL
lm_hogHTI = lm(log(Hedgehog + 1) ~ log(Human + 1)*HabitatType+Detection.Zone, data=FINAL) #simple linear model used as no random effects
summary(lm_hogHTI)
residplot(lm_hogHTI)
#Check Diagnostic plots using residplot
#Normal QQ plot is very curved at its right hand side, though fairly straight before this section
#Residual Fitted is funnelled and concentrated towards the bottom and centre of the plot. Scatter could be more uniform. Not a lot of homogeneity in the variances. 
#Red line is curved, but could be approximated overall by a linear model
#As linear models are reasonably resilient to some deviation from normality, normality is deemed sufficient for the model
#However, model lacks sufficient homogeneity of variances

#ALTERNATIVE ATTEMPT TO MEET ASSUMPTIONS OF THE MODEL
#The bias in the residuals may be due to the excessive zeroes, so model will be split into two (Hurdle Model):
#Binomial model with absence and presence of hedgehogs, 
#and a linear model showing the relationships between abundances with only those present 

#MODEL WITH ABSENCE AND PRESENCE
#BINARY RESPONSE VARIABLE CALCULATED ABOVE UNDER HABITAT COVER (I.E. 0 IF ABSENT, 1 IF PRESENT)
#MODEL STRUCTURE REMAINS THE SAME AS ABOVE MODEL BUT WITH BINARY RESPONSE VARIABLE
glm_hog_HT<-glmer(hog_binary~log(FINAL$Human + 1)*FINAL$HabitatType+FINAL$Detection.Zone+(1|FINAL$Month), family="binomial")
summary(glm_hog_HT)
#MONTH REMOVED AGAIN AS IT EXPLAINS NO VARIANCE, SO CAUSING SINGULARITY
glm_hog_HT1<-glm(hog_binary~log(FINAL$Human + 1)*FINAL$HabitatType+FINAL$Detection.Zone, family="binomial")
summary(glm_hog_HT1) #Limited detection zone almost affecting hedgehog presence
residplot(glm_hog_HT1) #As using a two level categorical variable, plots are not very clear.

#MODEL WITH ONLY PRESENCE (ZERO TRUNCATED MODEL)
lm_hog_presentHTI = lm(log(Hedgehog + 1) ~ log(Human + 1)*HabitatType+Detection.Zone, data=FINAL_with_hogs)
residplot(lm_hog_presentHTI) 
#Check diagnostics and compare to original model with all datapoints
#Normal QQ plot: little difference in the plots. 
#Residuals vs Fitted plot: still evidence of funneling, but red line is slightly flater in the zero-truncated data model

#OVERALL CONCLUSION: ASSUMPTIONS OF NORMALITY AND HOMOGENIETY OF VARIANCES ARE NOT BETTER MET BY EITHER MODEL
#THOUGH ZERO-TRUNCATED HAS SLIGHTLY MORE LINEARLY DISTRIBUTED RESIDUALS, I DECIDED TO USE THE FULL MODEL FOR THE ANALYSIS.  
#THIS IS BECAUSE THE IMPROVEMENT BY THE ZERO-TRUNCATED MODEL IS ONLY SLIGHT AND BELIEVE IT IS BETTER TO KEEP AS MUCH CONSISTENCY
#BETWEEN MODELS OF FOXES AND HEDGEHOGS TO MAKE MORE DIRECT COMPARISONS.

#NEXT DETERMINE IF THE SLOPES BETWEEN HABITAT TYPES ARE SIGNIFICANTLY DIFFERENT FROM EACHOTHER, AND NOT JUST THE REFERENCE, AMENITY GRASSLANDS
#PAIRWISE COMPARISONS
#examine gradients between logged relative abundance indices of humans and hedgehogs in different habitat types
#and compare them pairwise to see if any are significantly different from eachother
emtrends(lm_hogHTI, pairwise ~ HabitatType, var = "log(Human+1)")
#results typed up and placed in thesis

#### PLOTTING GRAPHS ####
install.packages("sjPlot")
library("sjPlot")

#Log abundance terms before putting in model for sjPlot to plot logs
#As sjplot tries to backtransform values when logged after being put in the model
FINAL$HumanLog<-log(FINAL$Human+1)
FINAL$FoxLog<-log(FINAL$Fox+1)
FINAL$HedgehogLog<-log(FINAL$Hedgehog+1)

#REDEFINE MODELS WITH NEW LOGGED TERMS FOR SJPLOT
#HABITAT COVER
lm_foxHCI_Plot = lmer(FoxLog ~ HumanLog*HabitatCover+Detection.Zone+(1|Month), data=FINAL)
lm_hogHCI_Plot = lmer(HedgehogLog ~ HumanLog*HabitatCover+Detection.Zone+(1|Month), data=FINAL)

#HABITAT TYPE
lm_foxHTI_Plot = lmer(FoxLog ~ HumanLog*HabitatType+Detection.Zone+(1|Month), data=FINAL)
lm_hogHTI_Plot = lm(HedgehogLog ~ HumanLog*HabitatType+Detection.Zone, data=FINAL) #simple linear model used as no random effects

#IF THE ERROR BELOW OCCURS WHEN PLOTTING THE GRAPHS BELOW, TYPE "dev.off()"
#Error in .Call.graphics(C_palette2, .Call(C_palette2, NULL)) : invalid graphics state
#Plot Habitat Cover
plot_model(lm_foxHCI_Plot, type = "pred", terms = c("HumanLog", "HabitatCover"), axis.title = c("Human Relative Abundance Index", "Fox Relative Abundance Index"), title="Predicted Fox Relative Abundance Index")
plot_model(lm_hogHCI_Plot, type = "pred", terms = c("HumanLog", "HabitatCover"), axis.title = c("Human Relative Abundance Index", "Hedgehog Relative Abundance Index"), title="Predicted Hedgehog Relative Abundance Index")

#Habitat Type
plot_model(lm_foxHTI_Plot, type = "pred", terms = c("HumanLog", "HabitatType"), axis.title = c("Logged Human Relative Abundance Index", " Logged Fox Relative Abundance Index"),title="Predicted Fox Relative Abundance Index")
plot_model(lm_hogHTI_Plot, type = "pred",  terms = c("HumanLog", "HabitatType"), axis.title = c("Logged Human Relative Abundance Index", "Logged Hedgehog Relative Abundance Index"), title="Predicted Hedgehog Relative Abundance Index")
#Saved by exporting with height of 400



#### MISCELLANEOUS ####
#EXAMINE CORRELATION BETWEEN HUMANS AND DOGS
#Not logged as sufficiently large RAIs unlike nocturnal mammals
#Detection Zone included as may affect detection of dogs, which are lower to the ground than humans.
#Month deemed unlikely to have affected human sightings as there are limited green spaces in London
#and people regular commute or walk dogs on Hampstead Heath
DOG_HUMAN = lm(Dog ~ Human+Detection.Zone, data=FINAL)
summary(DOG_HUMAN)
#58% of the variation in dog abundance is affected by human abundance
#Although lower than expected, this can be explained how dogs go off path more frequently than humans
#meaning they are more likely to encounter a camera trap away from a path
#Also there are far more humans generally than dogs on the heath at any one time, which 
#reduces the correlation between dog and human abundances 
#i.e. not all humans walk dogs, and dogs are often off leash where they can trigger a camera without it also detecting a human


#DISTRIBUTION OF HABITAT COVER BY HIGH AND LOW HUMAN ABUNDANCE
Low<-subset(FINAL, Human<1) #Only sites with low human abundance
High<-subset(FINAL, Human>=1) #only sites with high human abundance
table(Low$HabitatCover) #28 Closed #6 Open ~82% closed
table(High$HabitatCover) #33 Closed #47 Open ~41% closed
#Not unexpectedly higher proportion of closed habitat in low human abundance sites
#and higher proportion of open habitats in high human abundance

table(Low$HabitatType) #Predominantly broadleaved woodland and dense scrub. No Amenity grasslends
table(High$HabitatType) #More evenly distributed, with only dense scrub having less than ten sightings

