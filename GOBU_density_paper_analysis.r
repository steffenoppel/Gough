### ##############################################
### # GOUGH BUNTING DENSITY ASSESSMENT ####
### ##############################################

# written by steffen.oppel@rspb.org.uk
# 10 May 2019
# updated 15 May 2019 - used gdistsamp data from 2018 and 2019 - all other approaches yielded nonsensical density data
# discarded dsm model
# discarded distsamp models of all data (2017 and before Aug 2018)


### LOAD THE REQUIRED PACKAGES 
library(tidyverse)
library(data.table)
library(readxl)
library(geosphere)
library(data.table)
library(lubridate)
library(sp)
library(raster)
library(sf)
library(unmarked)
filter<-dplyr::filter
select<-dplyr::select





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD LANDBIRD TRANSECT DATA FROM DATABASE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## RUN SCRIPT IN 32-bit R TO re-run after changes to database
system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\DATA\\Landbirds\\RODBC_GOBU_import_with_TRAL.R")), wait = TRUE, invisible = FALSE)

try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\Landbird_abundance"), silent=T)
#try(setwd("C:\\Users\\Gough Conservation\\Dropbox\\Gough Conservation team folder\\Gough Conservation databases"), silent=T)
load("RODBC_GOBU_TRAL_input_data.RData")

head(surv)
head(trans)
head(counts)
head(TRAL2016)
trans_orig<-trans

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD ALBATROSS BREEDING COUNT DATA FROM DATABASE TO CALCULATE 'MOUSE INDEX' (=TRAL breeding success)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### SO FAR THIS IS DATA FROM 2018 - need to think about whether we want long-term average or temporal matching


mouseindex2016 <- TRAL2016 %>% select(Colony, INCU, CHIC) %>%
  mutate(Success=(CHIC/INCU)) %>%
  mutate(year=2016)
mouseindex2017 <- TRAL2017 %>% select(Colony, INCU, CHIC) %>% ### no incubator counts in Green Hill North and False Peak in Feb 2017 due to low cloud
  mutate(Success=(CHIC/INCU)) %>%
  mutate(year=2017)
mouseindex2018 <- TRAL2018 %>% select(Colony, INCU, CHIC) %>%
  mutate(Success=(CHIC/INCU)) %>%
  mutate(year=2018)


### calculate average index over time (averaged breeding success, inverted to indicate where there are more mice)

mouseindex<-rbind(mouseindex2016, mouseindex2017, mouseindex2018) %>%
  group_by(Colony) %>%
  summarise(index=1-(mean(Success, na.rm=T)))




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# INSERT ZERO COUNTS FOR SURVEYS THAT DID HAPPEN BUT RECORDED NO BIRDS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
allcounts<-surv$LandbirdSurveyID
countsWithGOBUData<-unique(counts$LandbirdSurveyID[counts$Species=="GOBU"])
zeroGOBUCounts<-allcounts[!(allcounts %in% countsWithGOBUData)]
addZerosB<- data.frame(LandbirdSurveyID=zeroGOBUCounts, Species="GOBU",Age="Adult",Sex="Male",Detection="Heard",Behaviour="on_ground",Distance="d0_20m",N_birds=0)
counts<-rbind(counts, addZerosB)
head(counts)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DECIDING WHICH DATA TO USE AND HOW TO STRUCTURE THEM TEMPORALLY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### plot of the number of transect counts per week in every year
surv %>% mutate(month=month(Date), year=year(Date), week=week(Date), count=1) %>%
  group_by(year, month, Transect) %>%
  summarise(N_surveys=sum(count)) %>%
  arrange(year,month) %>%
  
  ggplot()+ geom_point(aes(x=month,y=Transect, col=as.factor(N_surveys))) + facet_wrap(~year, ncol=2)

### insufficient temporal replicates to use gdistsamp for most data - hence we use distsamp instead





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATING A MATRIX WITH ALL THE SITE METADATA (siteCov)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Simplify habitat into lowlands and uplands and calculate transect length
trans <- trans %>% mutate(Length = diag(distm(cbind(Long_start, Lat_start), cbind(Long_end, Lat_end), fun = distVincentyEllipsoid))) %>%
  mutate(Length = ifelse(is.na(Length), 100, Length)) %>%
  #mutate(Habitat=ifelse(Habitat_description %in% c("Sphagnum bogs","Moorland", "Wet heath"),"upland", "lowland")) %>%
  #mutate(Habitat_description=droplevels(Habitat_description)) %>%  
  select(Transect, Length, Habitat_description, Biome, Colony, IslandPart,Lat_start,Long_start) %>%
  mutate(MICE=mouseindex$index[match(Colony, mouseindex$Colony)]) %>%
  mutate(MICE=ifelse(is.na(MICE),1,MICE))       #### set all areas outside of TRAL colony areas to mouse index =1, indicating higher abundance of mice

### convert lat and long into UTM grid
trapSP<-SpatialPoints(trans[,8:7], proj4string=CRS("+proj=longlat +datum=WGS84"))
trapTransformed<-spTransform(trapSP, CRSobj =CRS("+proj=aeqd"))
trans[,8:7]<-trapTransformed@coords

head(trans)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# READING IN ISLAND HABITAT LAYER AND PROJECTING
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### HABITAT AREAS POTENTIALLY AVAILABLE TO GOBU
# https://remap-app.org/about
# Transects start and end point with associated habitat description was used as a training set 
# there was not enough information on the classes 'coastal tussock' and 'bog' so those were excluded
# water was identified via google aerial imagery (i.e. some steep areas were erroneously identified as water)
# predictors: reflectances using LANDSAT8 (NDVI, NDWI, BG, Blue, Green, red, near infrared) and elevation and slope using SRTM

### read in geotiff
## unique(trans$Habitat_description)
habitatmap<-raster("export.classification.tif")

### transform into spatial projection as for transects
habTransformed<-projectRaster(habitatmap,crs=CRS("+proj=aeqd"), method="ngb") 
plot(habTransformed)




### read in coastline data
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\DATA\\MAPS\\GIS layers\\Topography")
coastline <- st_read('Gough Island Shoreline.shp')



### create transect line simple feature for plotting on map
head(trans_orig)
transects<- trans_orig %>% select(Transect, Lat_start,Long_start,Lat_end,Long_end) %>%
  gather(key="Point", value="deg",-Transect) %>%
  mutate(Coordinate=ifelse(grepl("Long",Point)==T,"Longitude","Latitude")) %>%
  mutate(Point=ifelse(grepl("start",Point)==T,"start","end")) %>%
  spread(key=Coordinate, value=deg) %>%
  sf::st_as_sf(coords = c("Longitude","Latitude")) %>% 
  sf::st_set_crs(4326) %>%
  group_by(Transect) %>%
  summarise(Point=first(Point)) %>%
  st_cast("LINESTRING")








### #######################################################################
### #######################################################################
### GOUGH BUNTING density estimation USING HIERARCHICAL DISTANCE SAMPLING - FOR COUNTS IN SEPT 2018 and FEB 2019
### #######################################################################
### #######################################################################

library(unmarked)
head(surv)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATING A MATRIX WITH ALL THE SURVEY METADATA (obsCov)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### CREATE BLANK ARRAY FOR INPUT DATA BECAUSE NOT ALL TRANSECTS WERE SURVEYED 3 times
SURVEYDATA<-data.frame(Transect=rep(unique(surv$Transect),6),Year=rep(c(2018,2019,2018,2019,2018,2019), each=length(unique(surv$Transect))),
                       Count=rep(1:3, each=2*length(unique(surv$Transect))))


### REMOVE OLD SURVEYS and FORMAT DATE AND TIME 
input <- surv %>%
  mutate(month=month(Date), Year=year(Date)) %>%
  mutate(seqday = yday(Date)) %>%
  mutate(time=hour(Start_time)) %>%
  filter(Date>ymd("2018-08-15")) %>%
  select(LandbirdSurveyID,Transect,Year,seqday,time,Wind,Visibility) %>%
  arrange(Transect, Year, seqday)
head(input)


### assign sequential Survey Nr (1-3) TO ALLOW EFFICIENT CASTING OF DATA
input$Count<-NA
for (y in 2018:2019) {
  for (t in 1:50){
    x<-input %>% filter(Transect==t) %>% filter(Year==y) %>%
      arrange(seqday)
    if(dim(x)[1]>0){
      x$Count<-seq(1:dim(x)[1])
      input$Count[match(x$LandbirdSurveyID,input$LandbirdSurveyID)]<-x$Count}
  }
}

input<-input %>% filter(Count<4) %>% arrange(Transect, Year, Count) ### remove excessive surveys
dim(input)

### REMOVE TRANSECTS THAT WERE NEVER SURVEYED AND THEN MERGE BLANK MATRIX WITH DATA
SURVEYDATA<-SURVEYDATA %>% filter(Transect %in% unique(input$Transect)) %>%
  arrange(Transect, Year, Count)

obscovs<-merge(y=input, x=SURVEYDATA, by=c("Transect", "Year","Count"), all.x=T) %>% arrange(Transect,Year,Count)
dim(obscovs)
head(obscovs)  ### these are the observation covariates passed to the unmarkedFrame

## find missing surveys
obscovs %>% filter(is.na(Wind))



### SITE COVARIATE DATA FRAME for unmarked Frame class
### needs one value per transect for vegetation etc.
head(trans)
sitecovs <- rbind(trans,trans) %>%
  filter(Transect %in% unique(obscovs$Transect)) %>%
  mutate(Year=rep(c(2018,2019), each=nrow(trans))) %>%
  #mutate(Habitat_description=ifelse(Habitat_description=="Sphagnum bogs","Moorland",as.character(Habitat_description))) %>%   ### MODEL BREAKS IF WE DO THIS! combine moorland and sphagnum as they are the same on the aerial image classification
  arrange(Transect,Year)

head(sitecovs)
dim(sitecovs)*3 ### this should match with dim(obscovs)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BASIC SUMMARY OF RAW OBSERVATIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
counts %>% filter(Species=="GOBU") %>% 
  left_join(input, by="LandbirdSurveyID") %>%
  filter(Year>2017) %>%
  left_join(trans, by="Transect") %>%
  group_by(LandbirdSurveyID, Habitat_description) %>%
  summarise(N=sum(N_birds)) %>%
  ungroup() %>%
  group_by(Habitat_description) %>%
  summarise(N=mean(N, na.rm=T), vary=sd(N, na.rm=T))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATING A MATRIX WITH the count data in distance bins
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### select species and summarise observations within the same category
### DETECTION AND AGE AND SEX WERE DISCARDED HERE AS TOO FEW OBSERVATIONS PER GROUP

simplecounts<-counts %>% filter(Species=="GOBU") %>%
  group_by(LandbirdSurveyID, Distance) %>%
  summarise(N=sum(N_birds)) 
dim(simplecounts)
head(simplecounts)



### CAST THE MOLTEN DATA FRAME INTO MATRIX WITH 1 COLUMN PER COUNT

GOBU_y_prelim<- simplecounts  %>% filter(LandbirdSurveyID %in% obscovs$LandbirdSurveyID) %>%
  mutate(Transect=obscovs$Transect[match(LandbirdSurveyID,obscovs$LandbirdSurveyID)]) %>%
  mutate(Year=obscovs$Year[match(LandbirdSurveyID,obscovs$LandbirdSurveyID)]) %>%
  mutate(Count=obscovs$Count[match(LandbirdSurveyID,obscovs$LandbirdSurveyID)]) %>%
  select(LandbirdSurveyID,Transect,Year,Count,Distance,N) %>%
  spread(key=Distance, value=N, fill=0) %>%
  arrange(Transect,Year,Count)
head(GOBU_y_prelim)
summary(GOBU_y_prelim)


### ADD MISSING SURVEYS TO ENSURE THAT SIZE OF DATA FRAME IS CORRECT
dim(GOBU_y_prelim)
GOBU_y_prelim<-merge(obscovs[,1:4], GOBU_y_prelim, by=c("LandbirdSurveyID", "Transect", "Year","Count"), all.x=T) %>% arrange(Transect,Year,Count)
names(GOBU_y_prelim)[5:8]<-c('d1','d2','d3','d4')
head(GOBU_y_prelim)
dim(GOBU_y_prelim)
summary(GOBU_y_prelim)

GOBU_y_prelim %>% filter (is.na(LandbirdSurveyID))



### CAST DATA FOR UNMARKED INPUT

umfgds<-GOBU_y_prelim %>% select(-LandbirdSurveyID) %>%
  gather(key="Distance", value="N",-Transect, -Year, -Count) %>%
  mutate(CountDist=paste(Count,Distance, sep="_")) %>%                             ### for gdistsamp spread by this value
  select(Transect,Year,CountDist,N) %>%
  spread(key=CountDist, value=N) %>%
  arrange(Transect, Year)

umfgds %>% filter(is.na(`1_d1`))


### REMOVE TRANSECTS 39 and 47 from 2018 as they were not surveyed
umfgds<- umfgds %>% filter(!is.na(`1_d1`))%>%
  arrange(Transect, Year)

GOBU_y <- as.matrix(umfgds[,c(3:14)], dimnames=NULL)	                  ### convert to matrix which is required for 'unmarked', consider column order

sitecovs<- sitecovs %>% filter(!(Year==2018 & Transect %in% c(39,47)))%>%
  arrange(Transect, Year)

obscovs<- obscovs %>% filter(!(Year==2018 & Transect %in% c(39,47)))%>%
  arrange(Transect, Year)



### CHECK THAT DATA FRAMES ARE RIGHT SIZE
dim(GOBU_y)
dim(sitecovs)
dim(obscovs)


### CHECK FOR MISSING VALUES AND REPLACE THEM - THIS IS IRRELEVANT BECAUSE THE DATA INDICATE THAT THE SURVEYS DID NOT OCCUR
obscovs %>% filter(is.na(seqday))
obscovs %>% filter(is.na(Visibility))
obscovs$Visibility[is.na(obscovs$Visibility)]<-'good'
obscovs$Wind[is.na(obscovs$Wind)]<-'moderate'
obscovs$time[is.na(obscovs$time)]<-mean(obscovs$time, na.rm=T)
#sitecovs$Habitat_description<-factor(sitecovs$Habitat_description, levels=unique(sitecovs$Habitat_description))

### COMBINE RESPONSE AND OBSERVATION COVARIATES TO UNMARKED FRAME AND STANDARDIZE NUMERIC COVARIATES
GOBU_UMF <- unmarkedFrameGDS(y=GOBU_y,
                             siteCovs=sitecovs,
                             yearlySiteCovs=obscovs,
                             dist.breaks=c(0, 20, 100, 300,600),
                             survey='line',
                             numPrimary=3,
                             tlength = sitecovs$Length,
                             unitsIn='m')

### standardize site covariates
yearlySiteCovs(GOBU_UMF)[,c(5,6)] <- scale(yearlySiteCovs(GOBU_UMF)[,c(5,6)])
summary(GOBU_UMF)



### ##############################################
### ANALYSIS OF DATA IN GDISTSAMP
### ##############################################

### in gdistsamp, first formula is for abundance, second for availability, third for detection


### FIRST WE TEST FOR THE MOST PARSIMONIOUS MODEL FOR DETECTION
### some of these models do not converge because the three repeat surveys are done under identical conditions

null <- gdistsamp(~1, ~1, ~1, data = GOBU_UMF, keyfun="hazard", output="density", unitsOut="kmsq")
null.n <- gdistsamp(~1, ~1, ~1, data = GOBU_UMF, keyfun="halfnorm", output="density", unitsOut="kmsq")
null.e <- gdistsamp(~1, ~1, ~1, data = GOBU_UMF, keyfun="exp", output="density", unitsOut="kmsq")
null.wind <- gdistsamp(~1, ~Wind, ~1, data = GOBU_UMF, keyfun="hazard", output="density", unitsOut="kmsq")
null.e.wind <- gdistsamp(~1, ~1, ~Wind, data = GOBU_UMF, keyfun="exp", output="density", unitsOut="kmsq")
null.time <- gdistsamp(~1, ~1, ~time, data = GOBU_UMF, keyfun="hazard", output="density", unitsOut="kmsq")
null.e.time <- gdistsamp(~1, ~1, ~time, data = GOBU_UMF, keyfun="exp", output="density", unitsOut="kmsq")

detfl <- fitList(null,null.n,null.e,null.wind,null.e.wind, null.time, null.e.time)
detms <- modSel(detfl, nullmod="null")
detms



### THE NULL MODEL HAZARD RATE IS THE MOST SENSIBLE DETECTION MODEL
### using hazard rate with various predictors for density differences across time and space

habitat <- gdistsamp(~Habitat_description, ~1, ~1, data = GOBU_UMF, keyfun="hazard", output="density", unitsOut="kmsq")

### THESE MODELS ALL YIELDED NONSENSICAL DENSITIES
# biome <- gdistsamp(~Biome, ~1, ~1, data = GOBU_UMF, keyfun="hazard", output="density", unitsOut="kmsq")
# island <- gdistsamp(~IslandPart, ~1, ~1, data = GOBU_UMF, keyfun="hazard", output="density", unitsOut="kmsq")
# mice <- gdistsamp(~MICE, ~1, ~1, data = GOBU_UMF, keyfun="hazard", output="density", unitsOut="kmsq")
# year <- gdistsamp(~as.factor(Year), ~1, ~1, data = GOBU_UMF, keyfun="hazard", output="density", unitsOut="kmsq")


### AIC MODEL TABLE COMPARING THESE MODELS

fl <- fitList(null, habitat)
ms <- modSel(fl, nullmod="null")
ms


### PREDICTING GOBU DENSITY ACROSS ISLAND

## the top model yields a rather ludicrous density of >600 birds/sq km!!
# GOBURESULT <- predict(biome, type='lambda', newdat=sitecovs, appendData=T, SE=T)
# GOBUsummary<-GOBURESULT %>% group_by(Biome) %>%
#   summarise(Density=mean(Predicted), lcl=mean(lower), ucl=mean(upper))
# GOBUsummary


## habitat model is the only model with a sensible density estimate
GOBURESULT <- predict(habitat, type='lambda', newdat=sitecovs, appendData=T, SE=T)
GOBUsummary<-GOBURESULT %>% group_by(Habitat_description) %>%
  summarise(Density=mean(Predicted), lcl=mean(lower), ucl=mean(upper))
GOBUsummary

try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\Landbird_abundance"), silent=T)
fwrite(GOBUsummary,"GOBU_density_estimates.csv")


### ##############################################
### PLOTTING A GRAPH showing mean density across habitats
### ##############################################
pdf("GOBU_density_by_habitat.pdf", width=9, height=6)
ggplot(GOBUsummary,
       aes(x=Habitat_description, y=Density)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1) +
  xlab("") +
  ylab(expression(paste("Gough Bunting density (ind ", km^-2,")",sep=""))) +
  ylim(0,80) +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=16, color="black"), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=16, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PREDICT DENSITY ACROSS ISLAND AND CALCULATE TOTAL POPULATION SIZE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### convert raster to flat prediction data frame
preddata<-rasterToPoints(habitatmap)
preddata<-as.data.frame(preddata)
preddata$Habitat_description<-ifelse(
  preddata$export.classification==4,"Water", ifelse(
    preddata$export.classification==3,"Moorland", ifelse(
      preddata$export.classification==2,"Wet heath", "Fern bush")))
names(preddata)[3]<-"habitat"
preddata$area<-(24.8*30)
preddata <- preddata %>% filter(Habitat_description!="Water")

### calculate extent of habitat
preddata %>% group_by(Habitat_description) %>%
  summarise(extent=sum(area)/1000000)

# complement data with prediction of habitat-specific density
preddata<-preddata %>% left_join(GOBUsummary, by="Habitat_description")


#### SUMMARY TOTAL POPULATION 
GOBUout<-preddata %>% mutate(N=Density*(area/1000000),
                             lcl=lcl*(area/1000000),
                             ucl=ucl*(area/1000000)) %>%
  group_by(Habitat_description) %>%
  summarise(n=sum(N),lcl=sum(lcl), ucl=sum(ucl))
GOBUout
sum(GOBUout$lcl)
sum(GOBUout$n)
sum(GOBUout$ucl)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT DENSITY ACROSS ISLAND ON A MAP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### NOTE: to remove panel grid from plot we need to colour the lines white - element.blank() does NOT WORK: https://github.com/tidyverse/ggplot2/issues/2071

pdf("GOBU_density_habitat_map.pdf", width=9, height=6)

ggplot(preddata)+
  
  ### ADD PREDICTED BUNTING DENSITY
  geom_rect(aes(xmin=x-0.00015,ymin=y-0.00015,xmax=x+0.00015,ymax=y+0.00015, fill = Density),inherit.aes=FALSE) +
  scale_fill_gradient(name = expression(paste("Gough Bunting \n density (ind",km^-2,")",sep="")), low="white", high="red", guide = "colourbar", limits=c(0, 50))+

  ### ADD COASTLINE
  geom_sf(data=coastline, color = "darkolivegreen", lwd=1.5, fill=NA) +
  
  ### ADD SURVEY TRANSECTS
  geom_sf(data=transects, color = "black", lwd=1) +
  
  ### FORMAT AXES
  xlab("Longitude") +
  ylab("Latitude") +
  
  ### ADJUST FONT SIZE
  theme_classic()+
  theme(panel.background=element_rect(fill="white", colour="black"),
        plot.background = element_rect(fill = "white"),
        axis.text=element_text(size=14, color="black"),
        legend.text=element_text(size=14, color="black"),
        legend.title=element_text(size=18, color="black"),  
        axis.title=element_text(size=18), 
        panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(colour = "white"), 
        panel.border = element_blank())

dev.off()







### ##############################################
### GOF TEST OF TOP MODEL
### ##############################################

# test indicates that model is inadequate (p=0)!!

# Function adopted from Sillet et al. 2012
freeTuke <- function(fm) {
  observed <- getY(fm@data)
  expected <- fitted(fm)
  sum((sqrt(observed) - sqrt(expected))^2)
}

pb <- parboot(island, freeTuke, nsim=200, report=1)
pb











####################~~~~~~~~~~~~ ABANDONED 15 MAY 2019 ~~~~~~~~~~~~~~~~~~~~~~~~###############################











### #######################################################################
### #######################################################################
### GOUGH BUNTING density estimation USING DENSITY SURFACE MODEL
### #######################################################################
### #######################################################################
### see example analysis: http://distancesampling.org/R/vignettes/mexico-analysis.html
library(dsm)

### CREATE SEGMENT DATA

head(surv)

segment.data <- surv %>% mutate(month=month(Date), year=year(Date), seqday=yday(Date)) %>%
  mutate(time=hour(Start_time)) %>%
  left_join(trans, by='Transect') %>%
  mutate(habitat=ifelse(Habitat_description=='Fern bush',1,ifelse(Habitat_description=='Wet heath',2,ifelse(Habitat_description=='Moorland',3,ifelse(Habitat_description=='Coastal tussock',4,5))))) %>%
  select(LandbirdSurveyID,Transect,Length,habitat,Habitat_description,Biome,IslandPart,Lat_start,Long_start,MICE,year,month, seqday,time,Wind,Visibility) %>%
  dplyr::rename(Effort=Length, Sample.Label=LandbirdSurveyID,x=Long_start, y=Lat_start)
head(segment.data)

head(counts)

obs.data <- counts %>% filter(Species=="GOBU") %>% filter(N_birds>0) %>%
  mutate(distance=ifelse(Distance=='d0_20m',10,ifelse(Distance=='d21_100m',60,ifelse(Distance=='d101_300m',200,450)))) %>%
  dplyr::rename(size=N_birds, Sample.Label=LandbirdSurveyID) %>%
  mutate(object=seq_along(Sample.Label)) %>%
  select(object, Sample.Label, size,distance)

distdata <- obs.data %>% left_join(segment.data, by='Sample.Label') %>%
  select(object, Sample.Label, size,distance,Transect,Effort,Habitat_description,Biome,IslandPart,x,y,MICE,year,month, seqday,time,Wind,Visibility)
head(distdata)

### checking lowland data distribution
distdata %>% filter(Habitat_description=="Fern bush") %>% filter(size>0) %>%
  ggplot()+ geom_point(aes(x=month,y=Transect, col=as.factor(size))) + facet_wrap(~year, ncol=2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXPLORATORY DATA ANALYSIS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
table(counts$Distance)
par(mfrow=c(2,2))

# histograms
hist(obs.data$distance,main="",xlab="Distance (m)")
hist(obs.data$size,main="",xlab="Cluster size")

# plots of distance vs. cluster size
plot(obs.data$distance, obs.data$size, main="", xlab="Distance (m)",
     ylab="Group size", pch=19, cex=0.5, col=gray(0.7))

# lm fit
l.dat <- data.frame(distance=seq(0,500,len=1000))
lo <- lm(size~distance, data=obs.data)
lines(l.dat$distance, as.vector(predict(lo,l.dat)))





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ESTIMATING DETECTION FUNCTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(Distance)

### INITIAL COMPARISON WHAT FUNCTION WORKS BETTER

detfc.hr.null<-ds(distdata, key="hr", cutpoints=c(0, 20, 100, 300,600), adjustment=NULL)
detfc.hn.null<-ds(distdata, key="hn", cutpoints=c(0, 20, 100, 300,600), adjustment=NULL)


## Diagnostics
par(mfrow=c(1,2))
plot(detfc.hr.null, showpoints=FALSE, pl.den=0, lwd=2)
plot(detfc.hn.null, showpoints=FALSE, pl.den=0, lwd=2)
ddf.gof(detfc.hr.null$ddf)
ddf.gof(detfc.hn.null$ddf)
AIC(detfc.hr.null)
AIC(detfc.hn.null)
### hazard rate appears to fit data much better than half-normal!


### INCLUDE COVARIATES IN DETECTION FUNCTION
# wind, visibility, and day failed to fit any function
detfc.hr.time<-ds(distdata, cutpoints=c(0, 20, 100, 300,600), formula=~time, key="hr", adjustment=NULL)
detfc.hr.year<-ds(distdata, cutpoints=c(0, 20, 100, 300,600), formula=~year, key="hr", adjustment=NULL)
detfc.hr.month<-ds(distdata, cutpoints=c(0, 20, 100, 300,600), formula=~month, key="hr", adjustment=NULL)


### DETECTION FUNCTION VARIES BY YEAR (probably observer related)
ddf.gof(detfc.hr.year$ddf) ## yields a lower total sum of squares than hr-null model
AIC(detfc.hr.year)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FITTING DENSITY SURFACE MODEL USING THE SELECTED DETECTION FUNCTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
unique(segment.data$habitat)
head(segment.data)
dsm.xy <- dsm(count~s(x,y), detfc.hr.year, segment.data, obs.data, method="REML")
dsm.xy.biome <- dsm(count~s(x,y)+Biome, detfc.hr.year, segment.data, obs.data, method="REML")
dsm.xy.mice <- dsm(count~s(x,y)+MICE, detfc.hr.year, segment.data, obs.data, method="REML")
dsm.xy.island <- dsm(count~s(x,y)+IslandPart, detfc.hr.year, segment.data, obs.data, method="REML")
#dsm.xy.habitat.num <- dsm(count~s(x,y)+as.factor(habitat), detfc.hr.year, segment.data, obs.data, method="REML")
dsm.xy.habitat <- dsm(count~s(x,y)+Habitat_description, detfc.hr.year, segment.data, obs.data, method="REML")
dsm.biome <- dsm(count~Biome, detfc.hr.year, segment.data, obs.data, method="REML")
dsm.mice <- dsm(count~MICE, detfc.hr.year, segment.data, obs.data, method="REML")
dsm.island <- dsm(count~IslandPart, detfc.hr.year, segment.data, obs.data, method="REML")
dsm.habitat <- dsm(count~as.factor(habitat), detfc.hr.year, segment.data, obs.data, method="REML")
dsm.xy.tweedie <- dsm(count~s(x,y), detfc.hr.year, segment.data, obs.data, family=tw(), method="REML")

## including a correlation structure fails
#dsm(count~s(x,y), detfc.hr.year, segment.data, obs.data, engine="gamm",correlation=corSymm(form=~1| Transect), method="REML")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MODEL CHECKING
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### PLOT DENSITY SURFACE
par(mfrow=c(2,2))
vis.gam(dsm.xy, plot.type="contour", view=c("x","y"), type="response", contour.col="black", n.grid=100, main="Spatial")
vis.gam(dsm.xy.mice, plot.type="contour", view=c("x","y"), type="response", contour.col="black", n.grid=100, main="Mouse index")
vis.gam(dsm.xy.island, plot.type="contour", view=c("x","y"), type="response", contour.col="black", n.grid=100, main="Island part")
vis.gam(dsm.xy.habitat, plot.type="contour", view=c("x","y"), type="response", contour.col="black", n.grid=100, main="Habitat numeric")


### CHECK RESIDUALS 
gam.check(dsm.xy)
rqgam.check(dsm.xy.tweedie)


mod_results <- data.frame("Model name" = c("`dsm.xy`", "`dsm.xy.biome`", "`dsm.xy.tweedie`", "`dsm.xy.habitat`",
                                           "`dsm.xy.mice`", "`dsm.xy.island`",
                                           "`dsm.biome`", "`dsm.mice`", "`dsm.habitat`",
                                          "`dsm.island`"),
                          "Description" = c("Bivariate smooth of location, quasipoisson",
                                            "Bivariate smooth of location, upland vs. lowland, quasipoisson",
                                            "Bivariate smooth of location, Tweedie",
                                            "Bivariate smooth of location, habitat, quasipoisson",
                                            "Bivariate smooth of location, mouse index, quasipoisson",
                                            "Bivariate smooth of location, island part, quasipoisson",
                                            "upland vs. lowland, quasipoisson",
                                            "mouse index, quasipoisson",
                                            "habitat, quasipoisson",
                                            "island, quasipoisson"),
                          "Deviance explained" = c(unlist(lapply(list(dsm.xy,
                                                                      dsm.xy.biome,
                                                                      dsm.xy.tweedie,
                                                                      dsm.xy.habitat,
                                                                      dsm.xy.mice,
                                                                      dsm.xy.island,
                                                                      dsm.biome,
                                                                      dsm.mice,
                                                                      dsm.habitat,
                                                                      dsm.island),
                                                                 function(x){paste0(round(summary(x)$dev.expl*100,2),"%")}))))

mod_results




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROJECTING DENSITY ACROSS ISLAND USING HABITAT MAP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### convert raster to flat prediction data frame
preddata<-rasterToPoints(habTransformed)
preddata<-as.data.frame(preddata)
preddata$Habitat_description<-ifelse(
  preddata$export.classification==4,"Water", ifelse(
    preddata$export.classification==3,"Moorland", ifelse(
      preddata$export.classification==2,"Wet heath", "Fern bush")))
names(preddata)[3]<-"habitat"
preddata$area<-(24.8*30)
preddata <- preddata %>% filter(Habitat_description!="Water")

### this is the likely classification - need to reclassify transect data
# 1="Fern bush"
# 2="Wet heath"
# 3="Moorland"
# 4="Water"
#reclassify(habitatmap, rcl, include.lowest=FALSE, right=TRUE, ...)


# predict over a grid
preddata$GOBU<-predict(dsm.xy.habitat, preddata, preddata$area)

# calculate the predicted abundance of Gough Bunting over the entire island
sum(preddata$GOBU)
sum(preddata$area)/10000





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOTTING ABUNDANCE ACROSS ISLAND ON A MAP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pdf("GOBU_density_surface_model.pdf", width=10, height=6)
preddata %>% mutate(density=GOBU/(area/1000000)) %>%

ggplot()+geom_rect(aes(xmin=x-12.4,ymin=y-15,xmax=x+12.4,ymax=y+15, fill = density)) +
  scale_fill_gradient(name = 'Density Gough Buntings (n per sq km)', low="blue", high="red", guide = "colourbar", limits=c(0, 300))+
  #guides(fill=guide_legend(title="Density Gough Buntings (n per 0.1 ha)"))+

  #geom_point(aes(x=col_long,y=col_lat),size=3)+
  
  xlab("Longitude") +
  ylab("Latitude") +
  
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=14, color="black"),
        legend.text=element_text(size=14, color="black"),
        legend.title=element_text(size=18, color="black"),  
        axis.title=element_text(size=18), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())




dev.off()












### #######################################################################
### #######################################################################
### GOUGH BUNTING density estimation USING DISTANCE SAMPLING FROM ROYLE ET AL. 2004
### #######################################################################
### #######################################################################

# revised in September 2018 to provide index of abundance for sampling episode (Jan/Feb or Sept)
# reverted to gdistsamp because repeat visits were within 1 week
# used only habitat as covariate, no detection-level covariates were supported.

#GDISTSAMP APPROACH - FROM SEPT 2018 ONWARDS
library(unmarked)


### FORMAT THE SURVEY DATA in the required R*J MATRIX

GOBUds<-counts %>% filter(Species=="GOBU") %>%
  group_by(LandbirdSurveyID,Distance) %>%
  summarise(N=sum(N_birds,na.rm=T)) %>%
  spread(key=Distance, value=N, fill=0) %>%
  arrange(LandbirdSurveyID)
dim(GOBUds)
dim(segment.data)

head(GOBUds)

GOBU_y <- as.matrix(GOBUds[,c(2,4,3,5)], dimnames=NULL)  ## names of distance band columns is not in proper distance order






##### COMBINE RESPONSE AND OBSERVATION COVARIATES TO UNMARKED FRAME

GOBUdis<-unmarkedFrameDS(y=GOBU_y, siteCovs=segment.data, dist.breaks=c(0, 20, 100, 300,600), survey='line',tlength=segment.data$Effort, unitsIn='m')
summary(GOBUdis)



### ##############################################
### ANALYSIS OF DATA IN DISTSAMP
### ##############################################
names(segment.data)
### in distsamp, first formula is for detection, second for abundance
### BIOLOGICALLY PLAUSIBLE MODELS EVALUATED
### some of these models do not converge because the three repeat surveys are done under identical conditions

year <- distsamp(~year ~1, data = GOBUdis, keyfun="hazard", output="density", unitsOut="kmsq")
null <- distsamp(~1 ~1, data = GOBUdis, keyfun="hazard", output="density", unitsOut="kmsq")
month <- distsamp(~month ~1, data = GOBUdis, keyfun="hazard", output="density", unitsOut="kmsq")
day <- distsamp(~seqday ~1, data = GOBUdis, keyfun="hazard", output="density", unitsOut="kmsq")
time <- distsamp(~time ~1, data = GOBUdis, keyfun="hazard", output="density", unitsOut="kmsq")
year.n <- distsamp(~year ~1, data = GOBUdis, keyfun="halfnorm", output="density", unitsOut="kmsq")
null.n <- distsamp(~1 ~1, data = GOBUdis, keyfun="halfnorm", output="density", unitsOut="kmsq")
month.n <- distsamp(~month ~1, data = GOBUdis, keyfun="halfnorm", output="density", unitsOut="kmsq")
day.n <- distsamp(~seqday ~1, data = GOBUdis, keyfun="halfnorm", output="density", unitsOut="kmsq")
time.n <- distsamp(~time ~1, data = GOBUdis, keyfun="halfnorm", output="density", unitsOut="kmsq")
year.e <- distsamp(~year ~1, data = GOBUdis, keyfun="exp", output="density", unitsOut="kmsq")
null.e <- distsamp(~1 ~1, data = GOBUdis, keyfun="exp", output="density", unitsOut="kmsq")
month.e <- distsamp(~month ~1, data = GOBUdis, keyfun="exp", output="density", unitsOut="kmsq")
day.e <- distsamp(~seqday ~1, data = GOBUdis, keyfun="exp", output="density", unitsOut="kmsq")
time.e <- distsamp(~time ~1, data = GOBUdis, keyfun="exp", output="density", unitsOut="kmsq")

detfl <- fitList(null,year, month, day, time,null.n,year.n, month.n, day.n, time.n,null.e,year.e, month.e, day.e, time.e)
detms <- modSel(detfl, nullmod="null")
detms


habitat <- distsamp(~1 ~Habitat_description, data = GOBUdis, keyfun="halfnorm", output="density", unitsOut="kmsq")
biome <- distsamp(~1 ~Biome, data = GOBUdis, keyfun="halfnorm", output="density", unitsOut="kmsq")
island <- distsamp(~1 ~IslandPart, data = GOBUdis, keyfun="halfnorm", output="density", unitsOut="kmsq")
mice <- distsamp(~1 ~MICE, data = GOBUdis, keyfun="halfnorm", output="density", unitsOut="kmsq")
transect <- distsamp(~1 ~as.factor(Transect), data = GOBUdis, keyfun="halfnorm", output="density", unitsOut="kmsq")


fl <- fitList(null,habitat, biome, island, mice,transect)
ms <- modSel(fl, nullmod="null")
ms

GOBURESULT <- predict(habitat, type='state', newdat=segment.data, appendData=T, SE=T)
GOBUsummary<-GOBURESULT %>% group_by(Habitat_description) %>%
  summarise(Density=mean(Predicted), lcl=mean(lower), ucl=mean(upper))
GOBUsummary





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PREDICT AND PLOT DENSITY ACROSS ISLAND ON A MAP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### convert raster to flat prediction data frame
preddata<-rasterToPoints(habTransformed)
preddata<-as.data.frame(preddata)
preddata$Habitat_description<-ifelse(
  preddata$export.classification==4,"Water", ifelse(
    preddata$export.classification==3,"Moorland", ifelse(
      preddata$export.classification==2,"Wet heath", "Fern bush")))
names(preddata)[3]<-"habitat"
preddata$area<-(24.8*30)
preddata <- preddata %>% filter(Habitat_description!="Water")

# predict over a grid
preddata<-preddata %>% left_join(GOBUsummary, by="Habitat_description")



pdf("GOBU_density_habitat_map.pdf", width=10, height=6)

  
  ggplot(preddata)+geom_rect(aes(xmin=x-12.4,ymin=y-15,xmax=x+12.4,ymax=y+15, fill = Density)) +
  scale_fill_gradient(name = 'Density Gough Buntings \n (ind per sq km)', low="blue", high="red", guide = "colourbar", limits=c(0, 150))+
  #guides(fill=guide_legend(title="Density Gough Buntings (n per 0.1 ha)"))+
  
  #geom_point(aes(x=col_long,y=col_lat),size=3)+
  
  xlab("Longitude") +
  ylab("Latitude") +
  
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=14, color="black"),
        legend.text=element_text(size=14, color="black"),
        legend.title=element_text(size=18, color="black"),  
        axis.title=element_text(size=18), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())




dev.off()


#### SUMMARY TOTAL POPULATION 
GOBUout<-preddata %>% mutate(N=Density*(area/1000000),
                             lcl=lcl*(area/1000000),
                             ucl=ucl*(area/1000000))
sum(GOBUout$lcl)
sum(GOBUout$N)
sum(GOBUout$ucl)


### ##############################################
### PLOTTING A GRAPH showing mean density across habitats
### ##############################################
pdf("GOBU_density_by_habitat_distsamp.pdf")
ggplot(GOBUsummary,
       aes(x=Habitat_description, y=Density)) +   ##facet_wrap(~year, ncol=1) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1) +
  xlab("") +
  ylab("Density of Gough buntings per sq km") +
  ylim(0,150) +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=16, color="black"), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=16, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

dev.off()
















