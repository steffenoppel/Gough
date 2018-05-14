###################################################################################################
########   GOUGH BUNTING AND MOORHEN PVA TO CALCULATE CAPTIVE POPULATION SIZE   ################
###################################################################################################
# written by Steffen Oppel on 8 May 2018 (steffen.oppel@rspb.org.uk)
# adjusted by Antje Steinfurth to include better demographic data
# major update on 14 May 2018 to include stochastic PVA in loop and create output documents
# simulation taken from https://naes.unr.edu/shoemaker/teaching/NRES-470/PVA1_421.html

library(popbio)
library(doParallel)
library(foreach)
library(tidyverse)
library(markdown)
library(rmarkdown)
library(knitr)

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Gough\\Eradication\\Aviculture")



##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######
#################### SIMULATION OF POPULATION TRAJECTORY ACROSS RANGE OF PARAMETERS ########################################
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######


### SPECIFY RANGE OF PARAMETERS FOR DEMOGRAPHIC MODEL ###

pop.size<-seq(30,120,2)				      ### population size in individuals - ENTER THE RANGE WE MAY KEEP IN CAPTIVITY
Sa<-seq(0.8,0.9,0.01)					      ### survival of adult birds
Sj<-seq(0.5,0.7,0.01)					      ### survival of first year birds
F<-seq(0.80,0.90,0.01)			        ### fecundity = number of fledglings raised per pair per year, based on Cuthbert and Hilton 2004 and unpublished data 2009


### SPECIFY PARAMETERS FOR POPULATION VIABILITY ANALYSIS ###

nyears <- 50                        ### number of years over which simulations are run 
nreps <- 500                        ### number of stochastic simulations
K <- 1000                           ### Carrying capacity of Gough Island for Gough Bunting
captfail <- 5                       ### number of birds that die during captivity or are not released in reproductive state [effectively removed from population]

### SPECIFY PARAMETERS FOR STOCHASTIC COMPONENTS OF PVA ###

SD_lambda <- 0.10                   ### standard deviation of lambda - how variable population growth rate is from year to year
Catastrophe_prob <- 0.001           ### probability of a major disease, fire or other catastrophe occurring - 0.1% chance of major catastrophe
Catastrophe_severity <- 0.25        ### magnitude of catastrophe, in terms of size of population that will survive (25%)

acceptable_risk <- 0.1              ### risk that is acceptable to managers that species will go extinct after eradication


#################### CREATING THE POPULATION MATRIX ##################################################
## this is just to test the base function in the loop
## SIMPLE RUN TO TEST FOR LATER INCORPORATION IN LOOP

bird.matrix<-expression(
  0,  F*0.5,
  Sj, Sa)

bird.vr<-list(F=0.80,Sa=0.75, Sj=0.5)
A<-matrix(sapply(bird.matrix, eval,bird.vr , NULL), nrow=sqrt(length(bird.matrix)), byrow=TRUE)
projections<-pop.projection(A,n=c(50,100),iterations=50)
projections$lambda



##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######
#################### DEFINING FUNCTIONS FOR STOCHASTIC POPULATION VIABILITY ANALYSIS ########################################
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######


####  Density-dependent Ricker model to simulate population growth

Ricker <- function(prev_abund){       # this is a function for computing next-year abundance -- includes env stochasticity
  prev_abund * exp(log(rnorm(1,R_max,SD_lambda))*(1-(prev_abund/K)))
}



####  Stochastic population viability analysis function

PVAdemo <- function(nreps,nyears,Init_N,R_max,K,Catastrophe_prob,Catastrophe_severity){
  PopArray2 <- array(0,dim=c((nyears+1),nreps))
  
  ## start looping through replicates
  
  for(rep in 1:nreps){
    
    # set initial abundance
    PopArray2[1,rep] <- Init_N - rpois(1,captfail)    ### initial abundance minus captive mortality and releases due to illness
    
    ### loop through years
    for(y in 2:(nyears+1)){
      ### stochasticity and density dependence
      nextyear <- max(0,trunc(Ricker(PopArray2[y-1,rep])))    ### calculates abundance based on Ricker model, rounded to integer and set to a min of 0
      
      ### catastrophe
      if(runif(1)<Catastrophe_prob) nextyear <- nextyear*Catastrophe_severity
      PopArray2[y,rep] <- nextyear 
    }
  }
  
  return(PopArray2)
}



#### CALCULATE PROPORTION OF SIMULATIONS WHERE SPECIES GOES EXTINCT

Extinction_bysim <- function(simdata){
  sum(apply(Default,2,function(t)  min(t)==0))/ncol(simdata)
}




##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######
#################### LOOP OVER ALL COMBINATIONS OF DEMOGRAPHY AND RUN STOCHASTIC PVA ########################################
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######


### COMPREHENSIVE TABLE OF ALL COMBINATIONS OF DEMOGRAPHIC PARAMETERS ###

simul_in<-expand.grid(pop.size, Sa, Sj,F)
dim(simul_in)
names(simul_in)<-c('pop.size','Sa','Sj','F')
SIM_OUT<-data.frame()


### setup parallel backend to use 8 processors

cl<-makeCluster(8)
registerDoParallel(cl, cores=8)



### START PARALLEL LOOP AND RETURN LAMBDA AND PROPORTION OF SIMULATIONS WHERE SPECIES GOES EXTINCT

SIM_OUT<-foreach(s=c(1:dim(simul_in)[1]), .packages='popbio',.combine=rbind) %dopar% {


### CREATE LESLIE MATRIX WITH SUBSET OF VITAL RATES

bird.vr<-list(F=simul_in[s,4], Sa=simul_in[s,2],Sj=simul_in[s,3])
A<-matrix(sapply(bird.matrix, eval,bird.vr , NULL), nrow=sqrt(length(bird.matrix)), byrow=TRUE)
pop.size<-c((simul_in[s,1]/3),(simul_in[s,1]/3)*2)
projections<-pop.projection(A,n=pop.size,iterations=50)


### STOCHASTIC RICKER POPULATION MODEL

R_max <- projections$lambda       # Maximum rate of growth (max lambda)
Init_N <- simul_in[s,1]
Default <- PVAdemo(nreps,nyears,Init_N,R_max,K,Catastrophe_prob,Catastrophe_severity)


### CALCULATING EXTINCTION PROBABILITY

out<-simul_in[s,]
out$lambda<-projections$lambda
out$OUTCOME<-Extinction_bysim(Default)

return(out)
}
stopCluster(cl)







##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######
#################### SUMMARISE OUTPUT AND PLOT CAPTIVE POPULATION SIZE WITH RISK VALUE ######################################
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######

head(SIM_OUT)


#### SIZE OF MINIMUM NUMBER OF BIRDS TO BE HELD IN CAPTIVITY GIVEN THE ACCEPTABLE RISK ####

POPEST<- SIM_OUT %>% group_by(pop.size) %>%
  summarise(risk=max(OUTCOME)) %>%
  filter(risk<acceptable_risk) %>%
  summarise(CAPTIVE_POP=min(pop.size))
  


#### PLOT SIMULATION OUTPUT ####

  
ggplot(SIM_OUT, aes(x=pop.size, y=OUTCOME)) + geom_point() +
  geom_hline(aes(yintercept=acceptable_risk), color='red', size=1) +
  xlab("Captive population size") +
  ylab("Prob. of Gough Bunting extinct within 50 years") +
  annotate("text", x=75, y=0.65, label= sprintf("Captive population of %s birds needed", POPEST), size=5, colour= 'red') + 
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=16, color="black"), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=16, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())




write.table(SIM_OUT,"Gough_Bunting_PVA_simulation_output.csv", sep=",", row.names=FALSE)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PRODUCE OUTPUT REPORT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



### create HTML report for overall summary report
Sys.setenv(RSTUDIO_PANDOC="C:/Program Files (x86)/RStudio/bin/pandoc")

rmarkdown::render('S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Gough\\Eradication\\Aviculture\\CaptivePopSizePVA.Rmd',
                  output_file = "RSPB_Captive_GoughBunting_assessment.html",
                  output_dir = 'S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Gough\\Eradication\\Aviculture')



