### ##############################################
### # GOUGH BUNTING NEST SURVIVAL ####
### ##############################################

# written by steffen.oppel@rspb.org.uk on 13 May 2019
# based on old code of Montserrat Oriole, St Helena Plover, and Aquatic Warbler
# Shaffer's logistic exposure method: https://rpubs.com/bbolker/logregexp



### LOAD THE REQUIRED PACKAGES 
library(tidyverse)
library(data.table)
library(geosphere)
library(lubridate)
library(sp)
filter<-dplyr::filter
select<-dplyr::select




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD NEST DATA FROM DATABASE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## RUN SCRIPT IN 32-bit R TO re-run after changes to database
system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\DATA\\Breeding_Database\\RODBC_nest_import.R")), wait = TRUE, invisible = FALSE)

try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\DATA\\Breeding_Database"), silent=T)
#try(setwd("C:\\Users\\Gough Conservation\\Dropbox\\Gough Conservation team folder\\Gough Conservation databases"), silent=T)
load("GOUGH_nest_visit_data.RData")

head(nests)
head(visits)


nests<-nests %>% filter(Species=="GOBU") %>% filter(Year>2005)
visits<-visits %>% filter(NestID %in% nests$NestID)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SIMPLE SUMMARIES FOR PAPER
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nests %>% mutate(count=1) %>% group_by(Year) %>%
  summarise(n=sum(count), fledged=sum(SUCCESS), succ=mean(SUCCESS))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FORMAT NEST VISIT DATA AND IDENTIFY FIRST AND LAST VISIT AND DATE LAST ALIVE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

alivedates<- visits %>% group_by(NestID, Status) %>%
  summarise(first=min(Date), last=max(Date)) %>%
  filter(Status=='Alive') %>%
  select(NestID,first, last) %>%
  rename(FirstAlive=first,LastAlive=last)

pastdates<- visits %>% group_by(NestID, Status) %>%
  summarise(first=min(Date)) %>%
  filter(Status!='Alive')%>%
  rename(LastDate=first)


### CALCULATE EXPOSURE DAYS ASSUMING NEST TERMINATION HALF-WAY BETWEEN LAST VISITS

nests<- nests %>% #
  left_join(alivedates, by='NestID') %>%
  left_join(pastdates, by='NestID') %>%
  mutate(LastDate=if_else(is.na(LastDate),DateLastChecked,LastDate)) %>%
  mutate(LastAlive=if_else(is.na(LastAlive),DateLastAlive,LastAlive)) %>%
  select(NestID,Nest_label,Year,Colony,Latitude,Longitude,DateFound,SUCCESS, FirstAlive, LastAlive, LastDate,Status) %>%
  mutate(exposure=as.integer((LastDate-0.5*(LastDate-LastAlive))-FirstAlive)) %>%
  mutate(DayInit=yday(FirstAlive))



### CALCULATE INITIATION DAY AS DAYS SINCE FIRST NEST DAY
firstnests<- nests %>% group_by(Year) %>%
  summarise(init=min(FirstAlive)-days(1))
nests<- nests %>% 
  mutate(DayInit=as.integer(FirstAlive-firstnests$init[match(Year, firstnests$Year)]))


### ASSIGN PAIR ID BASED ON NEST LABEL
head(nests)
unique(nests$Nest_label)

nests<- nests %>% 
  mutate(Pair=ifelse(Year==2009, NestID,substr(Nest_label,1,21))) %>%
  mutate(Pair=as.integer(as.factor(Pair)))



###### FIND DATA ERRORS AND FIX IN DATABASE ###

nests %>% filter(is.na(LastDate))
pastdates %>% filter(NestID==17214)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE MIN NESTING INTERVAL LENGTH
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## min nest exposure time is 37 days


nests %>% mutate(minAlive=as.integer(LastAlive-FirstAlive)) %>%
  group_by(Year) %>%
  summarise(length=max(minAlive)/(3600*24)) ## not sure why this is reported in seconds...



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPECIFYING LOGISTIC EXPOSURE LINK FUNCTION (Shaffer 2004)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define logistic exposure family
# from https://rpubs.com/bbolker/logregexp

logexp <- function(exposure = 1)
{
  linkfun <- function(mu) qlogis(mu^(1/exposure))
  ## FIXME: is there some trick we can play here to allow
  ##   evaluation in the context of the 'data' argument?
  linkinv <- function(eta)  plogis(eta)^exposure
  logit_mu_eta <- function(eta) {
    ifelse(abs(eta)>30,.Machine$double.eps,
           exp(eta)/(1+exp(eta))^2)
    ## OR .Call(stats:::C_logit_mu_eta, eta, PACKAGE = "stats")
  }
  mu.eta <- function(eta) {       
    exposure * plogis(eta)^(exposure-1) *
      logit_mu_eta(eta)
  }
  valideta <- function(eta) TRUE
  link <- paste("logexp(", deparse(substitute(exposure)), ")",
                sep="")
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta, 
                 name = link),
            class = "link-glm")
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FITTING NEST SURVIVAL MODELS TO TEST FOR YEAR OR COLONY EFFECT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(lme4)

m0<-glmer(SUCCESS~1+(1|Pair), data=nests,family=binomial(link=logexp(exposure=nests$exposure)))
m1<-glmer(SUCCESS~Colony+(1|Pair), data=nests,family=binomial(link=logexp(exposure=nests$exposure)))
m2<-glmer(SUCCESS~as.factor(Year)+(1|Pair), data=nests,family=binomial(link=logexp(exposure=nests$exposure)))

## LIKELIHOOD RATIO TEST
anova(m0,m1)  ## no evidence for difference between Gonydale and Tafelkop
anova(m0,m2)  ## no evidence for difference between 2009 and 2018





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SUMMARISE NEST SURVIVAL FOR REPORTING MEAN AND CONFIDENCE INTERVAL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GOBU_nest_success<-data.frame()
for (s in unique(nests$Colony)){
  sitenests<-nests %>% filter(Colony==s)
  for (y in unique(nests$Year)){
    yearnests<-sitenests %>% filter(Year==y) %>% mutate(Fate=ifelse(SUCCESS==1,0,1))
    out<-data.frame(Colony=s,
                    Year=y,
                    n_nests=dim(yearnests)[1],
                    succ=mean(yearnests$SUCCESS, na.rm=T),
                    Mayfield_dsr=(1-(sum(yearnests$Fate)/sum(yearnests$exposure))))
    GOBU_nest_success<-rbind(GOBU_nest_success,out)
    }			# end year loop
  }			# end site loop


GOBU_nest_success$nest_survival<-GOBU_nest_success$Mayfield_dsr^37


## add confidence intervals from a simple GLM
nests$Fate<-ifelse(nests$SUCCESS==1,0,1)
m0<-glm(cbind(Fate,exposure)~Colony+Year, fam=binomial, data=nests, na.action=na.omit)
GOBU_nest_success$predicted<-predict(m0, newdata=GOBU_nest_success, type='link', se=T)$fit
GOBU_nest_success$se<-predict(m0, newdata=GOBU_nest_success, type='link', se=T)$se.fit

ilink <- family(m0)$linkinv
GOBU_nest_success<-GOBU_nest_success %>% mutate(lcl = (1-(ilink(predicted + (1.96*se))))^37, ucl = (1-(ilink(predicted - (1.96*se))))^37) %>%
  mutate(mean_nest_survival=(1-(ilink(predicted)))^37) %>%
  filter(n_nests>0) %>%
  select(Colony,Year,n_nests,succ,Mayfield_dsr,nest_survival,mean_nest_survival,lcl,ucl)

write.table(GOBU_nest_success, "GOBU_nest_survival_summary.csv",row.names=F, sep=',')


