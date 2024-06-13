
### data processing 
library(reshape2)
library(plyr)
library(dplyr)
library(tidyverse)

#plots
library(ggplot2)
library(patchwork)
library(sjPlot)
#library(kableExtra)
library(ggthemes)
library(ggpubr)

#models

library(lmerTest)
library(fitdistrplus)
library(AER)
library(easystats)
library(energy)
library(raster)


# load spatial files 

EPSG4326<-CRS("+proj=longlat +datum=WGS84 +no_defs")

strata<-shapefile("C:/D1MPA trimestral/ASAM2022 updated strata/Strata_6932.shp")

plot(strata)

all.counts<-read.csv("C:/D1MPA trimestral/AllCounts_V_4_1.csv")

# first subset data within strata
head(all.counts)

pcounts<-SpatialPointsDataFrame(all.counts[4:5],all.counts[1:14],proj4string = EPSG4326)

pcounts_6932<-sp::spTransform(pcounts,CRS("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))

pcounts.strata<-raster::intersect(pcounts_6932,strata)

pcs<-data.frame(pcounts.strata)
head(pcs)

nests<-subset(pcs,count_type=="nests")

### some colonies had multiple counts over the same season
### this summarizes the count with the maximum nests


countsN<-ddply(nests, c("site_id","common_name","longitude_epsg_4326","latitude_epsg_4326","ID"), summarise,
               ncounts=length(penguin_count),
               interval=(max(season_starting)-min(season_starting)))

head(countsN)
write.csv(countsN,"countsN.csv")

ggplot(countsN,aes(ncounts))+geom_histogram()
ggplot(countsN,aes(interval))+geom_histogram()


nestM<-merge(nests,countsN) # identify number of counts for each population by merging

head(nestM)

poisson.mtest(nestM$penguin_count[nestM$ncounts>5 & nestM$penguin_count>0],R=199)  ##test for poisson distribution

nestm3<-subset(nestM,ncounts>5 & penguin_count>0) 


# --------chinstrap penguin glmm---------------- 

chp<-subset(nestm3,common_name=="chinstrap penguin")

summary(as.factor(subset(nestm3,site_id=="PING")$site_name))

summary(as.factor(chp$ID))

head(chp)

chp<-subset(chp, season_starting>1979,select=c("site_id","longitude_epsg_4326","latitude_epsg_4326",
                                               "ID","season_starting","penguin_count","accuracy"))
summary(as.factor(chp$site_id))

summary(subset(chp,site_id=="PCHA" |site_id=="WATE")$penguin_count)
summary(subset(chp,site_id=="SELV" )$penguin_count)
summary(subset(chp,site_id=="USEF" )$penguin_count)
summary(as.factor(subset(chp,site_id!="PING")$accuracy))
summary(as.factor(chp$accuracy))

chp<-subset(chp,site_id!="CUVE" & site_id!="PETE") # only 2 counts
chp<-subset(chp,site_id!="PCHA" & site_id!="WATE") # very small colonies (less than 100 nests)

ping<-subset(chp,site_id=="PING" & season_starting>2000)

ggplot(ping,aes(season_starting,penguin_count))+
  geom_smooth(method="gam",formula=y~s(x,k=3),se=F)+
  geom_point()+
  theme_bw()+xlab("breeding seasons")+ylab("penguin counts")+
  ggtitle(label="Pinguino island chinstrap penguins")


ggplot(chp,aes(season_starting, penguin_count))+geom_point()+facet_wrap(site_id~.)


chp.ac<-subset(chp,accuracy<3)

chp.ac<-subset(chp.ac,site_id!="RENI" & site_id!="RUGG" &
                 site_id!="WEEK") 

# remove outliers based on posterior diagnostics of the GLMM

outliers <- c(217, 307)

chp.ac <- chp.ac[-outliers, ]


ggplot(chp.ac,aes(season_starting, penguin_count))+geom_point()+facet_wrap(site_id~.)


poisson.mtest(chp.ac$penguin_count,R=199)  ##test for poisson distribution

# GLMM

chp.mm<-glmer(penguin_count~scale(season_starting)+(scale(season_starting)|site_id),
              family="poisson",data=chp.ac) 


chp.m<-glm(penguin_count~scale(season_starting),family="poisson",data=chp.ac)

summary(chp.mm)# significantly decreasing
anova(chp.mm,chp.m)# significant random effect

# variance explained by random effect 

chp.mm

(2.4108+0.9644)/((2.4108+0.9644)+(5.5098+0.4797))

# (random intercept + random slope)/ ((random intercept + random slope)+(fixed intercept+fixed slope))

pred.chp<-predict(chp.mm,newdata=chp.ac)

plot((chp.ac$penguin_count),exp(pred.chp))


cor(chp.ac$penguin_count,exp(pred.chp)) 

 # the model can produce fair predictions 

### ---------gentoo penguins glmm---------
gep<-subset(nestm3,common_name=="gentoo penguin")
gep<-subset(gep, season_starting>1979,select=c("site_id","longitude_epsg_4326","latitude_epsg_4326",
                                               "ID","season_starting","penguin_count","accuracy"))
summary(as.factor(gep$site_id))
gep.ac<-subset(gep,accuracy<3)

ggplot(gep.ac,aes(season_starting, penguin_count))+geom_point()+facet_wrap(site_id~.)

gep.ac<-subset(gep.ac,site_id!="BENE" & site_id!="HUMP" &
                 site_id!="ROBE" & site_id!="SPTR") 

poisson.mtest(gep.ac$penguin_count,R=199)  ##test for poisson distribution

gep.mm<-glmer(penguin_count~scale(season_starting)+(scale(season_starting)|site_id),
              family="poisson",data=gep.ac) 
gep.m<-glm(penguin_count~scale(season_starting),family="poisson",data=gep.ac)

summary(gep.mm)# significantly increasing
anova(gep.mm,gep.m)# significant random effect


# variance explained by random effect 

gep.mm

(1.4112+0.6347)/((1.4112+0.6347)+(6.5691+0.3839))

# (random intercept + random slope)/ ((random intercept + random slope)+(fixed intercept+fixed slope))

pred.gep<-predict(gep.mm,newdata=gep.ac)

plot((gep.ac$penguin_count),exp(pred.gep) )

cor(gep.ac$penguin_count,exp(pred.gep))

# the model can produce fair predictions 

###----adelie penguin ---------


adp<-subset(nestm3,common_name=="adelie penguin")
adp<-subset(adp, season_starting>1979,select=c("site_id","longitude_epsg_4326","latitude_epsg_4326",
                                               "ID","season_starting","penguin_count","accuracy"))
summary(as.factor(adp$site_id))

adp.ac<-subset(adp,accuracy<3)

adp.ac<-subset(adp.ac,site_id!="HOPE" & site_id!="DETA"&
                 site_id!="PLEN" & site_id!="CORM"& 
                 site_id!="CHIS"& site_id!="TORG"&
                 site_id!="TURR"& site_id!="BONG"&
                site_id!="PCHA") 
                                 # some colonies with less than 5 counts were still included
                                    # CORM, CHIS and TORG have only data in the 1980s and
                                    # BONG and TURR only in the 2010s

ggplot(adp.ac,aes(season_starting, penguin_count))+geom_point()+facet_wrap(site_id~.)


poisson.mtest(adp.ac$penguin_count,R=199)  ##test for poisson distribution

adp.mm<-glmer(penguin_count~scale(season_starting)+(scale(season_starting)|site_id),
              family="poisson",data=adp.ac) 
adp.m<-glm(penguin_count~scale(season_starting),family="poisson",data=adp.ac)

summary(adp.mm)# significantly decreasing
anova(adp.mm,adp.m)# significant random effect

adp.mm

(1.0366+0.4218)/((1.0366+0.4218)+(7.4838+0.6784))

pred.adp<-predict(adp.mm,newdata=adp.ac)

plot((adp.ac$penguin_count),exp(pred.adp) )

cor(adp.ac$penguin_count,exp(pred.adp))

# the model results in accurate predictions

###  -----------plots ------------



(plot_model(chp.mm,type="emm",terms=c("season_starting"),pred.type = "re",ci.lvl=0.3,
            show.data=T,se=T,transform="exp",grid=T,vcov.fun = std)+
    theme_bw()+xlab("breeding season")+ylab("number of nests")+
    ggtitle(label="a. Chinstrap penguin colony size"))/
  
  
  
  (plot_model(adp.mm,type="emm",terms=c("season_starting"),pred.type = "re",ci.lvl=0.3,
              show.data=T,se=T,transform="exp",grid=T,vcov.fun = std)+
     theme_bw()+xlab("breeding season")+ylab("number of nests")+
     ggtitle(label="b. Adelie penguin colony size"))/
  
  
  (plot_model(gep.mm,type="emm",terms=c("season_starting"),pred.type = "re",ci.lvl=0.3,
              show.data=T,se=T,transform="exp",grid=T,vcov.fun = std)+
     theme_bw()+xlab("breeding season")+ylab("number of nests")+
     ggtitle(label="e. Gentoo penguin colony size"))


# random effect plot

plot_model(chp.mm,type="re",grid=F,sort.est="sort.all")[1]

#### --------- models diagnostics for outliers --------------


check_outliers(chp.mm) # two outliers. Get back and do it again
check_outliers(adp.mm) # no outliers
check_outliers(gep.mm) # no outliers

#### random effects

# in a GLMM the slope random effect is the difference from the mean slope
# so lets calculate the difference to know whether the tendency of the pop is decreasing or increasing

chp.re<-data.frame(ranef(chp.mm),species=c("CHP"))

summary(chp.mm)
chp.re$trend<-chp.re$condval-0.4797



adp.re<-data.frame(ranef(adp.mm),species=c("ADP"))
summary(adp.mm)
adp.re$trend<-adp.re$condval-0.6784


gep.re<-data.frame(ranef(gep.mm),species=c("GEP"))
summary(gep.mm)
gep.re$trend<-gep.re$condval+0.3839



re<-rbind(chp.re,
          adp.re,
          gep.re)



re<-subset(re,term=="scale(season_starting)")
redf<-data.frame(site_id=re$grp,slope=re$trend,sd=re$condsd,species=re$species)


coords<-ddply(nests, c("site_id"), summarise,
              lon=mean(longitude_epsg_4326),
              lat=mean(latitude_epsg_4326))

ranef<-merge(redf,coords)

head(ranef)


ggplot(ranef,aes(lat,slope,colour=species))+geom_point()
# latitudinal trend is clear

write.csv(ranef,"randomeffect.csv")
### full plot

(plot_model(chp.mm,type="emm",terms=c("season_starting"),pred.type = "re",ci.lvl=0.3,
            show.data=T,se=T,transform="exp",grid=T,vcov.fun = std)+
    theme_bw()+xlab("breeding season")+ylab("number of nests")+
    ggtitle(label="a. Chinstrap penguin colony size")+
  
  
ggplot(subset(ranef,species=="CHP"),aes(lat,slope))+
  geom_hline(yintercept = 0,linetype="dashed")+
geom_smooth(method="lm",se=F)+
  geom_errorbar(aes(ymin=slope-sd,ymax=slope+sd))+
  geom_point()+
  theme_bw()+xlab("Colony latitude")+ylab("Colony trend")+
  ggtitle(label="b. Chinstrap penguin latitudinal trend"))/

  


(plot_model(gep.mm,type="emm",terms=c("season_starting"),pred.type = "re",ci.lvl=0.3,
            show.data=T,se=T,transform="exp",grid=T,vcov.fun = std)+
    theme_bw()+xlab("breeding season")+ylab("number of nests")+
    ggtitle(label="c. Gentoo penguin colony size")+


ggplot(subset(ranef,species=="GEP"),aes(lat,slope))+
  geom_hline(yintercept = 0,linetype="dashed")+
  geom_smooth(method="lm",se=F)+
  geom_errorbar(aes(ymin=slope-sd,ymax=slope+sd))+
  geom_point()+
  theme_bw()+xlab("Colony latitude")+ylab("Colony trend")+
  ggtitle(label="d. Gentoo penguin latitudinal trend"))/


  (plot_model(adp.mm,type="emm",terms=c("season_starting"),pred.type = "re",ci.lvl=0.3,
              show.data=T,se=T,transform="exp",grid=T,vcov.fun = std)+
     theme_bw()+xlab("breeding season")+ylab("number of nests")+
     ggtitle(label="e. Adelie penguin colony size")+
  

ggplot(subset(ranef,species=="ADP"),aes(lat,slope))+
  geom_hline(yintercept = 0,linetype="dashed")+
  geom_smooth(method="lm",se=F)+
  geom_errorbar(aes(ymin=slope-sd,ymax=slope+sd))+
  geom_point()+
  theme_bw()+xlab("Colony latitude")+ylab("Colony trend")+
  ggtitle(label="f. Adelie penguin latitudinal trend"))  
  

### -------- lagged comparison How have changes occurred through 5-year periods?
library(lubridate)
library(tidyr)
library(tidyquant)
library(dplyr)
library(broom)
library(purrr)
library(stringr)
library(knitr)
library(timetk)

head(chp.ac)

summary(as.factor(chp.ac$site_id))



### ------Chinstrap Penguin-------------


year <- data.frame(season_starting = 1980:2022)

# Merge to include all years in the range
chpy <- merge(year, chp.ac, by = "season_starting", all = TRUE)

# Create a time stamp for year
chpy$ts <- as.POSIXct(strptime(paste(chpy$season_starting, "01-01", sep = "-"), format = "%Y-%m-%d", tz = "GMT"))

# List of unique site IDs
site_ids <- unique(chpy$site_id)

# Initialize an empty list to store lagged data frames for each site
lagged_dfs <- list()

# Loop through each site ID
for (site in site_ids) {
  # Filter data for the current site
  site_data <- chpy %>% filter(site_id == site)
  
  # Create xts objects for penguin count
  chpts1 <- xts(site_data$penguin_count, order.by = site_data$ts)
  
  # Determine valid lag values based on the number of observations
  max_lag <- nrow(site_data) - 1
  lags <- c(0, 10, 20, 30,40)
  valid_lags <- lags[lags <= max_lag]
  
  # Create lagged data frames
  lagged_df <- na.omit(data.frame(year = site_data$season_starting, site_id = site, 
                                  penguin_count_lag = lag.xts(chpts1, k = valid_lags)))
  
  # Append to the list
  lagged_dfs[[as.character(site)]] <- lagged_df
}

# Combine all lagged data frames 
chp_lagged_df <- bind_rows(lagged_dfs)


print(chp_lagged_df)


# gentoo penguin


year <- data.frame(season_starting = 1980:2022)

# Merge to include all years in the range
gepy <- merge(year, gep.ac, by = "season_starting", all = TRUE)

# time stamp for year
gepy$ts <- as.POSIXct(strptime(paste(gepy$season_starting, "01-01", sep = "-"), format = "%Y-%m-%d", tz = "GMT"))

# unique site IDs
site_ids <- unique(gepy$site_id)

# Initialize an empty list to store lagged data frames for each site
lagged_dfs <- list()

# Loop through each site ID
for (site in site_ids) {
  # Filter data for the current site
  site_data <- gepy %>% filter(site_id == site)
  
  # Create xts objects for penguin count
  gepts1 <- xts(site_data$penguin_count, order.by = site_data$ts)
  
  # Determine valid lag values based on the number of observations
  max_lag <- nrow(site_data) - 1
  lags <- c(0, 10, 20, 30,40)
  valid_lags <- lags[lags <= max_lag]
  
  # Create lagged data frames
  lagged_df <- na.omit(data.frame(year = site_data$season_starting, site_id = site, 
                                  penguin_count_lag = lag.xts(gepts1, k = valid_lags)))
  
  # Append to the list
  lagged_dfs[[as.character(site)]] <- lagged_df
}

# Combine all lagged data frames into one
gep_lagged_df <- bind_rows(lagged_dfs)

# View the final combined lagged data frame
print(gep_lagged_df)



###------ adelie penguin-----------


year <- data.frame(season_starting = 1980:2022)

# Merge to include all years in the range
adpy <- merge(year, adp.ac, by = "season_starting", all = TRUE)

# time stamp for year
adpy$ts <- as.POSIXct(strptime(paste(adpy$season_starting, "01-01", sep = "-"), format = "%Y-%m-%d", tz = "GMT"))

# unique site IDs
site_ids <- unique(adpy$site_id)


lagged_dfs <- list()

# Loop through each site ID
for (site in site_ids) {
  # Filter data for the current site
  site_data <- adpy %>% filter(site_id == site)
  
  # Create xts objects for penguin count
  adpts1 <- xts(site_data$penguin_count, order.by = site_data$ts)
  
  # Determine valid lag values based on the number of observations
  max_lag <- nrow(site_data) - 1
  lags <- c(0, 10, 20, 30,40)
  valid_lags <- lags[lags <= max_lag]
  
  # Create lagged data frames
  lagged_df <- na.omit(data.frame(year = site_data$season_starting, site_id = site, 
                                  penguin_count_lag = lag.xts(adpts1, k = valid_lags)))
  
  # Append to the list
  lagged_dfs[[as.character(site)]] <- lagged_df
}

# Combine lagged data frames 
adp_lagged_df <- bind_rows(lagged_dfs)


print(adp_lagged_df)

###------- calculation of change --------------

##---chinstrap change-------------
head(chp_lagged_df)
summary(as.factor(chp_lagged_df$year))

chpld<-subset(chp_lagged_df,year>2012)


summary(chpld$penguin_count_lag.lag0)

head(chpld)

chpld$change_10_years<-((chpld$penguin_count_lag.lag0/chpld$penguin_count_lag.lag10)-1)
chpld$change_20_years<-((chpld$penguin_count_lag.lag0/chpld$penguin_count_lag.lag20)-1)
chpld$change_30_years<-((chpld$penguin_count_lag.lag0/chpld$penguin_count_lag.lag30)-1)


chp.ld<-data.frame(chpld[1:2],chpld[9:11])
chpldm<-na.omit(melt(chp.ld,id.vars=c("year","site_id")))

# remove outliers

quartiles <- quantile(na.omit(chpldm$value), probs=c(.25, .75), na.rm = FALSE)

IQR <- IQR(chpldm$value)

Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR 

chp.out<- subset(chpldm, chpldm$value> Lower & chpldm$value < Upper)


ggplot(chp.out,aes(variable,value))+
  geom_hline(yintercept = 0,linetype="dashed")+
  geom_boxplot()



##---gentoo change-------------


head(gep_lagged_df)
summary(as.factor(gep_lagged_df$year))

gepld<-subset(gep_lagged_df,year>2012)

head(gepld)

gepld$change_10_years<-((gepld$penguin_count_lag.lag0/gepld$penguin_count_lag.lag10)-1)
gepld$change_20_years<-((gepld$penguin_count_lag.lag0/gepld$penguin_count_lag.lag20)-1)
gepld$change_30_years<-((gepld$penguin_count_lag.lag0/gepld$penguin_count_lag.lag30)-1)


gep.ld<-data.frame(gepld[1:2],gepld[9:11])
gepldm<-na.omit(melt(gep.ld,id.vars=c("year","site_id")))



# remove outliers

quartiles <- quantile(na.omit(gepldm$value), probs=c(.25, .75), na.rm = FALSE)

IQR <- IQR(gepldm$value)

Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR 

gep.out<- subset(gepldm, gepldm$value> Lower & gepldm$value < Upper)


ggplot(gep.out,aes(variable,value))+
  geom_hline(yintercept = 0,linetype="dashed")+
  geom_boxplot()



###----adelie change--

head(adp_lagged_df)
summary(as.factor(adp_lagged_df$year))

adpld<-subset(adp_lagged_df,year>2009)

head(adpld)

adpld$change_10_years<-((adpld$penguin_count_lag.lag0/adpld$penguin_count_lag.lag10)-1)
adpld$change_20_years<-((adpld$penguin_count_lag.lag0/adpld$penguin_count_lag.lag20)-1)
adpld$change_30_years<-((adpld$penguin_count_lag.lag0/adpld$penguin_count_lag.lag30)-1)


adp.ld<-data.frame(adpld[1:2],adpld[9:11])
adpldm<-na.omit(melt(adp.ld,id.vars=c("year","site_id")))


# remove outliers

quartiles <- quantile(na.omit(adpldm$value), probs=c(.25, .75), na.rm = FALSE)

IQR <- IQR(adpldm$value)

Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR 

adp.out<- subset(adpldm, adpldm$value> Lower & adpldm$value < Upper)

ggplot(adp.out,aes(variable,value))+
  geom_hline(yintercept = 0,linetype="dashed")+
  geom_boxplot()



### three species 
chp.out$species<-c("CHP")
gep.out$species<-c("GEP")
adp.out$species<-c("ADP")


change<-rbind(chp.out,gep.out,adp.out)

head(change)

ggplot(change,aes(species,value,fill=variable,linetype=variable))+
  geom_hline(yintercept=0,linetype="dashed")+
  geom_boxplot()+
  theme_bw()+
  ylab("proportional change in colony size")+
  scale_fill_manual(values=c("#0066CC","#339900","#FF3300"))
write.csv(change,"penguin_change.csv")

shapiro.test(chp.out$value)

bartlett.test(chp.out$value,chp.out$variable)

ggplot(chp.out,aes(scale(value)))+geom_histogram()

### ---------- Environmental Variables -------------
library(raster)
library(terra)
library(sf)
library(dplyr)
library(ncdf4)
library(tidyverse)

library(RNetCDF)

###---------CHL-----------------
ncpath <- "C:/Pygoscelis trends/"
ncname <- "cmems_mod_glo_bgc_my_0.25_P1M-m_1717610493617"
ncfname <- paste(ncpath, ncname, ".nc", sep="")
dname <- "chl"
chl_raster <- brick(ncfname, varname=dname)

# Extract time dimension values (assuming time dimension is the first dimension)
time_values <- getZ(chl_raster)


# Convert the raster brick to a dataframe
chldf <- as.data.frame(chl_raster, xy=TRUE, na.rm=TRUE)

# Display the first few rows of the dataframe
summary(chldf$x)

chldfm<-melt(chldf,id.vars=c("x","y"))

tail(chldfm)

# time stamp for year
chldfm$ts <- as.POSIXct(strptime(substring(chldfm$variable,first=2,last=11), format = "%Y.%m.%d", tz = "GMT"))

chldfm$quarter<-quarter(chldfm$ts)
chldfm$year<-year(chldfm$ts)

chldfm$Xint<-as.integer(chldfm$x)
chldfm$Yint<-as.integer(chldfm$y)



chlM<-ddply(subset(chldfm,quarter=="1"|quarter=="4"), c("year","Xint","Yint"), summarise,
             chlm=mean(value),
            chlsd=sd(value))  


ggplot(chlM,aes(Yint,chlm))+geom_point(alpha=0.15)+facet_wrap(year~.)+
  theme_bw()+xlab("Latitude")+ylab("Chl-a concentration (mg/m-3)")+
  coord_flip()


chl_high<-subset(chlM,chlm>1)

chlH<-ddply(chl_high, c("year"), summarise,
        chl.lat=mean(Yint),latse=sd(Yint)/sqrt(length(Yint)-1))  


ggplot(chlH,aes(year,chl.lat))+
  geom_errorbar(aes(ymin=chl.lat-latse,ymax=chl.lat+latse))+
  geom_smooth(method="lm")+
  geom_point() +
  ylab("mean latitude of chlorophyll-a concentration > 1mg/m-3")+
  theme_bw()


### ----------sea icwe-------------




ncpath <- "C:/Pygoscelis trends/"
ncname <- "cmems_mod_glo_phy_my_0.083deg-climatology_P1M-m_1717610549793"
ncfname <- paste(ncpath, ncname, ".nc", sep="")
dname <- "siconc"
sic_raster <- brick(ncfname, varname=dname)

# Extract time dimension values (assuming time dimension is the first dimension)
time_values <- getZ(chl_raster)


# Convert the raster brick to a dataframe
sicdf <- as.data.frame(chl_raster, xy=TRUE, na.rm=TRUE)

# Display the first few rows of the dataframe
summary(sicdf$x)

sicdfm<-melt(sicdf,id.vars=c("x","y"))

tail(sicdfm)


# time stamp for year
sicdfm$ts <- as.POSIXct(strptime(substring(sicdfm$variable,first=2,last=11), format = "%Y.%m.%d", tz = "GMT"))

sicdfm$quarter<-quarter(sicdfm$ts)
sicdfm$year<-year(sicdfm$ts)

sicdfm$Xint<-as.integer(sicdfm$x)
sicdfm$Yint<-as.integer(sicdfm$y)



sicM<-ddply(subset(sicdfm,quarter=="2"|quarter=="3"), c("year","Xint","Yint"), summarise,
            sicm=mean(value))  

ggplot(sicM,aes(Yint,sicm))+geom_point(alpha=0.15)+facet_wrap(year~.)+coord_flip()+


  theme_bw()+xlab("Latitude")+ylab("Winter sea ice cover")
  

sic_high<-subset(sicM,sicm>0.35)

sicH<-ddply(sic_high, c("year"), summarise,
            sic.lat=mean(Yint),latse=sd(Yint)/sqrt(length(Yint)-1))  


ggplot(sicH,aes(year,sic.lat))+
  geom_errorbar(aes(ymin=sic.lat-latse,ymax=sic.lat+latse))+
  geom_smooth(method="lm")+
  geom_point() +xlim(1990,2020)+
  ylab("mean latitude of winter sea ice cover >35%")+
  theme_bw()



##-------- krillbase-------------

kbase<-read.csv("C:/Pygoscelis trends/krillbase_data.csv")

head(kbase)


kbsp<-SpatialPointsDataFrame(kbase[6:5],kbase[1:28],proj4string = EPSG4326)

kb_6932<-sp::spTransform(kbsp,CRS("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))


kb.strata<-raster::intersect(kb_6932,strata)

kbdf<-data.frame(kb.strata)

head(kbdf)
kbdf$Yint<-as.integer(kbdf$LATITUDE)
kbdf<-subset(kbdf,SEASON>1989)



ggplot(kbdf,aes(Yint,STANDARDISED_KRILL_UNDER_1M2))+geom_point(alpha=0.15)+
  
  theme_bw()+xlab("Latitude")+ylab("Krill density")+
  facet_wrap(SEASON~.)+coord_flip()

summary(kbdf$STANDARDISED_KRILL_UNDER_1M2)

head(kbdf)

kbdf$Xint<-as.integer(kbdf$LONGITUDE)


meandens<-ddply(subset(kbdf,STANDARDISED_KRILL_UNDER_1M2>0.16), c("SEASON","ID","Yint","Xint"), summarise,
                stdk=mean(STANDARDISED_KRILL_UNDER_1M2))  
mean(meandens$stdk)
krillchange<-ddply((kbdf), c("ID","Yint","Xint"), summarise,
                kchange=mean(STANDARDISED_KRILL_UNDER_1M2-34.72))  

write.csv(krillchange,"krill_change.csv")

ggplot(meandens,aes(SEASON,stdk))+geom_point()+facet_wrap(ID~.) 
# BS, EI,GS, JOIN on a smooth decrease, what is compensated by MB y PB


krilllat<-ddply(subset(kbdf,STANDARDISED_KRILL_UNDER_1M2>1.37), c("SEASON"), summarise,
                lat=mean(LATITUDE),latse=sd(LATITUDE)/sqrt(length(LATITUDE)-1))  


ggplot(krilllat,aes(SEASON,lat))+
  geom_smooth(method="gam")+
  geom_errorbar(aes(ymin=lat-latse,ymax=lat+latse))+
  geom_point()

## ----fishing (CCAMLR statistical bulletin)--------

ccamlr<-read.csv("C:/Pygoscelis trends/ccamlr_fishing_catch.csv")

catch<-subset(ccamlr,effort_or_catch=="catch")
summary(as.factor(catch$taxon_code))

kcatch<-subset(catch,taxon_code=="KRI")

summary(as.factor(kcatch$asd_code))

k481<-subset(kcatch,asd_code=="481")


head(k481)

ggplot(k481,aes(year,greenweight_caught_tonne))+
  xlim(1990,2020)+
  geom_point()


kcm<-ddply(k481, c("year"), summarise,
           Tcatch=sum(greenweight_caught_tonne))  

ggplot(kcm,aes(year,Tcatch))+
  xlim(1990,2020)+
  geom_point()


###-----all plots ----------------



krillat2<-subset(krilllat,SEASON<2012)


ggplot()+
  geom_vline(xintercept=2010,linetype="dotdash")+
  geom_smooth(data=krilllat,aes(x=SEASON,y=lat),method="gam",
              fullrange=T,colour="blue",linetype="solid",se=F)+
  geom_smooth(data=krillat2,aes(x=SEASON,y=lat),method="gam",
              fullrange=T,colour="red",linetype="dashed",se=F)+
  geom_errorbar(data=krilllat,aes(x=SEASON,y=lat,ymin=lat-latse,ymax=lat+latse))+
  geom_point(data=krilllat,aes(x=SEASON,y=lat),)+theme_bw()+xlim(1990,2020)+
  ylab("")+
  xlab("year")+ggtitle(label="a.mean latitude when krill density > 27-y median")+
  
  
  
  
  
  ggplot(chlH,aes(year,chl.lat))+
  geom_vline(xintercept=2010,linetype="dotdash")+
  geom_errorbar(aes(ymin=chl.lat-latse,ymax=chl.lat+latse))+
  geom_smooth(method="lm",se=F)+
  geom_point() +
  ylab("")+
  theme_bw()+ggtitle(label="b.mean latitude of chlorophyll-a concentration > 1mg/m-3")+
  
  
  
  ggplot(sicH,aes(year,sic.lat))+
  geom_vline(xintercept=2010,linetype="dotdash")+
  geom_errorbar(aes(ymin=sic.lat-latse,ymax=sic.lat+latse))+
  geom_smooth(method="lm",se=F)+
  geom_point() +xlim(1990,2020)+
  ylab("")+
  theme_bw()+ggtitle(label="c.mean latitude of winter sea ice cover >35%")+
  
  ggplot(kcm,aes(year,Tcatch))+
  geom_smooth(se=F)+
  theme_bw()+
  ylab("ton")+
  xlim(1990,2020)+
  geom_point()+
  ggtitle(label="d. accumulated fishing catch")



### -------- Path Analysis ---------------
library(lavaan)
library(semPlot)
head(chlH)

head(sicH)
head(krilllat)
head(kcm)

kl<-data.frame(year=krilllat$SEASON,krill.lat=krilllat$lat)

df1<-merge(sicH[1:2],chlH[1:2],by="year")
df2<-merge(df1,kcm,by="year")
df3<-merge(df2,kl,by="year")


head(df3)





# time stamp for year
df3$ts <- as.POSIXct(strptime(paste(df3$year, "01-01", sep = "-"), format = "%Y-%m-%d", tz = "GMT"))

# unique site IDs
years <- unique(df3$year)


lagged_dfs <- list()

# Loop through each site ID
for (year in years) {
  # Filter data for the current site
  year_data <- df3 %>% filter(years == year)
  
  # Create xts objects for penguin count
  df3.1 <- xts(year_data$Tcatch, order.by =year_data$ts)
  
  # Determine valid lag values based on the number of observations
  max_lag <- nrow(year_data) - 1
  lags <- c(0, 1)
  valid_lags <- lags[lags <= max_lag]
  
  # Create lagged data frames
  lagged_df <- na.omit(data.frame(year = year_data$year, 
                                  krill=year_data$krill.lat,
                                  chl=year_data$chl.lat, 
                                  sicL0=year_data$sic.lat,
                                  catch_lag = lag.xts(df3.1, k = valid_lags)))
  
  # Append to the list
  lagged_dfs[[as.character(site)]] <- lagged_df
}
kf_lagged_df <- bind_rows(lagged_dfs)


head(kf_lagged_df)

# unique site IDs
years <- unique(df3$year)


lagged_dfs <- list()

# Loop through each site ID
for (year in years) {
  # Filter data for the current site
  year_data <- df3 %>% filter(years == year)
  
  # Create xts objects for penguin count
  df3.1 <- xts(year_data$sic.lat, order.by =year_data$ts)
  
  # Determine valid lag values based on the number of observations
  max_lag <- nrow(year_data) - 1
  lags <- c(0, 1)
  valid_lags <- lags[lags <= max_lag]
  
  # Create lagged data frames
  lagged_df <- na.omit(data.frame(year = year_data$year, 
                                  sic_lag = lag.xts(df3.1, k = valid_lags)))
  
  # Append to the list
  lagged_dfs[[as.character(site)]] <- lagged_df
}
s_lagged_df <- bind_rows(lagged_dfs)


head(s_lagged_df)

mdf<-merge(s_lagged_df,kf_lagged_df)

head(mdf)

ggplot(mdf,aes(catch_lag.lag1,krill))+geom_point()
ggplot(mdf,aes(sic_lag.lag1,krill))+geom_point()
ggplot(mdf,aes(chl,krill))+geom_point()


sdf<-data.frame(krill=scale(mdf$krill),chl=scale(mdf$chl),
                sicL0=scale(mdf$sic_lag.lag0),sicL1=scale(mdf$sic_lag.lag1),
                catchL0=scale(mdf$catch_lag.lag0),catchL1=scale(mdf$catch_lag.lag1))


head(sdf)





## path analysis example 
model1 <- '

# regressions
krill ~ catchL1+sicL1+chl+catchL0+sicL0

catchL1~sicL1

catchL0~chl+sicL1+catchL1+sicL0

chl~sicL1+sicL0

# residual correlations


sicL1 ~~ sicL0

'

sem1<-sem(model = model1, data = sdf)
varTable(sem1)

summary(sem1, fit.measures = TRUE)

splot1<-semPlot::semPlotModel(sem1)



cfa1<-cfa(model = model1, data = sdf)

summary(cfa1)

semPaths(splot1,
         intercepts = TRUE, residuals = TRUE, thresholds = TRUE, intStyle = "multi", 
         rotation = 1,  curvature = 1, nCharNodes = 3, nCharEdges = 3, structural = FALSE, 
         ThreshAtSide = FALSE,  
         thresholdSize = 0.5, fixedStyle = 2, freeStyle = 1, 
         as.expression = character(0), optimizeLatRes = FALSE, inheritColor = TRUE, 
         pastel = FALSE, rainbowStart = 0, 
         springLevels = FALSE, nDigits = 2, exoCov = TRUE, centerLevels = TRUE, 
         panelGroups = FALSE, layoutSplit = FALSE, measurementLayout = "tree",
         subRes = 4, modelOpts = list(mplusStd = "std"), 
         curveAdjacent = '<->', edge.label.cex = 0.6,  cardinal = "none", 
         equalizeManifests = FALSE, covAtResiduals = TRUE,  optimPoints = 1:8 * (pi/4))




