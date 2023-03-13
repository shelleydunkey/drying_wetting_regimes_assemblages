#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: nmds w/full environ (final)
#Coder: Michelle Busch
#Date: 20221006
#Purpose: how do community compositions change with different events? What 
# environmental variables are driving those compositional changes?
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 1: Setup workspace -------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Clear memory
remove(list=ls())

# setting working directory

# load libraries
library(tidyverse)
library(vegan)
library(lubridate)
library(factoextra)
library(wesanderson)
library(RColorBrewer)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 2: Format Environ Data ---------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1. FULL SET OF SAMPLING EVENTS (30) 
hydro_metrics_30 <- read.csv("FINAL/final_metrics_full_events.csv", 
                             header = TRUE) %>% 
  rename(stationcode = site) %>% 
  dplyr::select(-c(hydro_metrics, notes))

# load hab / chem data 
# not lab chem data! can't have missing data, non consistent across all sites
hab <- read.csv("FINAL/final_habitat_matrix.csv", header = TRUE) %>% 
  dplyr::select(-c(X))

chem <- read.csv("FINAL/final_field_chem_matrix.csv", header = TRUE) %>% 
  dplyr::select(-c(X))

# need to convert date metrics into Julian Day - Hydro Metrics 30
  # yday()
## UPDATE - converting dates to WATER YEAR

# first to date format
hydro_metrics_30 <- hydro_metrics_30 %>% 
  separate(SampleDate, into = c("month", "day", "year"), sep = "/", 
           remove = FALSE) %>% 
  select(-c(month, day))


hydro_metrics_30$dry_date <- (as.Date(hydro_metrics_30$dry_date,
                                                 "%m/%d/%Y"))

hydro_metrics_30$peak_date <- (as.Date(hydro_metrics_30$peak_date,
                                             "%m/%d/%Y"))

hydro_metrics_30$rewet_date <- (as.Date(hydro_metrics_30$rewet_date,
                                              "%m/%d/%Y"))

hydro_metrics_30$SampleDate <- (as.Date(hydro_metrics_30$SampleDate,
                                              "%m/%d/%Y"))

hydro_metrics_30$first_wet_date <- (as.Date(hydro_metrics_30$first_wet_date,
                                            "%m/%d/%Y"))

# convert to day of water year
hydro.day.new = function(x, start.month = 10L){
  start.yr = year(x) - (month(x) < start.month)
  start.date = make_date(start.yr, start.month, 1L)
  as.integer(x - start.date + 1L)
  }

hydro_metrics_30$dry_date <- hydro.day.new(hydro_metrics_30$dry_date)

hydro_metrics_30$peak_date <- hydro.day.new(hydro_metrics_30$peak_date)

hydro_metrics_30$rewet_date <- hydro.day.new(hydro_metrics_30$rewet_date)

hydro_metrics_30$SampleDate <- hydro.day.new(hydro_metrics_30$SampleDate)

hydro_metrics_30$first_wet_date <- hydro.day.new(hydro_metrics_30$first_wet_date)

summary(hydro_metrics_30)
## gives you different values! so I guess it worked :D


# unite the stationcode and sample day for each sampling event - 1 row per event
hydro_metrics_30 <- hydro_metrics_30 %>% 
  unite(col = sample_event, c(stationcode, SampleDate), 
        sep = "_", remove = FALSE) 

## format Chem and Hab data to unite
# field chemistry
chem <- chem %>% 
  separate(event_id, into = c("stationcode", "SampleDate"), sep = "_")

chem$SampleDate <- hydro.day.new(as.Date(chem$SampleDate,
                                  "%Y-%m-%d"))
chem <- chem %>% 
  unite(col = sample_event, c(stationcode, SampleDate), 
        sep = "_", remove = FALSE) 

# habitat
hab <- hab %>% 
  separate(event_id, into = c("stationcode", "SampleDate"), sep = "_")

hab$SampleDate <- hydro.day.new(as.Date(hab$SampleDate,
                                 "%Y-%m-%d"))
hab <- hab %>% 
  unite(col = sample_event, c(stationcode, SampleDate), 
        sep = "_", remove = FALSE) 

## join it all together, based on hydro data
## this will make sure everything is lined up for future anaylses
env_variables_30 <- left_join(hydro_metrics_30, chem, 
                              by = c("stationcode", "SampleDate",
                                                             "sample_event"))

env_variables_30 <- left_join(env_variables_30, hab, 
                              by = c("stationcode", "SampleDate",
                                                            "sample_event"))

# need to remove any data with NA's
# expect to lose: 0 hydro metrics, 1 hab metrics (substrate.size.class), 
# 3 field chem (DO, Salinity, Turbidity) -> 4 total
env_variables_30 <- env_variables_30 %>% 
  select(-c(Oxygen..Dissolved, Salinity, Turbidity, Substrate.Size.Class)) 
# 30 events, 4 col lost :)


rm(chem)
rm(hab)
rm(hydro_metrics_30)


#### 2.c Correlations ####
library("Hmisc")
library("PerformanceAnalytics")


# full events:
env_variables_30_corr <- env_variables_30 %>% 
  dplyr:: select(-c(stationcode, SampleDate, sample_event,
                    sampling_date, year, antecedent_time_period_days, 
                    false_starts))

# look at correlations of hydro metrics
hydro_cor <- env_variables_30_corr[, 1:17]

corr <- chart.Correlation(hydro_cor, histogram=TRUE, pch=4)

# we know we want to keep total_rewet (H1), remvoe those correlated variables
hydro_cor <- hydro_cor %>% 
  select(-c(peak2sample_duration, prop_wet_days, prop_dry_days, rewet_date))

# try again
corr <- chart.Correlation(hydro_cor, histogram=TRUE, pch=4)
# looks good! all under |0.70|

# remove those from full environ
env_variables_30_corr <- env_variables_30_corr %>% 
  select(-c(peak2sample_duration, prop_wet_days, prop_dry_days, rewet_date))

corr <- chart.Correlation(env_variables_30_corr, histogram=TRUE, pch=4)

# remove station water depth
env_variables_30_corr <- env_variables_30_corr %>% 
  select(-StationWaterDepth)

# final check
corr <- chart.Correlation(env_variables_30_corr, histogram=TRUE, pch=4)

# MULTICOLLINARITY MATTERS! #
# for 30 sites: 
env_variables_30 <- env_variables_30 %>% 
  dplyr:: select(-c(peak2sample_duration, prop_wet_days, prop_dry_days, 
                    rewet_date, StationWaterDepth))

rm(env_variables_30_corr)
rm(hydro_cor)
rm(corr)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 3: Load / Format Biol Data -----------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### 3.a BMI (rel) ####
# format dates to combine
# rel abund bmi
bmi_sample <- read.csv("FINAL/final_bmispmatrix_rep1_duplicates_relabund.csv",
                       header = TRUE) %>% 
  dplyr::select(-X)

bmi_sample$SampleDate <- hydro.day.new(as.Date(bmi_sample$SampleDate,
                                        "%Y-%m-%d"))

bmi_sample <- bmi_sample %>% 
  unite(col = sample_event, c(stationcode, SampleDate), 
        sep = "_", remove = FALSE)  

# combine BMI with environ
full_bugs_30 <- left_join(env_variables_30, bmi_sample, by = c("stationcode", 
                                                               "SampleDate",
                                                               "sample_event"))

## write these!
# write.csv(full_bugs_30, "FINAL/final_bmi_env_30.csv")


#### 3.b Algae ####
# only quantitative soft bodied algae!
algae_sample <- read.csv("FINAL/final_SoftAlgaeMatrix_rep1_duplicates_relabund.csv",
                         header = TRUE) %>% 
  dplyr::select(-X) 

# format dates to combine
# algae
algae_sample$SampleDate <- hydro.day.new(as.Date(algae_sample$SampleDate,
                                          "%m/%d/%Y"))

algae_sample <- algae_sample %>% 
  unite(col = sample_event, c(stationcode, SampleDate), 
        sep = "_", remove = FALSE)  

# combine BMI with enviorn
full_algae_30 <- left_join(env_variables_30, algae_sample, by = c("stationcode", 
                                                                  "SampleDate",
                                                                  "sample_event"))


#### 3.b.1 qual Algae ####

algae_qual <- read.csv("FINAL/final_SoftAlgaeMatrix_rep1_dupliucates_full_qual.csv",
                       header = TRUE) %>% 
  dplyr::select(-X)

# format dates to combine
# algae
algae_qual <- algae_qual %>% 
  separate(SampleDate, into = c("month", "day", "year"), sep = "/", 
           remove = FALSE) %>% 
  select(-c(month, day))

algae_qual$SampleDate <- hydro.day.new(as.Date(algae_qual$SampleDate,
                                        "%m/%d/%Y"))

algae_qual <- algae_qual %>% 
  unite(col = sample_event, c(stationcode, SampleDate), 
        sep = "_", remove = FALSE)  

# combine algae with enviorn
algae_qual_30 <- left_join(env_variables_30, algae_qual, by = c("stationcode", 
                                                                  "SampleDate",
                                                                  "sample_event",
                                                                "year"))



#### 3.c Diatoms ####
diatom_sample <- read.csv("FINAL/final_Diatomspmatrix_rep1_duplicates_relabund.csv",
                          header = TRUE)%>% 
  dplyr::select(-X) 

# format dates to combine
# diatoms
diatom_sample$SampleDate <- hydro.day.new(as.Date(diatom_sample$SampleDate,
                                           "%m/%d/%Y"))

diatom_sample <- diatom_sample %>% 
  unite(col = sample_event, c(stationcode, SampleDate), 
        sep = "_", remove = FALSE)  

# combine BMI with enviorn
full_diatom_30 <- left_join(env_variables_30, diatom_sample, by = c("stationcode", 
                                                                    "SampleDate",
                                                                    "sample_event"))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 4: NMDS BMI (rel) --------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Clear memory
remove(list=ls())

# load
bmi_30 <- read.csv("FINAL/final_bmi_env_30.csv", header = TRUE) %>%
  select(-c(X, sampling_date, false_starts))

# remove May 2014 BELOLV sample
bmi_30 <- bmi_30 %>%
  filter(year != 2014) %>%
  filter(year != 2018)

bmi_30_num <- bmi_30[ , c(28:263)]

bmi_30_num <- bmi_30_num[,colSums(bmi_30_num/bmi_30_num, na.rm=T)>0] #152

# Running NMDS in vegan (metaMDS)
bio.spp_NMS <-
  metaMDS(log(bmi_30_num + 1),
          distance = "bray",
          k = 3,
          maxit = 999,
          trymax = 500,
          autotransform = FALSE)
# stress = 0.1272828  

stressplot(bio.spp_NMS)
# non metric r2 0.984
# linear r2 = 0.887


hydro <- bmi_30[, c(6:27)]
hydro <- hydro %>% 
  select(-c(recession_coef, recession_rsq, rewet_duration))

## trying to rotate
bio.spp_NMS <- with(hydro, MDSrotate(bio.spp_NMS, total_rewet))

data.scores <- as.data.frame(scores(bio.spp_NMS$points))
species.scores <- as.data.frame(scores(bio.spp_NMS$species))

#add columns to data frame
data.scores$site = bmi_30$stationcode
data.scores$date = bmi_30$SampleDate
data.scores$year = bmi_30$year

head(data.scores)


hydro.envfit <- envfit(bio.spp_NMS, hydro, permutations = 999, na.rm = TRUE)

# running environ fit #extracts relevant scores from envifit
env.scores <- as.data.frame(scores(hydro.envfit, display = "vectors"))
#and then gives them their names
env.scores <- cbind(env.scores, env.variables = rownames(env.scores))

# add pvalues to dataframe
env.scores <- cbind(env.scores, pval = hydro.envfit$vectors$pvals)
# add r2
env.scores <- cbind(env.scores, r2 = hydro.envfit$vectors$r)
# comment from co-authors - display r value rather than r2, testing
env.scores <- cbind(env.scores, r = hydro.envfit$vectors$arrows)

hydro.envfit$vectors$control


#### Plot ####

# make colored as a factor so it doesn't read as continuous
data.scores$year <- as.factor(data.scores$year)


sig.env.scrs <- subset(env.scores, pval<0.015)

# year color
yearColor <- wes_palette("FantasticFox1", 3, "discrete")
names(yearColor) <- levels(data.scores$year)
custom_colors <- scale_colour_manual(name = "Year", values = yearColor)

# plot by year 
ggplot(data.scores, aes(x=MDS1, y=MDS2)) + #sets up the plot
  geom_point(aes(MDS1, MDS2, colour = factor(data.scores$year)), size = 8) + #adds site points to plot, shape determined by Landuse, colour determined by Management
  custom_colors +
  theme_classic()+ 
  theme(axis.text.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 20),
        legend.text = element_text(size = 12, face ="bold", colour ="black"),
        legend.position = "right", axis.title.y = element_text(face = "bold", 
                                                               size = 14),
        axis.title.x = element_text(face = "bold", size = 18, colour = "black"),
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) +
  xlim(-1, 1) +
  ylim(-1, 1) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), #add environ data!
               data = sig.env.scrs, size = 1.5, alpha = 0.5, colour = "grey19") # +
  # geom_text(data = sig.env.scrs, aes(x = NMDS1, y = NMDS2), colour = "black",
  #           fontface = "bold", size = 6, label = sig.env.scrs$env.variables)


### MRPP - influence of year
year = bmi_30[ , 4]

mrpp(log(bmi_30_num + 1), group = year, distance = "bray")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 5: NMDS Algae ------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# quant samples
alg_30 <- read.csv("FINAL/final_algae_env_30.csv", header = TRUE) %>%
   select(-c(X, sampling_date, false_starts))
#
# # remove May 2014 BELOLV sample
alg_30 <- alg_30 %>%
  filter(year != 2014) %>% 
  filter(year != 2018) # %>% 
  # mutate(total_rewet = peak2sample_duration + rewet_duration)

alg_30_num <- alg_30[ , c(28:396)]

alg_30_num <- alg_30_num[,colSums(alg_30_num/alg_30_num, na.rm=T)>0] #207

# Running NMDS in vegan (metaMDS)
bio.spp_NMS <-
  metaMDS(log(alg_30_num + 1),
          distance = "bray",
          k = 3,
          maxit = 999,
          trymax = 500,
          autotransform = FALSE)
# # stress = 0.1803479   
#
stressplot(bio.spp_NMS)
# # non metric r2 0.967
# # linear r2 = 0.701


# # environ data - need to edit!
hydro <- alg_30[, c(6:27)]
hydro <- hydro %>% 
  select(-c(recession_coef, recession_rsq, rewet_duration))


## trying to rotate
bio.spp_NMS <- with(hydro, MDSrotate(bio.spp_NMS, total_rewet))

data.scores = as.data.frame(scores(bio.spp_NMS$points))
species.scores = as.data.frame(scores(bio.spp_NMS$species))

#add columns to data frame
data.scores$site = alg_30$stationcode
data.scores$date = alg_30$SampleDate
data.scores$year = alg_30$year

head(data.scores)

hydro.envfit <- envfit(bio.spp_NMS, hydro, permutations = 999, na.rm = TRUE)

# running environ fit #extracts relevant scores from envifit
env.scores <- as.data.frame(scores(hydro.envfit, display = "vectors"))
#and then gives them their names
env.scores <- cbind(env.scores, env.variables = rownames(env.scores))

# add pvalues to dataframe
env.scores <- cbind(env.scores, pval = hydro.envfit$vectors$pvals)
# add r2 to dataframe
env.scores <- cbind(env.scores, r2 = hydro.envfit$vectors$r)


#### Plot ####

# join with data.scores
data.scores <- left_join(data.scores, spatial)

# make colored as a factor so it doesn't read as continuous
data.scores$year <- as.factor(data.scores$year)



env.scores <- read.csv("FINAL/extra_info/env.scores.alg30_NO14_NO18.csv", header = TRUE) %>% 
  select(-X)

sig.env.scrs <- subset(env.scores, pval<0.015)


# year color
yearColor <- wes_palette("FantasticFox1", 3, "discrete")
names(yearColor) <- levels(data.scores$year)
custom_colors <- scale_colour_manual(name = "Year", values = yearColor)

# plot by year 
ggplot(data.scores, aes(x=MDS1, y=MDS2)) + #sets up the plot
  geom_point(aes(MDS1, MDS2, colour = factor(data.scores$year)), size = 8) + #adds site points to plot, shape determined by Landuse, colour determined by Management
  custom_colors +
  theme_classic()+ 
  theme(axis.text.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 20),
        legend.text = element_text(size = 12, face ="bold", colour ="black"),
        legend.position = "right", axis.title.y = element_text(face = "bold", 
                                                               size = 14),
        axis.title.x = element_text(face = "bold", size = 18, colour = "black"),
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) +
  xlim(-1, 1) +
  ylim(-1, 1) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), #add environ data!
               data = sig.env.scrs, size = 1, alpha = 0.5, colour = "grey30")# + 
# geom_text(data = sig.env.scrs, aes(x = NMDS1, y = NMDS2), colour = "black",
#           fontface = "bold", size = 10, label = sig.env.scrs$env.variables)

### MRPP - influence of year
year = alg_30[ , 4]

mrpp(log(alg_30_num + 1), group = year, distance = "bray")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 6: NMDS Diatom -----------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dia_30 <- read.csv("FINAL/final_diatom_env_30.csv", header = TRUE) %>%
  select(-c(X, sampling_date, false_starts))

# remove May 2014 BELOLV sample
dia_30 <- dia_30 %>%
   filter(year != 2014) %>% 
   filter(year != 2018)

dia_30_num <- dia_30[ , c(28:403)]

dia_30_num <- dia_30_num[,colSums(dia_30_num/dia_30_num, na.rm=T)>0] #225


# Running NMDS in vegan (metaMDS)
bio.spp_NMS <-
  metaMDS(log(dia_30_num + 1),
          distance = "bray",
          k = 3,
          maxit = 999,
          trymax = 500,
          autotransform = FALSE)
# stress = 0.1483523  


stressplot(bio.spp_NMS)
# non metric r2 0.978
# linear r2 = 0.818


## get hydro ready
# # environ data
hydro <- dia_30[, c(6:27)]
# # exclude sample event, stationcode, sampledate, sampling_date
hydro <- hydro %>% 
  select(-c(recession_coef, recession_rsq, rewet_duration))

## trying to rotate
bio.spp_NMS <- with(hydro, MDSrotate(bio.spp_NMS, total_rewet))


# make the files
data.scores = as.data.frame(scores(bio.spp_NMS$points))
species.scores = as.data.frame(scores(bio.spp_NMS$species))

#add columns to data frame
data.scores$site = dia_30$stationcode
data.scores$date = dia_30$SampleDate
data.scores$year = dia_30$year

head(data.scores)

hydro.envfit <- envfit(bio.spp_NMS, hydro, permutations = 999, na.rm = TRUE)

# # running environ fit #extracts relevant scores from envifit
env.scores <- as.data.frame(scores(hydro.envfit, display = "vectors"))
#and then gives them their names
env.scores <- cbind(env.scores, env.variables = rownames(env.scores))

# add pvalues to dataframe
env.scores <- cbind(env.scores, pval = hydro.envfit$vectors$pvals)
# add r2
env.scores <- cbind(env.scores, r2 = hydro.envfit$vectors$r)


#### Plot ####

# join with data.scores
data.scores <- left_join(data.scores, spatial)

# make colored as a factor so it doesn't read as continuous
data.scores$year <- as.factor(data.scores$year)


sig.env.scrs <- subset(env.scores, pval<0.015)


# year color
yearColor <- wes_palette("FantasticFox1", 3, "discrete")
names(yearColor) <- levels(data.scores$year)
custom_colors <- scale_colour_manual(name = "Year", values = yearColor)

# plot by year 
ggplot(data.scores, aes(x=MDS1, y=MDS2)) + #sets up the plot
  geom_point(aes(MDS1, MDS2, colour = factor(data.scores$year)), size = 8) + #adds site points to plot, shape determined by Landuse, colour determined by Management
  custom_colors +
  theme_classic()+ 
  theme(axis.text.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 20),
        legend.text = element_text(size = 12, face ="bold", colour ="black"),
        legend.position = "right", axis.title.y = element_text(face = "bold", 
                                                               size = 14),
        axis.title.x = element_text(face = "bold", size = 18, colour = "black"),
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) +
  xlim(-1, 1) +
  ylim(-1, 1) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), #add environ data!
               data = sig.env.scrs, size = 1, alpha = 0.5, colour = "grey30") # +
# geom_text(data = sig.env.scrs, aes(x = NMDS1, y = NMDS2), colour = "black",
#           fontface = "bold", size = 10, label = sig.env.scrs$env.variables)

### MRPP - influence of year
year = dia_30[ , 4]

mrpp(log(dia_30_num + 1), group = year, distance = "bray")

