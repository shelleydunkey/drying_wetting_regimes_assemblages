#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: modeling_final
#Coder: Michelle Busch
#Date: 20221006
#Purpose: What hydrologic metrics are most predictive of assemblage diversity?
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## final thoughts - removing 2014, 2018 samples
## 27 events total ##

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 1: Setup workspace -------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Clear memory
remove(list=ls())

#setting working directory

#loading libraries
require(tidyverse)
require(vegan)
require(hillR)
library(lme4)
require(lmerTest)
library(performance)
require(psych)
library(lubridate)
library(MuMIn)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 2: Calculate metrics for models: BMI (rel abund) -------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bmi_full <- read.csv("FINAL/final_bmi_env_30.csv", header = TRUE) %>%
  select(-c(X, sampling_date, false_starts, recession_coef, recession_rsq))

# remove May 2014 BELOLV sample
bmi_full <- bmi_full %>%
  filter(year != 2014) %>%
  filter(year != 2018)

bmi_sp <- bmi_full[ , c(26:261)]

bmi_sp <- bmi_sp[,colSums(bmi_sp/bmi_sp, na.rm=T)>0] #152

# create new matrix for metrics
bmi_metrics <- bmi_full[ , c(1:3)]

# 1. species richness
bmi_metrics$specrich <- specnumber(bmi_sp, MARGIN = 1)
# MARGIN = 1 b/c sp in columns ( = 2 if in rows)

# 2. shannons
bmi_metrics$shannon <- diversity(bmi_sp, index = "shannon", MARGIN = 1, 
                                 base = exp(1))

# 3. Evenness. There isn't a function in vegan, but it is just shannon 
# divided by the log of richness so we can calculate it by hand:
bmi_metrics$evenness <- bmi_metrics$shannon / (log(bmi_metrics$specrich))

summary(bmi_metrics)
# evenness all less than 1!

# q = 1 ; shannon
bmi_metrics$hill_shannon <- hill_taxa(bmi_sp, q = 1, MARGIN = 1, base = exp(1))

# q = 2 ; simpson
bmi_metrics$hill_simp <- hill_taxa(bmi_sp, q = 2, MARGIN = 1, base = exp(1))

# recombine with hydro metrics to have everything in one place
hydro <- bmi_full[ , c(1:25)]

bmi_metrics <- left_join(hydro, bmi_metrics, by = c("stationcode",
                                                    "SampleDate",
                                                    "sample_event"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 3: Calculate metrics for models: Algae -----------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# quant samples
alg_full <- read.csv("FINAL/final_algae_env_30.csv", header = TRUE) %>%
  select(-c(X, sampling_date, false_starts))
#
# # remove May 2014 BELOLV sample
alg_full <- alg_full %>%
  filter(year != 2014) %>% 
  filter(year != 2018) # %>% 
# mutate(total_rewet = peak2sample_duration + rewet_duration)

alg_sp <- alg_full[ , c(28:396)]

alg_sp <- alg_sp[,colSums(alg_sp/alg_sp, na.rm=T)>0] #207

# create new matrix for metrics
alg_metrics <- alg_full[ , c(1:3)]

# 1. species richness
alg_metrics$specrich <- specnumber(alg_sp, MARGIN = 1)
# MARGIN = 1 b/c sp in columns ( = 2 if in rows)

# 2. shannons
alg_metrics$shannon <- diversity(alg_sp, index = "shannon", MARGIN = 1, 
                                 base = exp(1))

# 3. Evenness. There isn't a function in vegan, but it is just shannon 
# divided by the log of richness so we can calculate it by hand:
alg_metrics$evenness <- alg_metrics$shannon / (log(alg_metrics$specrich))

summary(alg_metrics)
# evenness all less than 1!

# q = 1 ; shannon
alg_metrics$hill_shannon <- hill_taxa(alg_sp, q = 1, MARGIN = 1, base = exp(1))

# q = 2 ; simpson
alg_metrics$hill_simp <- hill_taxa(alg_sp, q = 2, MARGIN = 1, base = exp(1))

# recombine with hydro metrics to have everything in one place
hydro <- alg_full[ , c(1:27)]

alg_metrics <- left_join(hydro, alg_metrics, by = c("stationcode",
                                                    "SampleDate",
                                                    "sample_event"))

## qualitative data ##
# load matrix
qual_matrix <- read.csv("FINAL/final_SoftAlgaeMatrix_rep1_dupliucates_full_qual.csv",
                        header = TRUE) %>% 
  select(-X)

qual_sp <- qual_matrix[ , c(3:385)]
qual_sp <- qual_sp[,colSums(qual_sp/qual_sp, na.rm=T)>0] # 383

# calc spec rich
alg_qual_metric <- qual_matrix[ , c(1:2)]

# 1. species richness
alg_qual_metric$specrich_qual <- specnumber(qual_sp, MARGIN = 1)

# 2. need to convert date to julian day to combine, this was made with the full
  # 85 observations
# but first need to pull our year
alg_qual_metric <- alg_qual_metric %>% 
  separate(SampleDate, into = c("month", "day", "year"), remove = FALSE) %>% 
  select(-c(month, day))

alg_qual_metric$year <- as.integer(alg_qual_metric$year)  

## day of water year
hydro.day.new = function(x, start.month = 10L){
  start.yr = year(x) - (month(x) < start.month)
  start.date = make_date(start.yr, start.month, 1L)
  as.integer(x - start.date + 1L)
}

alg_qual_metric$SampleDate <- hydro.day.new(as.Date(alg_qual_metric$SampleDate,
                                        "%m/%d/%Y"))

alg_metrics <- left_join(alg_metrics, alg_qual_metric)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 4: Calculate metrics for models: Diatoms ---------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dia_full <- read.csv("FINAL/final_diatom_env_30.csv", header = TRUE) %>%
  select(-c(X, sampling_date, false_starts))

# remove May 2014 BELOLV sample
dia_full <- dia_full %>%
  filter(year != 2014) %>% 
  filter(year != 2018)

dia_sp <- dia_full[ , c(28:403)]

dia_sp <- dia_sp[,colSums(dia_sp/dia_sp, na.rm=T)>0] #225

# create new matrix for metrics
dia_metrics <- dia_full[ , c(1:3)]

# 1. species richness
dia_metrics$specrich <- specnumber(dia_sp, MARGIN = 1)
# MARGIN = 1 b/c sp in columns ( = 2 if in rows)

# 2. shannons
dia_metrics$shannon <- diversity(dia_sp, index = "shannon", MARGIN = 1, 
                                 base = exp(1))

# 3. Evenness. There isn't a function in vegan, but it is just shannon 
# divided by the log of richness so we can calculate it by hand:
dia_metrics$evenness <- dia_metrics$shannon / (log(dia_metrics$specrich))

summary(dia_metrics)
# evenness all less than 1!

# q = 1 ; shannon
dia_metrics$hill_shannon <- hill_taxa(dia_sp, q = 1, MARGIN = 1, base = exp(1))

# q = 2 ; simpson
dia_metrics$hill_simp <- hill_taxa(dia_sp, q = 2, MARGIN = 1, base = exp(1))

# recombine with hydro metrics to have everything in one place
hydro <- dia_full[ , c(1:27)]

dia_metrics <- left_join(hydro, dia_metrics, by = c("stationcode",
                                                    "SampleDate",
                                                    "sample_event"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 5: Testing the effects of year -------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### 5.a relative abundance BMI ####
bmi <- read.csv("FINAL/metrics_hydro_allsp_relBMI_30.csv", header = TRUE) %>% 
  select(-X) %>% 
  filter(year != 2014) %>% 
  filter(year != 2018)

# get repeated sites (16 data points >> 14)
bmi_rep <- bmi %>% 
  filter(stationcode %in% c("901BELOLV", "901NP9FLC", "911NP9UCW", "903NP9SLR",
                            "909SWCASR", "911TJKC1x"))

bmi_rep$stationcode <- as.factor(bmi_rep$stationcode)

model.rep <- lm(hill_shannon ~  
                  # dry_date +
                  # dry_duration +
                  # peak_date +
                  # peak_depth +
                  # peak2sample_slope +
                  # recession_slope +
                  # false_starts_per_duration +
                  # first_wet_date +
                  # rewet_slope +
                  # total_rewet +
                  year, 
                  data = bmi_rep)
summary(model.rep) 
  # year p value 0.0552 
  # station code as random effect p = 0.122
check_singularity(model.rep)
  # FALSE
model_performance(model.rep)
AICc(model.rep)
anova(model.rep)
  # year p value = 0.1206


#### 5.b Algae ####
alg <- read.csv("FINAL/metrics_hydro_allsp_Algae_30.csv", header = TRUE) %>% 
  select(-X)%>% 
  filter(year != 2014) %>% 
  filter(year != 2018)

# get those sites (16 data points)
alg_rep <- alg %>% 
  filter(stationcode %in% c("901BELOLV", "901NP9FLC", "911NP9UCW", "903NP9SLR",
                            "909SWCASR", "911TJKC1x"))

alg_rep$stationcode <- as.factor(alg_rep$stationcode)

model.rep <- lmer(hill_shannon ~  
                    # dry_date +
                    # dry_duration +
                    # peak_date +
                    # peak_depth +
                    # peak2sample_slope +
                    # recession_slope +
                    # false_starts_per_duration +
                    # first_wet_date +
                    # rewet_slope +
                    # total_rewet +
                    year + (1|stationcode)
                  ,
                  data = alg_rep)
summary(model.rep) 
  # year p value 0.304
  # station code random effect year p = 0.226
check_singularity(model.rep)
  # FALSE
model_performance(model.rep)
AICc(model.rep)
anova(model.rep)
`# year p value = 0.3043


#### 5.c Diatom ####
dia <- read.csv("FINAL/metrics_hydro_allsp_Diatom_30.csv", header = TRUE)%>% 
  select(-X)%>% 
  filter(year != 2014) %>% 
  filter(year != 2018)

# get those sites (16 data points)
dia_rep <- dia %>% 
  filter(stationcode %in% c("901BELOLV", "901NP9FLC", "911NP9UCW", "903NP9SLR",
                            "909SWCASR", "911TJKC1x"))

dia_rep$stationcode <- as.factor(dia_rep$stationcode)

model.rep <- lm(hill_shannon ~  
                  # dry_date +
                  # dry_duration +
                  # peak_date +
                  # peak_depth +
                  # peak2sample_slope +
                  # recession_slope +
                  # false_starts_per_duration +
                  # first_wet_date +
                  # rewet_slope +
                  # total_rewet +
                  year,
                data = dia_rep)
summary(model.rep) 
  # year p value 0.761
  # station code as random effect 0.772
check_singularity(model.rep)
  # FALSE
model_performance(model.rep)
AICc(model.rep)
anova(model.rep)
  # year p value = 0.7610


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 6: BMI (rel abund) models ------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bmi <- read.csv("FINAL/metrics_hydro_allsp_relBMI_30.csv", header = TRUE) %>% 
  select(-X) %>% 
  filter(year != 2014) %>% 
  filter(year != 2018)

model.glo <- lm(specrich ~  
                  dry_date +
                  dry_duration +
                  peak_date +
                  peak_depth +
                  peak2sample_slope +
                  recession_slope +
                  false_starts_per_duration +
                  first_wet_date +
                  rewet_slope +
                  total_rewet,
                data = bmi)
summary(model.glo) 
# total_rewet (0.210) and year (0.0407)
check_singularity(model.glo)
# FALSE
model_performance(model.glo)

par(mfrow = c(2, 2))
plot(model.glo)
# ok.. the only bad thing is linear response

AICc(model.glo)
# 214.6341
anova(model.glo)
# year (0.04075)

# null model
model.null <- lm(hill_shannon ~ 1, data = bmi)
check_singularity(model.null)
summary(model.null) 
AICc(model.null)


# univariate models
model <- lm(hill_shannon ~  
                   dry_date 
                  # dry_duration 
                  # false_starts_per_duration
                  # first_wet_date 
                  # peak_date 
                  # peak_depth 
                  # peak2sample_slope 
                  # recession_slope 
                  # rewet_slope 
                  # total_rewet,  
                  # year,
                , data = bmi)
check_singularity(model)
par(mfrow = c(2, 2))

plot(model)
summary(model) 

AICc(model)
# anova(model.null, model.tw, model.fs)
# ## NO DIFFERENCE

# plot the "best" models
ggplot() +
  geom_point(data = bmi,
             aes(x = dry_duration,
                 y = specrich)) +
  theme_classic() +
  ylab("Species Richness (q = 0)")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 7: Algae models ----------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

alg <- read.csv("FINAL/metrics_hydro_allsp_Algae_30.csv", header = TRUE) %>% 
  select(-X) %>% 
  filter(year != 2014) %>% 
  filter(year != 2018)

model.glo <- lm(hill_shannon ~  
                  dry_date +
                  dry_duration +
                  peak_date +
                  peak_depth +
                  peak2sample_slope +
                  recession_slope +
                  false_starts_per_duration +
                  first_wet_date +
                  rewet_slope +
                  total_rewet,
                data = alg)
summary(model.glo) 
# nothing
check_singularity(model.glo)
# FALSE

par(mfrow = c(2,2))
plot(model.glo)
# spec - slight heteroscedasticity 
# h_shann - slight heteroscedasticity 
 # log / sqrt transform makes it worse...

summary(model.glo) 
model_performance(model.glo)
# R2 = 0.303 
# adj R2 = -0.341
AICc(model.glo)
anova(model.glo)

# null model
model.null <- lm(hill_shannon ~ 1, data = alg)
check_singularity(model.null)
summary(model.null) 
AICc(model.null)


# univariate models
model <- lm(hill_shannon ~  
                dry_date 
              #  dry_duration 
              #  false_starts_per_duration
              #  first_wet_date 
              # peak_date 
              #  peak_depth 
              #  peak2sample_slope 
              #  recession_slope 
              #  rewet_slope 
              #  total_rewet  
               # year
            , data = alg)
check_singularity(model)
summary(model) 
AICc(model)

par(mfrow = c(2,2))
plot(model)

anova(model.null, model.pd, model.rd)


# plot the "best" models
ggplot() +
  geom_point(data = alg,
             aes(x = peak_depth,
                 y = specrich_qual)) +
  theme_classic() +
  ylab("Species Richness (q = 0)")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 8: Diatom models ---------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dia <- read.csv("FINAL/metrics_hydro_allsp_Diatom_30.csv", header = TRUE) %>% 
  select(-X) %>% 
  filter(year != 2014) %>% 
  filter(year != 2018)

model.glo <- lm(specrich ~  
                  dry_date +
                  dry_duration +
                  peak_date +
                  peak_depth +
                  peak2sample_slope +
                  recession_slope +
                  false_starts_per_duration +
                  first_wet_date +
                  rewet_slope +
                  total_rewet,
                data = dia)
summary(model.glo) 
# almost - dry_date (0.0577) and recession slope (0.0735)
check_singularity(model.glo)
# FALSE

par(mfrow = c(2,2))
plot(model.glo)
# spec - slight heteroscedasticity 
# h_shann - slight heteroscedasticity , some outliers
# log / sqrt transform doesn't really help

model_performance(model.glo)
# R2 = 0.562  
# adj R2 = 0.158 
AICc(model.glo)
anova(model.glo)
# dry date (0.01014) and recession sleop (0.06665)


# null model
model.null <- lm(hill_shannon ~ 1, data = dia)
check_singularity(model.null)
summary(model.null) 
AICc(model.null)


# univariate models
model <- lm(hill_shannon ~  
                  dry_date 
                # dry_duration 
                # false_starts_per_duration
                # first_wet_date 
                # peak_date 
                # peak_depth 
                # peak2sample_slope 
                # recession_slope 
                # rewet_slope 
                # total_rewet,  
                # year,
            , data = dia)
check_singularity(model)
summary(model) 
AICc(model)

par(mfrow = c(2,2))
plot(model)


anova(model.null, model.p2s, model.dd, model.rd)

# plot the "best" models
ggplot() +
  geom_point(data = dia,
             aes(x = rewet_date,
                 y = specrich)) +
  theme_classic() +
  ylab("Species Richness (q = 0)")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 9: year on metrics ------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

hydro <- read.csv("FINAL/metrics_hydro_allsp_relBMI_30.csv", header = TRUE) %>% 
  select(-X) %>% 
  filter(year != 2014) %>% 
  filter(year != 2018)

hydro <- hydro[ , c(1:25)]

# univariate models
model <- lm(# dry_date 
            # dry_duration 
            # false_starts_per_duration
            # first_wet_date 
            # peak_date 
            # peak_depth 
            # peak2sample_slope 
            # recession_slope 
            # rewet_slope 
            # total_rewet
  # Alkalinity.as.CaCO3
  # pH
  # SpecificConductivity
  # Temperature
  # Canopy.Cover
  # Pool
  # Riffle
  # Run
   Wetted.Width 
            ~ year,
            , data = hydro)
check_singularity(model)

summary(model) 
AICc(model)

par(mfrow = c(2, 2))
plot(model)


