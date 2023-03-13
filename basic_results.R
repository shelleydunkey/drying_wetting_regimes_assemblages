#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: basic_results
#Coder: Michelle Busch
#Date: 20221107
#Purpose: Basic results to report for results section
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## 27 events total ##

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 1: Setup workspace -------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Clear memory
remove(list=ls())

#setting working directory
setwd("C:/Users/mhope/Desktop/Allen_Lab/RA/drying_regimes_2_github/processed_data/")

#loading libraries
require(tidyverse)
require(performance)
require(lmerTest)
require(MuMIn)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 2: BMI Data --------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load full matrix
bmi <- read.csv("FINAL/final_bmi_env_30.csv", header = TRUE) %>% 
  select(-X) %>% 
  mutate(total_rewet = peak2sample_duration + rewet_duration) %>% 
  filter(year != 2014) %>% 
  filter(year != 2018)

# q1. what are the most / least common species?
bmi_sp <- bmi[ , c(32:267)]

bmi_sp <- bmi_sp[,colSums(bmi_sp/bmi_sp, na.rm=T)>0] 
# 152 species

# change to presence / absence matrix, sum up to find most common and most rare
bmi_sp[bmi_sp > 0] <- 1

bmi_totals <- as.data.frame(colSums(bmi_sp))
bmi_totals$percent <- bmi_totals$`colSums(bmi_sp)`/27

table(bmi_totals$percent)

# species metrics
bmi_metrics <- read.csv("FINAL/metrics_hydro_allsp_relBMI_30.csv", 
                        header = TRUE) %>% 
  select(-X) %>% 
  mutate(total_rewet = peak2sample_duration + rewet_duration) %>% 
  filter(year != 2014) %>% 
  filter(year != 2018)

# q2. max / min of species richness?
summary(bmi_metrics$specrich)
metrics_year <- bmi_metrics %>% 
  select(year, specrich)

table(metrics_year$year)

# summary stats
bmi_summary <- bmi_metrics[ , c(1, 2, 3, 32:36)]

summary(bmi_summary)
sd(bmi_summary$hill_shannon, na.rm = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 3: Algae Data ------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load full matrix
alg <- read.csv("FINAL/final_Algae_env_30.csv", header = TRUE) %>% 
  select(-X) %>% 
  mutate(total_rewet = peak2sample_duration + rewet_duration) %>% 
  filter(year != 2014) %>% 
  filter(year != 2018)

# q1. what are the most / least common species?
alg_sp <- alg[ , c(32:400)]

alg_sp <- alg_sp[,colSums(alg_sp/alg_sp, na.rm=T)>0] 
# 207 species

# change to presence / absence matrix, sum up to find most common and most rare
alg_sp[alg_sp > 0] <- 1

alg_totals <- as.data.frame(colSums(alg_sp))
alg_totals$percent <- alg_totals$`colSums(alg_sp)`/27

table(alg_totals$percent)

# species metrics
alg_metrics <- read.csv("FINAL/metrics_hydro_allsp_Algae_30.csv", 
                        header = TRUE) %>% 
  select(-X) %>% 
  mutate(total_rewet = peak2sample_duration + rewet_duration) %>% 
  filter(year != 2014) %>% 
  filter(year != 2018)

# q2. max / min of species richness?
summary(alg_metrics$specrich)
metrics_year <- alg_metrics %>% 
  select(year, specrich)

table(metrics_year$year)

summary(alg_metrics)
sd(alg_metrics$specrich_qual, na.rm = TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 3: Diatom Data -----------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load full matrix
dia <- read.csv("FINAL/final_Diatom_env_30.csv", header = TRUE) %>% 
  select(-X) %>% 
  mutate(total_rewet = peak2sample_duration + rewet_duration) %>% 
  filter(year != 2014) %>% 
  filter(year != 2018)

# q1. what are the most / least common species?
dia_sp <- dia[ , c(32:407)]

dia_sp <- dia_sp[,colSums(dia_sp/dia_sp, na.rm=T)>0] 
# 225 species

# change to presence / absence matrix, sum up to find most common and most rare
dia_sp[dia_sp > 0] <- 1

dia_totals <- as.data.frame(colSums(dia_sp))
dia_totals$percent <- dia_totals$`colSums(dia_sp)`/27

table(dia_totals$percent)

# species metrics
dia_metrics <- read.csv("FINAL/metrics_hydro_allsp_Diatom_30.csv", 
                        header = TRUE) %>% 
  select(-X) %>% 
  mutate(total_rewet = peak2sample_duration + rewet_duration) %>% 
  filter(year != 2014) %>% 
  filter(year != 2018)

# q2. max / min of species richness?
summary(dia_metrics$specrich)
metrics_year <- dia_metrics %>% 
  select(year, specrich)

table(metrics_year$year)

summary(dia_metrics)
sd(dia_metrics$hill_shannon, na.rm = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 4: Hydrologic Metrics ----------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# upload a file that you know is accurate (used in other analyses)
hydro <- read.csv("FINAL/final_bmi_env_30.csv", header = TRUE) %>% 
  select(-X) %>% 
  filter(year != 2014) %>% 
  filter(year != 2018)

# filter out the variables you need
hydro <- hydro[ , c(1:29)]

hydro <- hydro %>%
  select(-c(recession_rsq, recession_coef, 
            rewet_duration, antecedent_time_period_days, 
            false_starts, sampling_date))

summary(hydro)
sd(hydro$first_wet_date, na.rm = TRUE)

# just do some plots
ggplot() +
  geom_point(data = hydro,
             aes(x = year,
                 y = Run)) +
  theme_classic() +
  xlab("Year")

# linear models to explain the variables by year?
# set 1: hydro
mod1 <- lm(dry_date ~ year, data = hydro)
mod2 <- lm(dry_duration ~ year, data = hydro)
mod3 <- lm(peak_date ~ year, data = hydro)
mod4 <- lm(peak_depth ~ year, data = hydro)
mod5 <- lm(peak2sample_slope ~ year, data = hydro)
mod6 <- lm(recession_slope ~ year, data = hydro)
mod7 <- lm(rewet_date ~ year, data = hydro)
mod8 <- lm(false_starts_per_duration ~ year, data = hydro)
mod9 <- lm(first_wet_date ~ year, data = hydro)
mod10 <- lm(rewet_slope ~ year, data = hydro)
mod11 <- lm(total_rewet ~ year, data = hydro)


summary(mod11) 
check_singularity(mod11)
model_performance(mod11)
AICc(mod11)


# set 2: environ
moda <- lm(Alkalinity.as.CaCO3 ~ year, data = hydro)
modb <- lm(pH ~ year, data = hydro)
modc <- lm(SpecificConductivity ~ year, data = hydro)
modd <- lm(Temperature ~ year, data = hydro)
mode <- lm(Canopy.Cover ~ year, data = hydro)
modf <- lm(Pool ~ year, data = hydro)
modg <- lm(Riffle ~ year, data = hydro)
modh <- lm(Run ~ year, data = hydro)
modi <- lm(Wetted.Width ~ year, data = hydro)


summary(modi) 
check_singularity(modi)
model_performance(modi)
AICc(modi)