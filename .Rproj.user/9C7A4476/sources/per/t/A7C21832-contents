#force install pacman and orchaRd

#install.packages("pacman")

#load all other packages
pacman::p_load(tidyverse, # tidy family and related pacakges below
               kableExtra, # nice tables
               ctmm,  
               fitdistrplus,
               metRology,
               DAAG, 
               ggdist,
               lme4,
               mgcv,
               ggplot2, # post-hoc tests
               gridExtra, 
               pander,   # nice tables
               metafor,  # package for meta-analysis
               ape,      # phylogenetic analysis
               ggtree,
               ggstance,
               MuMIn,  # multi-model inference
               patchwork,   # putting ggplots together - you need to install via devtool
               here,    # making reading files easy
               MASS,
               cowplot #multipanel plots
)


#-------------------------------------
# Data import and carpentry
#-------------------------------------
#this is from speed_analysis.R script

data <- read.csv(here("data", "nmc_tracking_data.csv"))
data$WLH_ID <- as.factor(data$WLH_ID)

#Define when they are moving/resting
data$Active <- 0
data[which(data$est > 0),"Active"] <- 1



#other dataset (not used yet)
data2 <- read.csv(here("data", "nmc_speeds2.csv"))