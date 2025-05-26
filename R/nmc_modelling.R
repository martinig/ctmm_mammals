#library(ctmm)
#setwd("~/Dropbox/UBC/Grants/2023_Caribou_Collaring/year_1_report")

#Data import
data <- read.csv("nmc_tracking_data.csv")
data <- data[data$outlier !=1,]


#Get herd IDs
pink <- unique(data[data$Pop_Unit == "Pink Mountain","WLH_ID"])
muskwa <- unique(data[data$Pop_Unit == "Muskwa","WLH_ID"])
#-------------------------------------
#Function for detrending the speeds
#-------------------------------------

detrend.speeds <- function(DATA, CTMM) {
  
  # estimate speeds conditional on the data
  EST <- speeds(DATA, CTMM, fast=TRUE, level=NULL)
  # null estimates of speed (no data, only model)
  EST.NULL <- speeds(CTMM, t=DATA$t, fast=TRUE, level=NULL)
  
  # reduce DOF by that of null distribution
  DOF <- EST$DOF - EST.NULL$DOF
  all(DOF>0) # check for validity... always has held so far
  DOF <- pmax(0,DOF) # haven't needed this, but just in case
  
  S2 <- EST$speed^2 + EST$VAR
  
  AVE <- (EST$DOF*S2)/DOF
  
  # Calculate CIs
  CI <- sapply(1:length(DATA$t), function(i){ctmm:::chisq.ci(AVE[i], DOF=DOF[i], level=0.95, robust=TRUE)})
  CI <- sqrt(CI)
  
  SPEEDS <- as.data.frame(t(CI))
  SPEEDS$time <- DATA$timestamp
  
  return(SPEEDS)
}

#-------------------------------------
#Movement modelling
#-------------------------------------

#Convert to telemetry object
DATA <- as.telemetry(data)
#key things to check for/keep in mind: 
#movebank format is usually okay - you need to have the column individual.local.identifier, this is grabbed first
#converted Acquisition_Date to timestamp - the timestamp column is grabbed first
#outliers are 1's, need to be dealt with before proceeding as they mess with the model output
#there will be warnings and/or errors - those are sanity checks
#mark.rm = TRUE will remove the outliers

#summary(DATA) gives the sampling interval and coordinates if you want to plot the summary data
#plot(DATA) is a nice plot function to visualize the data 
#plot(DATA, col=rainbow(length(DATA)))
#class(DATA[[1]])
#test<-DATA[[1]]
#head(test)
#HDOP is the measurement error (horizontal dilution of precision)
#test@info #this gives you the identity if you need the know the individual dataset
#test@UERE #each collar manufacturer will give a different value for this, it converts the DOP value into an error in meters (DOP of 1 is 10 meters)
#another example for errors: test<-oulie(turtle[[3]]) #you will see drifting in the movement, this is a vague spot where you have to make subjective decisions, cross that bridge when you get to it

#test<-oulie(DATA[[4]]) #the farther an animal gets away from it's core use area is more red and larger points (red is less informative - those are distance based outliers)
#head(test)
#plot(test) #the blue lines connect the locations, they get thicker and more blue as the speed increases (you want to see connections between paths, but watch for the angles, they shouldn't be too tight) (blue are more informative - these are the speed based outliers)
#plot(test, units=FALSE) #changes the units

#fit error parameters to calibration data
UERE<-uere.fit(DATA[[5]][2134:3099,])
#do not run uere.fit on data where they are actually moving

#estimated error model parameters
summary(UERE)

#Drop mortality event
DATA[[5]] <- DATA[[5]][1:2134,] #[2134:3099,] #when there are lots of animals, consider just dropping it 

#uere.fit(DATA[[5]][2134:3099,])

FITS <- list()
for(i in 1:length(DATA)){
  GUESS <- ctmm.guess(DATA[[i]],interactive=FALSE) #pre-step is guessing the initial parameters, if you run just this line alone for i=4 and interactive=TRUE you get a variogram you can zoom in on. You want to see curvature if you want to estimate speed. If you are estimating home range, you want to see it level off. Red line is model estimate, black line is actual data
  GUESS$error <- TRUE #this is the DOP of 10 unless you do calibration (see above - lines 71 to 76)
  FITS[[i]] <- ctmm.select(DATA[[i]], GUESS, trace = 2) #trace 2 shows you the model is running in the console
  save(FITS, file = "nmc_models.Rda")
  print(i)
}

summary(FITS[[4]]) #, units = FALSE standardizes in seconds
#google anisotropic - Describes a material, property, or process that behaves differently depending on the direction in which it is measured. E.g., anisotropic movement means that the probability of moving in one direction differs from the probability of moving in another—often due to landscape features like barriers or gradients
#area is not useful, use the home range estimate below instead
#DOF is degrees of freedom, estimated based on the autocorrelation structure, e.g., DOF mean = on the centroid, DOF area = effective sample size for home range estimation, DOF diffusion = how much info you have to estimate this rate, DOF speed = how much info you have to estimate this rate - all these values could go into a meta-analysis
#your sampling frequency needs to be ~3x the tau v to get an estimate for speed
#tau p positional autocorrelational time scale (how long until it reverts back to the mean - i.e., return to the center of its home range)
#the speed estimate from here is from the fitted model, you can get a better estimate elsewhere
#diffusion (think if you drop ink in water, how long does it take to spread out) - this should be a straight line - you can often get this even if you can't get speed - this is a good variable, correlates almost perfectly with speed

names(FITS) <- names(DATA)
save(FITS, file = "nmc_models.Rda")

#Range crossing time
meta(FITS[pink], variable = "tau position") #pink is the name of the herd
meta(FITS[muskwa], variable = "tau position")
meta(list(muskwa = UDS[muskwa],
          pink = UDS[pink]),
     variable = "tau position")


#tau_v
meta(FITS[pink], variable = "tau velocity")
meta(FITS[muskwa], variable = "tau velocity")


#Diffusion rates
meta(FITS[pink], variable = "diffusion")
meta(FITS[muskwa], variable = "diffusion")
meta(list(muskwa = UDS[muskwa],
          pink = UDS[pink]),
     variable = "diffusion")

#-------------------------------------
#Home range estimation
#-------------------------------------

# create aligned UDs
UDS <- akde(DATA,FITS, weights = TRUE) #estimates the home range conditional on the fitted model, weights = TRUE use this if the sampling is uneven, it upweights and downweights the sampling based 
save(UDS, file = "nmc_home_ranges.Rda")

meta(UDS[pink])

meta(UDS[muskwa])

meta(list(muskwa = UDS[muskwa], pink = UDS[pink]), variable = "area")

#-------------------------------------
#Population range estimation
#-------------------------------------

#Get Pink Mountain pkde
pink_data <- DATA[pink]
pink_FITS <- FITS[pink]
pink_UDs <- UDS[pink]

#Estimate pkde for the pink mountain herd
pink_PKDE <- pkde(pink_data,pink_UDs) #this likely won't be needed (but it estimates the population kernel density estimate)
save(pink_PKDE, file = "pink_PKDE.Rda")

#Get Muskwa pkde
muskwa <- unique(data[data$Pop_Unit == "Muskwa","WLH_ID"])
muskwa_data <- DATA[muskwa]
muskwa_FITS <- FITS[muskwa]
muskwa_UDs <- UDS[muskwa]

#Estimate pkde for the muskwa herd
muskwa_PKDE <- pkde(muskwa_data,muskwa_UDs)
save(muskwa_PKDE, file = "muskwa_PKDE.Rda")

#-------------------------------------
#Speed estimation
#-------------------------------------

#speed is the average speed
#test<-speed(DATA[[4]], FITS[[4]]) 
#surprisingly less data = it takes longer to run
#test #this gives you the speed estimate - this is going to be very different between the output it gives at the model level

#speeds is the instantaneous movement speed at each time stamp
#test<-speeds(DATA[[4]], FITS[[4]]) 
#head(test)
#test<-speeds(DATA[[4]], FITS[[4]], level=NA) #returns the DOF and VAR for each timestep instead of 95% CIs
#plot(est ~ t, data=test) #note that the speed never reaches zero #this is why there is the detrends speed function above (on line 13 above) 
#instead of using the speeds function, use the detrend.speeds function above (it pulls all the individuals down to the same start point, i.e., 0)

SPEEDS <- list()

for(i in 1:length(DATA)){
  if(grepl("OUF",summary(FITS[[i]])$name)){
    
    #Estimate the speeds with the model detrended
    nmc_speed <- detrend.speeds(DATA[[i]], FITS[[i]])
    
    #Carpentry and export
    nmc_speed$WLH_ID <- DATA[[i]]@info$identity
    nmc_speed$Longitude  <- DATA[[i]]$longitude
    nmc_speed$Latitude  <- DATA[[i]]$latitude
    nmc_speed$temp_c <- data[as.numeric(row.names(nmc_speed)),"temp_c"]
    SPEEDS[[length(SPEEDS)+1]] <- nmc_speed
  }
  print(i)
}

SPEEDS <- do.call(rbind, SPEEDS)

write.csv(SPEEDS, file = "nmc_speeds2.csv")

#Plot
plot(est ~ temp_c, data = SPEEDS)


png(filename="/Users/michaelnoonan/Dropbox/UBC/Grants/2023_Caribou_Collaring/year_1_report/NMC_data.png",
    width = 4, height = 6, units = "in",
    res = 600)

projection(DATA) <- median(DATA)
plot(DATA, col = viridis::viridis(9))
dev.off()
