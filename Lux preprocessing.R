## FlightR share code
## updated 5/5/20
## N.Cook

rm(list=ls())

# General_variables -------------------------------------------------------
# Enter the values here to name and draw from the correct folders
# load packages
# library(lubridate)

## Suggestions
## Put all geolocator files for each individual in a seperate folder titled with unique name (ie Species_Metal band number)

## May have to run several checks or individually check geolocator files for enough data to make it worth running

##choose calibration period as: before geolocator placed on bird (rooftop calibration) or after geolocator on bird
cal.period <- "b" # b- before on bird; a - after 

## Starting and ending dates for breeding and molting windows
breedStart <- "-05-01"
breedEnd <- "-06-15"
moltStart <- "-08-15"
moltEnd <- "-09-15"

##choose species (used in folder naming)
species <- "CANV" # Choose between CANV, MALL and WODU

## Set maindDir to folder with individual folders
mainDir <- file.path("lux_files")

#setwd(mainDir)

## Read in CSV file with summary information
#Columns needed Species, metal band number, sex, Turned.on.date,Date.of.Initial.fitting,
#shot.recovered.date,last.record.date,Capture.Lat,Capture.Lon
geo_df <- read.csv("CMW_geo.csv",stringsAsFactors = FALSE)

## List all the files or folders for individuals for which you want to run geolocator analysis
fileList <- list.files(mainDir)

## Creates index for matching geo_df to folders with geolocator files
index.bird <- which(fileList %in% paste0(geo_df$Species,"_",geo_df$Metal))
birdList <- fileList[index.bird]

## Check that geo_df and birdList have the same order - It just makes it easier for assumptions made later using unique ID
i=1

# Light_setup_process -----------------------------------------------------
# Runs the first loop to get setup data for each bird
{
  require(TwGeos)     # NOTE: must be installed from GitHub: devtools::install_github("SLisovski/TwGeos")
  for (i in 1:nrow(geo_df)) {
    #create file paths
    birdName <-fileList[index.bird[i]]
    birdFold <- file.path(mainDir,birdName)
    lightDir <- file.path(birdFold,"light")
    #Add folder for lux manipulation if it doesn't already exist
    if (!dir.exists(lightDir)) {dir.create(file.path(lightDir))}
    #read in Lux file data
    LUXfiles <- grep("lux",list.files(birdFold))
    LUXdata <- read.table(file.path(birdFold,list.files(birdFold)[LUXfiles[length(LUXfiles)]]),
                          skip=21,col.names=c("Date","Light"), sep = "\t",stringsAsFactors = FALSE)
    
    # format as a POSIX date-time -- Especially look for date matching day, month, year order and - or / seperators
    LUXdata$Date <- as.POSIXct(strptime(LUXdata$Date,format="%d/%m/%Y %H:%M:%S",tz="GMT"))
    
    # Set the lowest light value to 1.
    LUX.min <- min(LUXdata$Light)
    LUXdata$Light[LUXdata$Light == LUX.min] <- 1
    
    # Take the log of light data # 
    LUXdata$Light <- log(LUXdata$Light)
    threshold <- 1
    
    # Sets the capture and recapture locations based on spreadsheet
    cap.lat <- geo_df$Capture.Lat[i]
    cap.long <- geo_df$Capture.Lon[i]
    ### recap.lat <- geo_df$Recapture.Lat[i]
    ### recap.long <- geo_df$Recapture.Lon[i]
    ### capRecapDiff <- !(cap.lat == recap.lat & cap.long == recap.long)
    ###capRecapDiff <- TRUE # If the calculated difference is not working then use this line
    
    # Create a vector of the dates the geolocator was recording data
    tm <- seq(LUXdata$Date[1], LUXdata$Date[nrow(LUXdata)], by = "day")
    # write.csv(tm,"tm.csv",row.names = FALSE)
    
    # Create a vector of TRUE /FALSE to determine whether it's sunrise or sunset
    rise <- rep(c(TRUE, FALSE), length(tm))
    # making predicted twilight times at capture location
    c.dat <- data.frame(Twilight = SGAT::twilight(rep(tm, each = 2),
                                                  lon = cap.long,
                                                  lat = cap.lat,
                                                  rise = rise, zenith = 96),
                        Rise = rise)
    # set offset for plotting - this is just a plotting parameter
    offset=19
    
    #create plot to visualize data before preprocessLight function
    lightImage(LUXdata, offset = offset, zlim = c(0,12))
    tsimagePoints(c.dat$Twilight, offset = offset,
                  pch = 16, cex = 0.5,
                  col = ifelse(c.dat$Rise, "dodgerblue", "firebrick"))
    # adds line at two equinoxes. Dates run 2015 to 2018
    # Change the dates if necessary (can vary by year). This is for illustration only. Has no effect on analysis
    eqnx<-as.POSIXct(c("2015-3-20","2015-9-23","2016-3-20","2016-9-22","2017-3-20","2017-9-22","2018-3-20","2018-9-23"), tz = "GMT")
    abline(v = eqnx, lwd=3, lty=3, col="purple")
    abline(v = as.POSIXct(strptime(geo_df$Date.of.Initial.fitting[i],"%Y-%m-%d",tz="GMT")),lwd=1, lty=3, col="yellow")
    abline(v = as.POSIXct(strptime(geo_df$shot.recovered.date[i],"%Y-%m-%d",tz="GMT")),lwd=1, lty=3, col="yellow")
    title(paste(birdName," (",geo_df$Sex[i],") LUX plot",sep = ""))
    
    # Preprocess light with 4 step process described well in R documentation for function
    twl <- preprocessLight(LUXdata, threshold, offset = offset, zlim = c(0, 12))
    
    # Save LUX plot with markers for specific dates
    #       -yellow dashed dates for geo on and off bird
    #       -green solid for start and end dates of data used
    #       -purple dashed for equinoxes
    pdf(file.path(lightDir,"Add_dates_Lux_plot.pdf"),title = paste(birdFold,"LUX plot"))
    lightImage(LUXdata, offset = offset, zlim = c(0,12))
    tsimagePoints(c.dat$Twilight, offset = offset,
                  pch = 16, cex = 0.5,
                  col = ifelse(c.dat$Rise, "dodgerblue", "firebrick"))
    # adds line at two equinoxes. Dates run 2015 to 2018
    # Change the dates if necessary (can vary by year). This is for illustration only. Has no effect on analysis
    eqnx<-as.POSIXct(c("2015-3-20","2015-9-23","2016-3-20","2016-9-22","2017-3-20","2017-9-22","2018-3-20","2018-9-23"), tz = "GMT")
    abline(v = eqnx, lwd=3, lty=3, col="purple")
    useDates <- as.POSIXct(c(twl$Twilight[1],twl$Twilight[nrow(twl)]),tz = "GMT")
    abline(v = useDates,lwd=1, col="green")
    abline(v = as.POSIXct(strptime(geo_df$Date.of.Initial.fitting[i],"%Y-%m-%d",tz="GMT")),lwd=1, lty=3, col="yellow")
    abline(v = as.POSIXct(strptime(geo_df$shot.recovered.date[i],"%Y-%m-%d",tz="GMT")),lwd=1, lty=3, col="yellow")
    title(paste(birdName," (",geo_df$Sex[i],") LUX plot",sep = ""))
    dev.off()
    
    # function to check for and fill missing days in twilight dataframe (twl)---
    fill.twl.gaps <- function (twl) {
      # Get dates only from twilight times
      twlCheck <- as.Date(twl$Twilight)
      # get sunrise and sunset dates
      rise <- twlCheck[twl$Rise]
      set <- twlCheck[twl$Rise==FALSE]
      
      # if sunrise and sunset dates match then remove the duplicate sunrise and sunset dates
      if (all(set == rise)) {
        rise <- rise[!diff(rise) == 0]
        set <- set[!diff(set) == 0]
      } else {
        
      }
      
      # Determine if sunrise or sunset is first
      if (twl$Twilight[!twl$Rise][1] == twl$Twilight[1]) {
        ti <- set
        riseCol <- c(FALSE,TRUE)
      } else {
        stop("Sunset is not the first twilight time according to GMT based dates")
      }
      
      # Check for missing dates
      if (!all(diff(ti) <= 1)) {
        # Add the missing dates
        rn <- which(diff(ti) > 1)
        diff <- abs(difftime(ti[rn],ti[rn+1],units = "days"))
        #set up dates for missing days
        for (k in 1:length(rn)) {
          twlAdd <- seq(ti[rn][k]+1,ti[rn+1][k]-1,1)
          twlAdd <- as.POSIXct(strptime(twlAdd,format = "%Y-%m-%d", tz="GMT"))
          # based on sunrise or sunset first set the generic times
          twlAdd <- sort(rep(twlAdd,2))
          twlAdd <- twlAdd + c(3*3600,11*3600)
          lenAdd <- length(twlAdd)
          twlAdj <- data.frame(Twilight=twlAdd,
                               Rise=rep(riseCol,lenAdd/2),
                               Deleted=rep(TRUE,lenAdd),
                               Marker=0,
                               Inserted=rep(TRUE,lenAdd),
                               Twilight3=twlAdd,
                               Marker3=0)
          
          twl <- rbind(twl,twlAdj)
          twl <- twl[order(twl$Twilight),]
        }
      }
      
      return(twl)
    }
    
    # Run twilight gap fill function then record which lines were added to pass to twl after Automatic twilight editing
    twl <- fill.twl.gaps(twl)
    delCol <- which(twl$Deleted==TRUE)
    #save as csv file and workspace so it can be reused
    write.csv(twl, file = file.path(lightDir,"twlsampleLUX.csv"))
    # save files for further analysis and troubleshooting code errors
    save(tm,LUXdata,twl,threshold,cap.long,cap.lat,
         file = file.path(lightDir,"Initial_values.RData"))
    
    save.image(file = file.path(lightDir,"Initial_setup.RData"))
    rm(list=setdiff(ls(),c("mainDir","spDir","geo_df","fileList","index.bird", "lcheck",
                           "ptm","l.correct","cal.period", "species","recheck","errs")))
  }
  
  if (length(errs) > 0) {
    write(errs,paste("Light_Errors",Sys.Date(),".txt"))
    print(paste("These birds had errors",errs,sep = ""))
  } else {
    print("No errors recorded")
  }
}


# Automatic removal of worst twilight slopes--------------
run.twilight.lowLight.removal <- function(birdFold, speciesWD=mainDir,
                                          firstStep=TRUE, offset=19,window=4,outlier.mins=45,
                                          stationary.mins=15, zThreshold=3) {
  require(TwGeos)
  require(FLightR)
  runerror <- c()
  lightDir <- file.path(birdFold,"light")
  flPrep <- file.path(speciesWD,birdFold,"FLightR_prep.csv")
  
  oad(file.path(lightDir,"Initial_values.RData"))
  
  # First step twilight editing - TwGeos comparison to neighbors
  delCol <- which(twl$Deleted)
  if (firstStep) {
    tryCatch({
      twlAuto <- twilightEdit(twl,
                              offset = offset,
                              window = window, # two days before and two days after
                              outlier.mins = outlier.mins, # difference in minutes to determine an outlier
                              stationary.mins = stationary.mins) # are the surrounding twilights within 25 mins of one another
      twlAuto$Deleted[delCol] <- TRUE
      write.csv(twlAuto, file = file.path(lightDir,"twlsampleLUX_autoedit.csv"))
      dev.copy2pdf(file = file.path(lightDir,"editedLuxplot.pdf"))
      dev.off()
    }, warning = function(warn) {
      runerror <- c(runerror,paste(birdFold,warn))
    }, error = function(err) {
      runerror <- c(runerror,paste(birdFold,err))
    })
  }
  
  ## convert and write out the date to csv for FlightR format
  TAGS.twilight.raw <- twGeos2TAGS(LUXdata,twlAuto,threshold)
  # TAGS.twilight.raw <- BAStag2TAGS(LUXdata,twl,threshold)
  TAGS.twilight.raw$datetime <- format(TAGS.twilight.raw$datetime, format = "%Y-%m-%d%T.000Z")
  write.csv(TAGS.twilight.raw, file= flPrep, quote=FALSE, row.names = FALSE)
  
  ## Read in csv and set some of the variables for recording on geolocator
  Proc.data <- get.tags.data(flPrep, saves="max", measurement.period = 60)
  
  pDa <- Proc.data$Twilight.log.light.mat.dawn[26:49,]
  pDaDF <- as.data.frame(pDa)
  pDaMax <- apply(pDa,2,max)
  pDaMin <- apply(pDa,2,min)
  
  pDu <- Proc.data$Twilight.log.light.mat.dusk[1:24,]
  pDuDF <- as.data.frame(pDu)
  pDuMax <- apply(pDu,2,max)
  pDuMin <- apply(pDu,2,min)
  
  #get the twilights with more zeros than the zthreshold
  d1 <- c(apply(pDa,2,function(x) length(which(x == 0)) > zThreshold ),
          apply(pDu,2,function(x) length(which(x == 0)) > zThreshold ))
  d2 <- rep(FALSE,length(d1))
  # get which twilights are not excluded under either condition
  k3 <- !apply(data.frame(d1=d1,d2=d2),1,function(x) any(x))
  
  # create data frame of datetime of twilight in one column and a vector of keep logical values in the other
  kdf <- data.frame(twlTime=c(Proc.data$Twilight.time.mat.dawn[25,],Proc.data$Twilight.time.mat.dusk[25,]),
                    keep=k3)
  # change date time values from number to POSIX values
  kdf$time <- as.POSIXct(kdf$twlTime,origin = "1970-01-01",tz="GMT")
  # reorder data frame by date time
  kdf <- kdf[order(kdf$twlTime),]
  
  # number of twilights that should still be marked as undeleted after additional ones marked
  length(twlAuto$Deleted[twlAuto$Deleted == FALSE]) -
    length(which(!kdf$keep))
  
  twlAuto$Deleted[twlAuto$Deleted == FALSE][which(!kdf$keep)] <- TRUE
  
  save(tm,LUXdata,twl,twlAuto,threshold,cap.long,cap.lat,recap.long,recap.lat,capRecapDiff,
       file = file.path(lightDir,"Initial_valuesEdited.RData"))
  file.remove(flPrep)
  
  return(runerror)
}


#Run the function to remove low light twilights
birdList <- fileList[index.bird]
runError <- lapply(birdList, run.twilight.lowLight.removal,offset = 19,
                   window = 8,outlier.mins = 30,stationary.mins = 10,
                   zThreshold=3 #more than zThreshold quantity of 0 values in light data before sunset or after sunrise are marked as excluded
)
#checks for and writes out file with error messages
ulError <- unlist(runError)
if(!is.null(ulError)) {
  write(ulError,file.path(mainDir,paste0(species,"_removeTwilights_runErrors_",Sys.Date(),".txt")))
} else {
  print("No error messages were saved to file")
}

