## FlightR share code
## updated 5/5/20
## N.Cook


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

setwd(mainDir)

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


# Light_setup_process -----------------------------------------------------
# Runs the first loop to get setup data for each bird
{
  require(TwGeos)
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

## Single analysis of movement for FLightR ------------------
#A google key is now required for using google maps as the base map
#Register for a key (for free) and enter as the last value in the function

#Errors can occur in different parts of anaylsis so each subsection can be run seperately if needed
#Calibration must run correctly first, then the particle filter and finally the mapping
#To run calibration set runCalib to true, to run Particle filter set runPF to true, to run map creation set runMap to true
#calib.days is length of calibration after geo on bird and stopDur is minimum duration of stopovers to be used for calibration
run.FLightR.movement <- function(birdFold, passNum, birdList=fileList[index.bird],
                                 geolocatorDF=geo_df, species=unique(geo_df$Species),calp=cal.period,
                                 mainWD=mainDir, speciesWD=mainDir,
                                 nestingTimes=nestTimes, calib.days=7, stopDur=14, stopPeriodsExcluded=NULL, 
                                 nParticles=1e4,runCalib=TRUE,runPF=TRUE,runMap=TRUE,Threads=-1, google_key="") {
  require(FLightR)
  require(grid)
  require(lubridate)
  require(ggmap)
  register_google(key=google_key)
  
  # Error and warning vector to output all the warnings and errors experienced
  runError <- c()
  
  if (runCalib) {
    tryCatch({
      birdNum <- which(birdList %in% birdFold)
      lightDir <- file.path(birdFold,"light")
      moveDir <- file.path(speciesWD,birdFold,"movement_data")
      passFold <- paste0("Movement",passNum)
      passDir <- file.path(moveDir,passFold)
      
      # create sub folder if it doesn't exist already
      if (!dir.exists(moveDir)) {dir.create(moveDir)}
      #For later passes of movement analysis this loads previous passes' data
      if(file.exists(file.path(lightDir,"Initial_valuesEdited.RData"))) {
        load(file.path(lightDir,"Initial_valuesEdited.RData"))
        cat("Edited values loaded")
      } else {
        load(file.path(lightDir,"Initial_values.RData"))
        cat("Unedited values loaded")
      }
      
      # sets the folders based on pass number of movement analysis
      if(!dir.exists(passDir)) {dir.create(passDir)}
      
      if (passNum == 1) {
        # Set calibration period to when turned on until day before deployment
        # Full record from day of deployment to earlier of last recorded day or date recovered
        on.start <- as.POSIXct(strptime(geolocatorDF$Date.of.Initial.fitting[birdNum], format = "%Y-%m-%d",tz = "GMT"))+86400
        
        # Set the calibration period as the days before geo was put on bird, or set it as the
        # calibration days (calib.days) starting when put on bird
        if (calp == "b") {
          if (geolocatorDF$Turned.on.date[birdNum] == geolocatorDF$Date.of.Initial.fitting[birdNum]) {
            print("There is no calibration periods before attachment")
          } else {
            calib.start <- as.POSIXct(strptime(geolocatorDF$Turned.on.date[birdNum], format = "%Y-%m-%d",tz = "GMT"))
            calib.end <- on.start - 86400
          }
        } else if (calp == "a") {
          calib.start <- as.POSIXct(strptime(geolocatorDF$Date.of.Initial.fitting[birdNum], format = "%Y-%m-%d",tz = "GMT"))
          calib.end <- on.start + (calib.days * 86400) # First number is days and 86400 is seconds per day
        }
        if (is.na(geolocatorDF$shot.recovered.date[birdNum])) {
          # If there is no shot or recovery date then place it as the last day of light data
          shot.date <- as.POSIXct(twlAuto[nrow(twlAuto),1])
        } else {
          shot.date <- as.POSIXct(strptime(geolocatorDF$shot.recovered.date[birdNum],format = "%Y-%m-%d",tz = "GMT"))
        }
        
        lastRec.date <- as.POSIXct(twlAuto[nrow(twlAuto),1])
        on.end <- min(shot.date, lastRec.date) - (1 * 86400) # 1 day * 86400 seconds/day
        
        if (twlAuto[1,1] > on.start) {
          calib.start <- twlAuto[1,1] + 86400
          calib.end <- calib.start + (calib.days * 86400)
          runError <- c(runError,paste("warning on ",birdFold," (birdNum=",birdNum,") calibration period after first days of banding/still run",sep = ""))
        }
        
        ## convert and write out the date to csv for FlightR format
        TAGS.twilight.raw <- twGeos2TAGS(LUXdata,twlAuto,threshold)
        TAGS.twilight.raw$datetime <- format(TAGS.twilight.raw$datetime, format = "%Y-%m-%d%T.000Z")
        write.csv(TAGS.twilight.raw, file= file.path(moveDir,"FlightR_prep.csv"), quote=FALSE, row.names = FALSE)
        
        
        ## Read in csv and set some of the variables for recording on geolocator
        Proc.data <- get.tags.data(file.path(moveDir,"FlightR_prep.csv"), saves="max", measurement.period = 60)
        
        ## Calibration
        calibration.periods <- data.frame(
          calibration.start = calib.start,
          calibration.stop = calib.end,
          lon = cap.long, lat = cap.lat)
      } else if (passNum > 1) {
        prePassNum <- passNum - 1
        prepassDir <- file.path(moveDir,paste0("Movement",prePassNum))
        preStopoverDir <- file.path(prepassDir,paste0("stopover_",prePassNum,"analysis"))
        
        # gets variables from previous analysis - From selected variables saved to RData file
        combSoFile <- file.path(preStopoverDir,paste0("CombStopover",prePassNum,"_values",".RData"))
        if (file.exists(combSoFile)) {
          load(file.path(preStopoverDir,paste0("CombStopover",prePassNum,"_values",".RData")))
          cat("Combined stopovers loaded")
        } else {
          load(file.path(preStopoverDir,paste0("Stopover",prePassNum,"_values",".RData")))
          cat("Unombined stopovers loaded")
        }
        
        #load stopover file based on 14 or 21 days minimum duration
        if (!file.exists(combSoFile)){
          subStopovers <- do.call("match",args=list(eval(parse(text=paste0("stopoversNC",stopDur))),stopovers$stopNum))
          # remove na coordinates
          subStopovers <- subStopovers[sapply(outputCoords,FUN = function(x) any(!is.na(x)))]
          #If not loaded currently load the stopover periods saved previously
          if (!is.null(stopPeriodsExcluded)) {
            spExc <- stopPeriodsExcluded[[which(names(stopPeriodsExcluded) %in% birdFold)]]
            subStopovers <- subStopovers[!(subStopovers %in% spExc)]
          }
        } else {
          subStopovers <- 1:nrow(stopovers)
        }
        
        # Error check to verify stopover periods used for calibration exist
        if (length(subStopovers) == 0) {stop("There are no calibration periods - cannot calculate further")}
        #gather coordinates for using stopovers for calibration
        lons <- sapply(outputCoords,'[[',1)
        lats <- sapply(outputCoords,'[[',2)
        #remove NA from coordinates
        lons <- lons[!is.na(lons)]
        lats <- lats[!is.na(lats)]
        
        ## Calibration
        calibration.periods <- data.frame(
          calibration.start =  stopovers$arrival[subStopovers],
          calibration.stop = stopovers$departure[subStopovers],
          lon = lons[subStopovers],
          lat = lats[subStopovers]
          )
        }
      
      ## Calibration object compiling
      calibration <- make.calibration(Proc.data,calibration.periods)
      
      ##  check all errors based on a single location and create figure
      pdf(file.path(passDir,"Errors_1_location.pdf"))
      plot_slopes_by_location(Proc.data, location = c(cap.long,cap.lat))
      abline(v=calibration.periods$calibration.start) # start of calibration period
      abline(v=calibration.periods$calibration.end) # end of calibration period
      dev.off()
      
      ## Assign spatial extent for possible locations
      ## left,right, bottom and top are the extents in degrees that birds are allowed to go
      ## distance from land use and stay are in km from shoreline birds allowed to occur/stay
      ## each has a minimum and maximum
      Grid<-make.grid(left=-179, bottom=0, right=-80, top=85,
                      distance.from.land.allowed.to.use=c(-Inf, 15),
                      distance.from.land.allowed.to.stay=c(-Inf, 3))
      
      # save values for later analysis without rerun
      save(list=ls(),file=file.path(speciesWD,birdFold,paste0(passFold,"_Calib.RData")))
      save(nestingDates,Proc.data,Grid,calibration,known.last,last.loc,cap.lat,cap.long,passFold,passDir,birdNum,
           file=file.path(speciesWD,birdFold,paste0(passFold,"_CalibValues.RData")))
    }, error = function(err) {
      runError <- c(runError,paste("error on ",birdFold," (birdNum=",birdNum,")",err,sep = ""))
      print(paste("MY_ERROR: ",birdFold,err))
    })
  }
  #Run particle filter portion
  if (runPF) {
    tryCatch({
      birdNum <- which(birdList %in% birdFold)
      if (!exists("Grid")) {
        load(file=file.path(speciesWD,birdFold,paste0("Movement",passNum,"_CalibValues.RData")))
      }
      
      birdNum <- which(birdList %in% birdFold)
      lightDir <- file.path(birdFold,"light")
      moveDir <- file.path(speciesWD,birdFold,"movement_data")
      passFold <- paste0("Movement",passNum)
      passDir <- file.path(moveDir,passFold)
      
      ## Prepare the model for run
      all.in<-make.prerun.object(Proc.data, Grid, start=c(cap.long, cap.lat), end = last.loc,
                                 Calibration=calibration,threads = Threads)
      ## Particle filter run
      # Do outlier check on 3rd pass
      if (passNum < 3) {
        check.outliers <- FALSE
      } else {
        check.outliers <- TRUE
      }
      
      Result<-run.particle.filter(all.in, threads=Threads,
                                  nParticles=nParticles, known.last=known.last,
                                  precision.sd=25, check.outliers=check.outliers)
      #Save model and partical filter results to file
      save(all.in,Result, file=file.path(passDir,paste0("Result_outliers",passNum,".RData")))
      # look at stationary periods and locations (ie stopovers, breeding, molting)
      statMigSum <- stationary.migration.summary(Result, prob.cutoff = 0.1, min.stay = 3)
      #Save model, particle filter and stationary periods to file
      save(all.in,Result,statMigSum,
           file=file.path(speciesWD,birdFold,paste0(passFold,"_statValues.RData")))
    }, warning = function(warn) {
      runError <- c(runError,paste("warning on ",birdFold," (birdNum=",birdNum,")",warn,sep = ""))
      print(paste("MY_WARNING: ",birdFold,warn))
    }, error = function(err) {
      runError <- c(runError,paste("error on ",birdFold," (birdNum=",birdNum,")",err,sep = ""))
      print(paste("MY_ERROR: ",birdFold,err))
      save(all.in,Result,
           file=file.path(speciesWD,birdFold,paste0(passFold,"_statValuesNS.RData")))
    })
    tryCatch({
      save(Proc.data,known.last,last.loc,cap.long,cap.lat,statMigSum,Result,nParticles,
           file=file.path(passDir,paste0(passFold,"_values.RData")))
    }, error = function(err) {
      runError <- c(runError,paste("error on",birdFold,"(birdNum=",birdNum,")",err))
      print(paste("MY_ERROR: ",birdFold,err))
    })
  }
  #Run the mapping portion
  if (runMap) {
    tryCatch({
      birdNum <- which(birdList %in% birdFold)
      if (!exists("statMigSum")) {
        load(file=file.path(speciesWD,birdFold,paste0("Movement",passNum,"_CalibValues.RData")))
        sval <- file.path(speciesWD,birdFold,paste0("Movement",passNum,"_statValues.RData"))
        if (file.exists(sval)) {
          load(file=sval)
        } else {
          load(file=file.path(passDir,paste0("Result_outliers",passNum,".RData")))
          print(paste("missing statMigSum for",birdFold))
        }
      }
      
      birdNum <- which(birdList %in% birdFold)
      lightDir <- file.path(birdFold,"light")
      moveDir <- file.path(speciesWD,birdFold,"movement_data")
      passFold <- paste0("Movement",passNum)
      passDir <- file.path(moveDir,passFold)
      
      ## Mapping of the results
      mapObj <- map.FLightR.ggmap(Result,return.ggobj = TRUE,
                                  save.options = list(filename=file.path(passDir,
                                                                         paste0("FR.",birdFold,"_",passNum,"map.pdf"))))
      #Add a title
      mapObj <- mapObj + ggplot2::ggtitle(paste(birdFold," (",geolocatorDF$Sex[birdNum],") Locations plot",sep = ""))
      
      # plot map as pdf
      pdf(file = file.path(passDir,paste0(birdFold,".",passNum,"map.pdf")))
      plot(mapObj)
      dev.off()
      
      #plot lat/lon graph as pdf
      pdf(file.path(passDir,paste0("lonlat",passNum,".pdf")))
      plot_lon_lat(Result)
      dev.off()
      
      # # Plot space utilization for specific range
      # plot_util_distr(Result,
      #                 dates = data.frame(tm[1],tm[7]),
      #                 add.scale.bar=TRUE, percentiles=0.5)
      
      ## more specific plotting that I am working on
      # stopover <- 5
      # sRow <- statMigSum$Potential_stat_periods[stopover,1]
      # eRow <- statMigSum$Potential_stat_periods[stopover,2]
      
      # UTIL <- plot_util_distr(Result,
      # dates = data.frame(Proc.data$FLightR.data$twilights[sRow,1],Proc.data$FLightR.data$twilights[eRow,1]),
      # add.scale.bar = TRUE, percentiles = 0.5, zoom = 5)
    },error = function(err) {
      runError <- c(runError,paste("error on",birdFold,"(birdNum=",birdNum,")",err))
      print(paste("MY_ERROR: ",birdFold,err))
    })
  }
  if (runMap) {
    dt <- ""
  } else if (runPF) {
    dt <- "cp"
  } else if (runCalib) {
    dt <- "c"
  }
  ls1 <- ls(envir=.GlobalEnv)
  save(list=ls(),file=file.path(passDir,paste0(passFold,"analysis",dt,".RData")))
  return(runError)
}




# Single stopover code -----------------
# set threholds of distance (km) and bearing (degrees) that are used to determine distinct stopover locations
run.stopover.review <- function(birdFold, passNum, threshDistance=150,birdList=fileList[index.bird],
                                mainWD=mainDir, speciesWD=mainDir,
                                geolocatorDF=geo_df, species=NULL) {
  require(lubridate)
  runError <- c()
  tryCatch({
    #pull species from geolocator data frame
    if (is.null(species)) {
      species <- unique(geo_df$Species)
      print("species pulled from geolocator data frame")
    }
    birdNum <- which(birdList %in% birdFold)
    moveDir <- file.path(speciesWD,birdFold,"movement_data")
    passFold <- paste0("Movement",passNum)
    passDir <- file.path(moveDir,passFold)
    #names and creates stopover data folder if it doesn't exist already
    stopoverDir <- file.path(passDir,paste0("stopover_",passNum,"analysis"))
    if(!dir.exists(stopoverDir)){dir.create(stopoverDir)}
    
    # gets variables from previous analysis - From selected variables saved to RData file
    load(file.path(passDir,paste0(passFold,"_values.RData")))
    
    ## Setup a dataframe of stationary period information
    statPeriod <- statMigSum$Stationary.periods
    # get labels for the stopovers
    nvn <- c(letters,as.character(c(0:9)))
    nvn <- nvn[1:nrow(statPeriod)]
    quants <- Result$Results$Quantiles
    pStatPeriod <- statMigSum$Potential_stat_periods
    pStatPeriod <- cbind(stopoverID=nvn[1:nrow(pStatPeriod)],pStatPeriod)
    pMovePeriod <- statMigSum$Potential_movement_periods
    if(pStatPeriod$start[1] == 0) {
      pStatPeriod$startDate <- c(quants$time[1] - 86400,quants$time[pStatPeriod$start])
    } else if(pStatPeriod$start[1] > 0) {
      pStatPeriod$startDate <- quants$time[pStatPeriod$start]
    }
    pStatPeriod$endDate <- quants$time[pStatPeriod$end]
    pStatPeriod$MedstartDate <- c(quants$time[1] - 86400,statPeriod$Arrival.Q.50[2:length(pStatPeriod$start)])
    pStatPeriod$MedendDate <- statPeriod$Departure.Q.50[1:length(pStatPeriod$end)]
    
    # Add columns for checks for distance
    naColumn <- rep(NA,nrow(statPeriod))
    statPeriod <- cbind(location = nvn, newStopover = naColumn,
                        overDistance = naColumn, distance = statPeriod$Distance.Moved, statPeriod)
    # Mark the columns for those stopovers that are over distance threshold
    # also mark the distinct stopovers
    statPeriod$overDistance <- statPeriod$Distance.Moved > threshDistance
    statPeriod$newStopover <- apply(cbind(statPeriod$overDistance,statPeriod$overBearing),1,all)
    
    ## Plot the stopover locations
    # get length of stopover periods
    Stoplength <- pStatPeriod$Duration/2
    lengthGT14 <- Stoplength > 14
    lengthGT21 <- Stoplength > 21
    # add logical columns marked true for stopover durations longer than 2 weeks and 3 weeks
    pStatPeriod$lengthGT14 <- lengthGT14
    pStatPeriod$lengthGT21 <- lengthGT21
    
    stopoversNC14 <- which(lengthGT14)
    stopoversNC21 <- which(lengthGT21)
    
    initStops <- stopoversNC14
    
    # Plot the timeline with median start and end times for stopover locations
    lylim <- min(statPeriod$Distance.Moved)
    uylim <- max(statPeriod$Distance.Moved)
    #upper plot height (uph) gives the y value for plotting stopover IDs
    uph <- 500*floor(uylim/500)
    #create PDF of plot
    pdf(file=file.path(stopoverDir,paste0("stopoverTiming",passNum,".pdf")))
    plot(quants$time,rep(0,length(quants$time)),
         ylim=c(lylim,uylim), type="n", main = paste0(birdFold,"_stopoverTimeline"))
    # nestMatch <- names(nestTimes) %in% birdFold2
    if (!is.null(nestingDates)) {
      abline(v=as.POSIXct(nestingDates,tm="GMT"),col="purple",lwd=4)
    }
    rect(xleft = pStatPeriod$startDate, xright = pStatPeriod$endDate,
         ybottom = lylim - 1,  ytop = uylim + 1, col="green")
    
    # plot bars above the stopover periods that are blue for over 14 days and orange for over 21 days
    if(length(stopoversNC14) > 0) {
      rect(xleft = pStatPeriod$startDate[lengthGT14], xright = pStatPeriod$endDate[lengthGT14],
           ybottom = uylim + 1,  ytop = uylim + 30, col="blue")
    }
    if(length(stopoversNC21) > 0) {
      rect(xleft = pStatPeriod$startDate[lengthGT21], xright = pStatPeriod$endDate[lengthGT21],
           ybottom = uylim + 1,  ytop = uylim + 30, col="orange")
    }
    
    # get the date/time that is centered between mean start and end dates of stopovers
    centerDates <- pStatPeriod$startDate + (pStatPeriod$endDate - pStatPeriod$startDate)/2
    # Plot the stopover location names on the stopover locations rectangles
    points(centerDates,rep(uph,length(centerDates)),pch=nvn[1:(length(nvn)-1)])
    # Add the put on and take off dates
    abline(v=as.POSIXct(geolocatorDF$Date.of.Initial.fitting[birdNum],tz="GMT"),col="orange",lwd=2)
    # Plot the distances moved from previous stopover location to current location
    # mark threshold for considering distinct stopover locations
    points(centerDates,statPeriod$Distance.Moved[1:length(centerDates)],col="red",pch=19)
    abline(h=150,col="blue")
    dev.off()
    
    ### Get the estimated locations from stopover periods
    soNum <- nrow(pStatPeriod)
    stopNums <- c(1:soNum)
    stopovers <- data.frame(stopNum=stopNums,
                            meanlon=statPeriod$Meanlon[1:soNum],meanlat=statPeriod$Meanlat[1:soNum],
                            arrival=pStatPeriod$startDate,departure=pStatPeriod$endDate)[initStops,]
    #create list of stopovers
    ll <- split(stopovers,seq(nrow(stopovers)))
    
    getStopoverCoords <- function(ll, stopNums, Result) {
      tryCatch({
        so <- which(stopNums %in% ll$stopNum)
        ## get stopover coordinates from FLightR particle filter results
        {spG <- Result$Spatial$Grid
          
          twlSeq <- which(Result$Results$Quantiles$time >= ll$arrival &
                            Result$Results$Quantiles$time <= ll$departure)
          prtPost <- Result$Results$Points.rle[twlSeq]
          # get all coordinates
          prtLoc <- spG[unlist(sapply(prtPost,"[[",2)),]
          # get lengths
          prtLen <- unlist(sapply(prtPost,"[[",1))
          #mean lon
          mulon <- sum(prtLoc[,1] * prtLen)/(nParticles*length(twlSeq))
          #mean lat
          mulat <- sum(prtLoc[,2] * prtLen)/(nParticles*length(twlSeq))
          #mean coordinates
          muCoords <- c(mulon,mulat)
        }
        
        # look at log(slope) error plots for each calibration period
        #Function taken from FLightR named plot_slopes_by_location and returns parameter data as well as plotting
        plot_slopes_by_location_adjusted <- function (Proc.data, location = c(calibration.periods$lon,calibration.periods$lat), 
                                                log.light.borders = "auto", log.irrad.borders = "auto", ylim = NULL) {
          old.par <- graphics::par(no.readonly = TRUE)
          Calibration.period <- data.frame(calibration.start = as.POSIXct("1900-01-01"), 
                                           calibration.stop = as.POSIXct("2050-01-01"), lon = location[1], 
                                           lat = location[2])
          if (log.light.borders[1] == "auto") 
            log.light.borders <- Proc.data$log.light.borders
          if (log.irrad.borders[1] == "auto") 
            log.irrad.borders <- Proc.data$log.irrad.borders
          graphics::par(old.par)
          calibration.parameters <- suppressWarnings(FLightR:::get.calibration.parameters(Calibration.period, 
                                                                                          Proc.data, model.ageing = FALSE, log.light.borders = log.light.borders, 
                                                                                          log.irrad.borders = log.irrad.borders, plot.each = FALSE, 
                                                                                          plot.final = FALSE, suggest.irrad.borders = FALSE))
          suppressWarnings(FLightR:::plot_slopes(calibration.parameters$All.slopes, 
                                                 ylim = ylim))
          
          return(calibration.parameters)
        }
        # calulate and plot
        Calibration.parameters <- plot_slopes_by_location_adjusted(Proc.data,location = muCoords)
        abline(v=ll$arrival,col="blue")
        abline(v=ll$departure,col="blue")
        dev.copy2pdf(file=file.path(stopoverDir,paste0(so,"_calibration_Entire.pdf")))
        dev.off()
        #seperate dawns and dusk slopes
        slopeData <- Calibration.parameters$All.slopes$Slopes
        dawnslopeData <- slopeData[slopeData$Type == "Dawn",]
        duskslopeData <- slopeData[slopeData$Type == "Dusk",]
        # set limits on the x and y values of plotting
        xlim <- c(as.POSIXct(as.Date(ll$arrival)) - 5*86400,
                  as.POSIXct(as.Date(ll$departure)) + 5*86400)
        slopeSub <- apply(cbind(!is.na(slopeData$logSlope),is.finite(slopeData$logSlope),
                                slopeData$Time >= xlim[1],slopeData$Time <= xlim[2]),1,all)
        ylim <- c(min(slopeData$logSlope[slopeSub]),max(slopeData$logSlope[slopeSub]))
        
        # Plot dawns
        pdf(file=file.path(stopoverDir,paste0(so,"_calibration.pdf")))
        plot(x=as.POSIXct(dawnslopeData$Time,origin = "1970-01-01"),
             y=dawnslopeData$logSlope,
             xlim=xlim,ylim=ylim,pch="+",col="red",
             xlab="Date",ylab="log(slope)")
        if (so %in% match(stopoversNC21,stopoversNC14)) {
          title(paste(birdFold,"location errors","-over 21 days"))
        } else {
          title(paste(birdFold,"location errors","-over 14 days"))
        }
        lines(x=as.POSIXct(dawnslopeData$Time,origin = "1970-01-01"),
              y=dawnslopeData$logSlope,col="red")
        # Plot dusks
        points(x=as.POSIXct(duskslopeData$Time,origin = "1970-01-01"),
               y=duskslopeData$logSlope,pch="+",col="black")
        lines(x=as.POSIXct(duskslopeData$Time,origin = "1970-01-01"),
              y=duskslopeData$logSlope,col="black")
        # plot stationary period start and end
        abline(v=ll$arrival,col="blue")
        abline(v=ll$departure,col="blue")
        dev.off()
        
        rmSO <- rep(FALSE)
        return(list(coords=muCoords,calibration=Calibration.parameters,useStopover=rmSO))
        
      },error=function(err){
        runError <- c(runError,paste("error on getting new location",birdFold,"(birdNum=",birdNum,")",err))
        if (!is.null(dev.list())) {dev.off()}
        print(paste("MY_ERROR: ",birdFold,err))
      })
    }
    
    #apply function to get stopover coordinates for all stopovers
    allOut <- lapply(ll,getStopoverCoords,stopNums, Result)
    
    outputCoords <- lapply(allOut,"[[",1)
    calib.para <- lapply(allOut,"[[",2)
    rmStopovers <- lapply(allOut,"[[",3)
    
    # check for null values in estimated coordinates (outputCoords)
    if (length(which(sapply(outputCoords,is.null))) > 0) {
      outputCoords[[which(sapply(outputCoords,is.null))]] <- c(NA,NA)
    }
    
    # save coordinates into RData file
    save(outputCoords, calib.para, rmStopovers, stopovers, statPeriod, pStatPeriod,
         file=file.path(stopoverDir,paste0("LocationDetermine",passNum,".RData")))
    
    subStopovers <-  which(stopoversNC14 %in% stopoversNC21)
    
    lons <- sapply(outputCoords,'[[',1)
    lats <- sapply(outputCoords,'[[',2)
    
    lons <- lons[!is.na(lons)]
    lats <- lats[!is.na(lats)]
    
    xVal <- c(statPeriod$Meanlon,lons)
    yVal <- c(statPeriod$Meanlat,lats)
    
    
    ## Plot the mean locations of stopover sites
    pdf(file=file.path(stopoverDir,paste0("stopoverLocations",passNum,".pdf")))
    plot(statPeriod$Meanlon,statPeriod$Meanlat,
         xlim=c(min(xVal),max(xVal)),ylim=c(min(yVal),max(yVal)),
         col="red",pch=nvn, main=paste0(birdFold,"_stopoverLocations"))
    lines(statPeriod$Meanlon,statPeriod$Meanlat,col="blue")
    maps::map("world",add=TRUE)
    maps::map("state",add=TRUE)
    points(x=sapply(outputCoords,'[[',1),
           y=sapply(outputCoords,'[[',2),
           col="blue",pch=19)
    points(x=sapply(outputCoords,'[[',1)[subStopovers],
           y=sapply(outputCoords,'[[',2)[subStopovers],
           col="orange",pch=19)
    dev.off()
    
  },warning = function(warn) {
    runError <- c(runError,paste("warning on",birdFold,"(birdNum=",birdNum,")",warn))
    print(paste("MY_WARNING: ",birdFold,warn))
  },error = function(err) {
    runError <- c(runError,paste("error on",birdFold,"(birdNum=",birdNum,")",err))
    print(paste("MY_ERROR: ",birdFold,err))
    if (!is.null(dev.list())) {dev.off()}
  })
  # save workspace to RData file
  tryCatch({
    save(Proc.data,known.last,last.loc,stopovers,stopoversNC14,stopoversNC21,
         outputCoords,cap.long,cap.lat,quants,
         file=file.path(stopoverDir,paste0("Stopover",passNum,"_values",".RData")))
  }, warning = function(warn) {
    runError <- c(runError,paste("warning on",birdFold,"(birdNum=",birdNum,")",warn))
    print(paste("MY_WARNING: ",birdFold,warn))
  }, error = function(err) {
    runError <- c(runError,paste("error on",birdFold,"(birdNum=",birdNum,")",err))
    print(paste("MY_ERROR: ",birdFold,err))
  })
  ls1 <- ls(envir=.GlobalEnv)
  save(list=ls(),file=file.path(stopoverDir,paste0("Stopover",passNum,"analysis.RData")))
  return(runError)
}




# Single combined stopover code -----------------
# Review stopover locations# set threholds of distance (km) and bearing (degrees) that are used to determine distinct stopover locations
run.combstopover.review <- function(birdFold, passNum, stopovers, threshDistance=150, birdList=fileList[index.bird],
                                    mainWD=mainDir, speciesWD=mainDir, geolocatorDF=geo_df,
                                    species=NULL) {
  #require(lubridate) ?????
  tryCatch({
    runError <- c()
    #pull species from geolocator data frame
    if (is.null(species)) {
      species=unique(geo_df$Species)
      print("species pulled from geolocator data frame")
    }
    birdNum <- which(birdList %in% birdFold)
    moveDir <- file.path(speciesWD,birdFold,"movement_data")
    passFold <- paste0("Movement",passNum)
    passDir <- file.path(moveDir,passFold)
    # gets variables from previous analysis - From selected variables saved to RData file
    load(file.path(passDir,paste0(passFold,"_values.RData")))
    #names and creates folder for stopover if it doesn't existe yet
    stopoverDir <- file.path(passDir,paste0("stopover_",passNum,"analysis"))
    if(!dir.exists(stopoverDir)){dir.create(stopoverDir)}
    
    ## Setup a dataframe of stationary period information
    statPeriod <- statMigSum$Stationary.periods
    # get labels for the stopovers
    nvn <- c(letters,as.character(c(0:9)))
    nvn <- nvn[1:nrow(statPeriod)]
    quants <- Result$Results$Quantiles
    pStatPeriod <- statMigSum$Potential_stat_periods
    pStatPeriod <- cbind(stopoverID=nvn[1:nrow(pStatPeriod)],pStatPeriod)
    pMovePeriod <- statMigSum$Potential_movement_periods
    if(pStatPeriod$start[1] == 0) {
      pStatPeriod$startDate <- c(quants$time[1] - 86400,quants$time[pStatPeriod$start])
    } else if(pStatPeriod$start[1] > 0) {
      pStatPeriod$startDate <- quants$time[pStatPeriod$start]
    }
    pStatPeriod$endDate <- quants$time[pStatPeriod$end]
    pStatPeriod$MedstartDate <- c(quants$time[1] - 86400,statPeriod$Arrival.Q.50[2:length(pStatPeriod$start)])
    pStatPeriod$MedendDate <- statPeriod$Departure.Q.50[1:length(pStatPeriod$end)]
    
    # Add columns for checks for distance and bearing
    naColumn <- rep(NA,nrow(statPeriod))
    statPeriod <- cbind(location = nvn, newStopover = naColumn,
                        overDistance = naColumn, distance = statPeriod$Distance.Moved,
                        statPeriod)
    # Mark the columns for those stopovers that are over distance threshold
    # also mark the distinct stopovers
    statPeriod$overDistance <- statPeriod$Distance.Moved > threshDistance
    statPeriod$newStopover <- apply(cbind(statPeriod$overDistance,statPeriod$overBearing),1,all)
    
    ## Plot the stopover locations
    # get length of stopover periods
    Stoplength <- pStatPeriod$Duration/2
    lengthGT14 <- Stoplength > 14
    lengthGT21 <- Stoplength > 21
    # add logical columns marked true for stopover durations longer than 2 weeks and 3 weeks
    pStatPeriod$lengthGT14 <- lengthGT14
    pStatPeriod$lengthGT21 <- lengthGT21
    
    stopoversNC14 <- which(lengthGT14)
    stopoversNC21 <- which(lengthGT21)
    
    # Plot the timeline with median start and end times for stopover locations
    lylim <- min(statPeriod$Distance.Moved)
    uylim <- max(statPeriod$Distance.Moved)
    #upper plot height (uph) gives the y value for plotting stopover IDs
    uph <- 500*floor(uylim/500)
    #create PDF of plot
    pdf(file=file.path(stopoverDir,paste0("combstopoverTiming",passNum,".pdf")))
    plot(quants$time,rep(0,length(quants$time)),
         ylim=c(lylim,uylim), type="n", main = paste0(birdFold,"_stopoverTimeline"))
    if (!is.null(nestingDates)) {
      abline(v=as.POSIXct(nestingDates,tm="GMT"),col="purple",lwd=4)
    }
    rect(xleft = pStatPeriod$startDate, xright = pStatPeriod$endDate,
         ybottom = lylim - 1,  ytop = uylim + 1, col="green")
    #plot the combined stopovers
    rect(xleft = stopovers$arrival, xright = stopovers$departure,
         ybottom = lylim - 1,  ytop = uylim/2, col="black")
    
    # plot bars above the stopover periods that are blue for over 14 days and orange for over 21 days
    if(length(stopoversNC14) > 0) {
      rect(xleft = pStatPeriod$startDate[lengthGT14], xright = pStatPeriod$endDate[lengthGT14],
           ybottom = uylim + 1,  ytop = uylim + 30, col="blue")
    }
    if(length(stopoversNC21) > 0) {
      rect(xleft = pStatPeriod$startDate[lengthGT21], xright = pStatPeriod$endDate[lengthGT21],
           ybottom = uylim + 1,  ytop = uylim + 30, col="orange")
    }
    
    # get the date/time that is centered between mean start and end dates of stopovers
    centerDates <- pStatPeriod$startDate + (pStatPeriod$endDate - pStatPeriod$startDate)/2
    # Plot the stopover location names on the stopover locations rectangles
    points(centerDates,rep(uph,length(centerDates)),pch=nvn[1:(length(nvn)-1)])
    # Add the put on and take off dates
    abline(v=as.POSIXct(geolocatorDF$Date.of.Initial.fitting[birdNum],tz="GMT"),col="orange",lwd=2)
    # Plot the distances moved from previous stopover location to current location
    # mark threshold for considering distinct stopover locations
    points(centerDates,statPeriod$Distance.Moved[1:length(centerDates)],col="red",pch=19)
    abline(h=150,col="blue")
    dev.off()
    
    ## Get the estimated locations from stopover periods
    initStops <- stopoversNC14
    soNum <- nrow(pStatPeriod)
    stopNums <- stopovers$stopNum
    #create a list of stopovers data
    ll <- split(stopovers,seq(nrow(stopovers)))
    
    getStopoverCoords <- function(ll, stopNums, Result,gelmanCheck) {
      tryCatch({
        so <- which(stopNums %in% ll$stopNum)
        ## get stopover coordinates from FLightR particle filter results
        {spG <- Result$Spatial$Grid
          twlSeq <- which(Result$Results$Quantiles$time >= ll$arrival & 
                            Result$Results$Quantiles$time <= ll$departure)
          prtPost <- Result$Results$Points.rle[twlSeq]
          # get all coordinates
          prtLoc <- spG[unlist(sapply(prtPost,"[[",2)),]
          # get lengths
          prtLen <- unlist(sapply(prtPost,"[[",1))
          #mean lon
          mulon <- sum(prtLoc[,1] * prtLen)/(nParticles*length(twlSeq))
          #mean lat
          mulat <- sum(prtLoc[,2] * prtLen)/(nParticles*length(twlSeq))
          #mean coordinates
          muCoords <- c(mulon,mulat)
        }
        
        # look at log(slope) error plots for each calibration period
        #Function taken from FLightR named plot_slopes_by_location and returns parameter data as well as plotting
        plot_slopes_by_location_adjusted <- function (Proc.data, location = c(calibration.periods$lon,calibration.periods$lat), 
                                                      log.light.borders = "auto", log.irrad.borders = "auto", ylim = NULL) {
          old.par <- graphics::par(no.readonly = TRUE)
          Calibration.period <- data.frame(calibration.start = as.POSIXct("1900-01-01"), 
                                           calibration.stop = as.POSIXct("2050-01-01"), lon = location[1], 
                                           lat = location[2])
          if (log.light.borders[1] == "auto") 
            log.light.borders <- Proc.data$log.light.borders
          if (log.irrad.borders[1] == "auto") 
            log.irrad.borders <- Proc.data$log.irrad.borders
          graphics::par(old.par)
          calibration.parameters <- suppressWarnings(FLightR:::get.calibration.parameters(Calibration.period, 
                                                                                          Proc.data, model.ageing = FALSE, log.light.borders = log.light.borders, 
                                                                                          log.irrad.borders = log.irrad.borders, plot.each = FALSE, 
                                                                                          plot.final = FALSE, suggest.irrad.borders = FALSE))
          suppressWarnings(FLightR:::plot_slopes(calibration.parameters$All.slopes, 
                                                 ylim = ylim))
          
          return(calibration.parameters)
        }
        
        # calulate and plot
        Calibration.parameters <- plot_slopes_by_location_adjusted(Proc.data,location = muCoords)
        abline(v=ll$arrival,col="blue")
        abline(v=ll$departure,col="blue")
        #save plot as PDF
        dev.copy2pdf(file=file.path(stopoverDir,paste0("comb",so,"_calibration_Entire.pdf")))
        dev.off()
        
        #seperate dawns and dusk slopes
        slopeData <- Calibration.parameters$All.slopes$Slopes
        dawnslopeData <- slopeData[slopeData$Type == "Dawn",]
        duskslopeData <- slopeData[slopeData$Type == "Dusk",]
        # set limits on the x and y values of plotting
        xlim <- c(as.POSIXct(as.Date(ll$arrival)) - 5*86400,
                  as.POSIXct(as.Date(ll$departure)) + 5*86400)
        slopeSub <- apply(cbind(!is.na(slopeData$logSlope),is.finite(slopeData$logSlope),
                                slopeData$Time >= xlim[1],slopeData$Time <= xlim[2]),1,all)
        ylim <- c(min(slopeData$logSlope[slopeSub]),max(slopeData$logSlope[slopeSub]))
        
        #plot and save dawns as PDF
        pdf(file=file.path(stopoverDir,paste0("comb",so,"_calibration.pdf")))
        plot(x=as.POSIXct(dawnslopeData$Time,origin = "1970-01-01"),
             y=dawnslopeData$logSlope,
             xlim=xlim,ylim=ylim,pch="+",col="red",
             xlab="Date",ylab="log(slope)")
        title(paste(birdFold,"location errors"))
        lines(x=as.POSIXct(dawnslopeData$Time,origin = "1970-01-01"),
              y=dawnslopeData$logSlope,col="red")
        # Plot and save dusks as PDF
        points(x=as.POSIXct(duskslopeData$Time,origin = "1970-01-01"),
               y=duskslopeData$logSlope,pch="+",col="black")
        lines(x=as.POSIXct(duskslopeData$Time,origin = "1970-01-01"),
              y=duskslopeData$logSlope,col="black")
        # plot stationary period start and end
        abline(v=ll$arrival,col="blue")
        abline(v=ll$departure,col="blue")
        dev.off()
        
        rmSO <- rep(FALSE)
        return(list(coords=muCoords,calibration=Calibration.parameters,useStopover=rmSO))
        
      },error=function(err){
        runError <- c(runError,paste("error on getting new location",birdFold,"(birdNum=",birdNum,")",err))
        if (!is.null(dev.list())) {dev.off()}
        print(paste("MY_ERROR: ",birdFold,err))
      })
    }
    
    #apply function to get stopover coordinates for all stopovers
    allOut <- lapply(ll,getStopoverCoords,stopNums, Result)
    
    outputCoords <- lapply(allOut,"[[",1)
    calib.para <- lapply(allOut,"[[",2)
    rmStopovers <- lapply(allOut,"[[",3)
    
    # check for null values in estimated coordinates (outputCoords)
    if (length(which(sapply(outputCoords,is.null))) > 0) {
      outputCoords[[which(sapply(outputCoords,is.null))]] <- c(NA,NA)
    }
    
    subStopovers <-  which(stopoversNC14 %in% stopoversNC21)
    
    lons <- sapply(outputCoords,'[[',1)
    lats <- sapply(outputCoords,'[[',2)
    
    stopovers$meanlat <- lats
    stopovers$meanlon <- lons
    
    # save coordinates into RData file
    save(outputCoords, calib.para, rmStopovers, stopovers, statPeriod, pStatPeriod,
         file=file.path(stopoverDir,paste0("CombLocationDetermine",passNum,".RData")))
    
    lons <- lons[!is.na(lons)]
    lats <- lats[!is.na(lats)]
    
    xVal <- c(statPeriod$Meanlon,lons)
    yVal <- c(statPeriod$Meanlat,lats)
    
    ## Plot the mean locations of stopover sites
    pdf(file=file.path(stopoverDir,paste0("combStopoverLocations",passNum,".pdf")))
    plot(statPeriod$Meanlon,statPeriod$Meanlat,
         xlim=c(min(xVal),max(xVal)),ylim=c(min(yVal),max(yVal)),
         col="red",pch=nvn, main=paste0(birdFold,"_stopoverLocations"))
    lines(statPeriod$Meanlon,statPeriod$Meanlat,col="blue")
    maps::map("world",add=TRUE)
    maps::map("state",add=TRUE)
    points(x=sapply(outputCoords,'[[',1),
           y=sapply(outputCoords,'[[',2),
           col="orange",pch=19)
    text(x=sapply(outputCoords,'[[',1)+.25,
         y=sapply(outputCoords,'[[',2)+.25,
         labels = stopovers$stopNum)
    dev.off()
    
  },error = function(err) {
    runError <- c(runError,paste("error on",birdFold,"(birdNum=",birdNum,")",err))
    print(paste("MY_ERROR: ",birdFold,err))
    if (!is.null(dev.list())) {dev.off()}
  })
  # save workspace to RData file
  tryCatch({
    save(Proc.data,known.last,last.loc,stopovers,stopoversNC14,stopoversNC21,
         outputCoords,cap.long,cap.lat,quants,
         file=file.path(stopoverDir,paste0("CombStopover",passNum,"_values",".RData")))
  }, error = function(err) {
    runError <- c(runError,paste("error on",birdFold,"(birdNum=",birdNum,")",err))
    print(paste("MY_ERROR: ",birdFold,err))
  })
  ls1 <- ls(envir=.GlobalEnv)
  save(list=ls(),file=file.path(stopoverDir,paste0("CombStopover",passNum,"analysis.RData")))
  return(list(coords=outputCoords,rErr=runError))
}




# get Excluded stopover periods --------------------
getExcludedStopovers <- function(birdFold,passNum,birdList=fileList[index.bird],
                                 speciesWD=mainDir,geolocatorDF=geo_df,
                                 bs=breedStart,be=breedEnd,ms=moltStart,me=moltEnd) {
  try({
    require(lubridate)
    #create file path to specific folders
    moveDir <- file.path(speciesWD,birdFold,"movement_data")
    passFold <- paste0("Movement",passNum)
    passDir <- file.path(moveDir,passFold)
    stopoverFold <- paste0("stopover_",passNum,"analysis")
    stopoverDir <- file.path(passDir,stopoverFold)
    #load previous data files
    load(file.path(stopoverDir,paste0("LocationDetermine",passNum,".RData")))
    
    birdNum <- which(birdList %in% birdFold)
    #stopovers excluded from previous visual inspection of data
    soex1 <- unlist(rmStopovers)
    soLen <- length(soex1)
    #create list of stopovers
    soTimeList <- lapply(split(cbind(stopovers,years=year(stopovers$arrival)),seq(nrow(stopovers))),
                         function(x) list(times=as.Date(seq(x$arrival,x$departure,1*86400)),years=x$years))
    names(soTimeList) <- paste0("so",1:soLen)
    #stopovers that occured while geolocator is on bird
    soex2 <- sapply(soTimeList,function(x) all(x$times < geolocatorDF$Date.of.Initial.fitting[birdNum]))
    #stopovers occuring between start and end dates of breeding
    soex3 <- lapply(soTimeList,function(x) x$times >= as.Date(paste0(x$years,bs)) &
                      x$times <= as.Date(paste0(x$years,be)))
    names(soex3) <- paste0("so",1:soLen)
    soex3in <- sapply(soex3,function(x) length(which(x)) < .25*length(x))
    soex3 <- sapply(soex3,any)
    
    #stopovers occuring between start and end dates of molting - with IF statement for those birds that don't have lux data during that time
    if(!is.null(ms) | !is.null(me)) {
      soex4 <- lapply(soTimeList,function(x) x$times >= as.Date(paste0(x$years,ms)) &
                        x$times <= as.Date(paste0(x$years,me)))
      names(soex4) <- paste0("so",1:soLen)
      soex4in <- sapply(soex4,function(x) length(which(x)) < .25*length(x))
      soex4 <- sapply(soex4,any)
    } else {
      message("No molt period provided, so no stopovers from molt period used")
      soex4 <- !soex3
      names(soex4) <- paste0("so",1:soLen)
      soex4in <- soex4
    }
    
    # Checking for any conditions that the stopover should be excluded
    # (the conditions explained in order)
    # 1-marked remove from stopover analysis,
    # 2-stopovers before geo attached to bird,
    # 3- stopovers with any part inside of breeding or molting ,
    # 4- stopovers with less than 25% of duration inside breeding,
    # 5- stopovers with less than 25% of duration inside breeding,
    soex <- as.numeric(which(soex1 | soex2 | !soex3 & !soex4 | soex3 & soex3in | soex4 & soex4in))
    return(soex)
  })
  
}




# get Combined stopover periods --------------------
getCombinedStopovers <- function(birdFold,passNum,includedThreshold=NULL,birdList=fileList[index.bird],
                                 speciesWD=mainDir,geolocatorDF=geo_df,#stopovers,
                                 bs=breedStart,be=breedEnd,ms=moltStart,me=moltEnd) {
  require(lubridate)
  moveDir <- file.path(speciesWD,birdFold,"movement_data")
  passFold <- paste0("Movement",passNum)
  stopoverFold <- paste0("stopover_",passNum,"analysis")
  load(file.path(moveDir,passFold,stopoverFold,paste0("LocationDetermine",passNum,".RData")))
  
  birdNum <- which(birdList %in% birdFold)
  
  soLen <- nrow(pStatPeriod)
  nvn <- c(letters,as.character(c(0:9)))
  stopovers$stopNum <- nvn[stopovers$stopNum]
  soTimeList <- lapply(split(cbind(pStatPeriod,years=year(pStatPeriod$startDate)),seq(nrow(pStatPeriod))),
                       function(x) list(times=seq(x$startDate,x$endDate,1*86400),years=x$years))
  names(soTimeList) <- paste0("so",1:soLen)
  # logical vector for stopover after inital fitting date
  soex2 <- sapply(soTimeList,function(x) all(x$times < geolocatorDF$Date.of.Initial.fitting[birdNum]))
  # logical vector for stopover between breeding period dates
  soex3 <- lapply(soTimeList,function(x) x$times >= as.Date(paste0(x$years,bs)) &
                    x$times <= as.Date(paste0(x$years,be)))
  names(soex3) <- paste0("so",1:soLen)
  # logical vector for stopover between molt migration dates
  if(!is.null(ms) | !is.null(me)) {
    soex4 <- lapply(soTimeList,function(x) x$times >= as.Date(paste0(x$years,ms)) &
                      x$times <= as.Date(paste0(x$years,me)))
    names(soex4) <- paste0("so",1:soLen)
  } else {
    message("No molt period provided, so no stopovers from molt period used")
    soex4 <- list()
    for (ex in 1:soLen) {
      soex4[[ex]] <- FALSE
    }
    names(soex4) <- paste0("so",1:soLen)
  }
  
  
  if(is.null(includedThreshold)) {
    soex3in <- sapply(soex3,any)
    soex4in <- sapply(soex4,any)
  } else {
    soex3in <- sapply(soex3,function(x) length(which(x)) > .5*length(x))
    soex4in <- sapply(soex4,function(x) length(which(x)) > .5*length(x))
  }
  #create empty vectors
  stopNum <- c()
  meanlon <- c()
  meanlat <- c()
  arrival <- c()
  departure <- c()
  
  #combine soex3 - stopovers within breeding
  if(any(soex3in)) {
    subSoTimeList <- lapply(soTimeList[soex3in],function(x) x$times)
    soTimeVec <- unlist(subSoTimeList)#as.vector(unlist(subSoTimeList))
    soTimeYearList <- split(soTimeVec,year(as.POSIXct(soTimeVec,origin = "1970-01-01")))
    soNameYearList <- split(as.vector(pStatPeriod$stopoverID[which(soex3in)]),sapply(soTimeList[soex3in],function(x) unique(year(x$times))))
    stopNum <- c(stopNum,sapply(soNameYearList,function(x) paste(x,collapse=",")))
    arrival <- c(arrival,sapply(soTimeYearList,function(x) x[1])) #c(arrival,soTimeVec[1])
    departure<- c(departure,sapply(soTimeYearList,function(x) x[length(x)])) #c(departure,soTimeVec[length(soTimeVec)])
  }
  
  #combine soex4 - stopovers within molting
  if(any(soex4in)) {
    subSoTimeList <- lapply(soTimeList[soex4in],function(x) x$times)
    soTimeVec <- unlist(subSoTimeList)#as.vector(unlist(subSoTimeList))
    soTimeYearList <- split(soTimeVec,year(as.POSIXct(soTimeVec,origin = "1970-01-01")))
    soNameYearList <- split(pStatPeriod$stopoverID[which(soex4in)],sapply(soTimeList[soex4in],function(x) unique(year(x$times))))
    stopNum <- c(stopNum,sapply(soNameYearList,function(x) paste(x,collapse=",")))
    arrival <- c(arrival,sapply(soTimeYearList,function(x) x[1])) #c(arrival,soTimeVec[1])
    departure<- c(departure,sapply(soTimeYearList,function(x) x[length(x)])) #c(departure,soTimeVec[length(soTimeVec)])
  }
  
  # return(soex)
  return(data.frame(stopNum=as.vector(stopNum),
                    meanlon=rep(NA,length(arrival)),#meanlon,
                    meanlat=rep(NA,length(arrival)),#meanlat,
                    arrival=as.POSIXct(arrival,origin="1970-01-01"),
                    departure=as.POSIXct(departure,origin="1970-01-01")))
}



# set Combined stopover periods --------------------
setCombinedStopovers <- function(birdFold,passNum,soCombination=reviewSO,includedThreshold=NULL,
                                 birdList=fileList[index.bird],speciesWD=mainDir,
                                 geolocatorDF=geo_df,bs=breedStart,be=breedEnd,ms=moltStart,me=moltEnd) {
  require(lubridate)
  moveDir <- file.path(speciesWD,birdFold,"movement_data")
  passFold <- paste0("Movement",passNum)
  stopoverFold <- paste0("stopover_",passNum,"analysis")
  load(file.path(moveDir,passFold,stopoverFold,paste0("LocationDetermine",passNum,".RData")))
  
  #create empty vectors
  stopNum <- c()
  meanlon <- c()
  meanlat <- c()
  arrival <- c()
  departure <- c()
  
  stopNum <- c(stopNum,soCombination[1])
  rows <- match(unlist(strsplit(soCombination[1],",")),pStatPeriod$stopoverID)
  #include mean coordinates, arrival date and departure date for combined stopover
  meanlon <- c(meanlon,NA)
  meanlat <- c(meanlat,NA)
  arrival <- c(arrival,pStatPeriod$startDate[rows[1]])
  departure<- c(departure,pStatPeriod$endDate[rows[length(rows)]])
  
  if(length(soCombination) > 1) {
    stopNum <- c(stopNum,soCombination[2])
    rows <- match(unlist(strsplit(soCombination[2],",")),pStatPeriod$stopoverID)
    
    meanlon <- c(meanlon,NA)
    meanlat <- c(meanlat,NA)
    arrival <- c(arrival,pStatPeriod$startDate[rows[1]])
    departure<- c(departure,pStatPeriod$endDate[rows[length(rows)]])
  }
  
  return(data.frame(stopNum=as.vector(stopNum),
                    meanlon=meanlon,
                    meanlat=meanlat,
                    arrival=as.POSIXct(arrival,origin = "1970-01-01"),
                    departure=as.POSIXct(departure,origin = "1970-01-01")))
}



# Review the results -------------------
# A function to open stopover map and timing figures for each bird - to ease visual verification of automatic stopover assumptions
openGraphics <- function (birdFold,passNum,showCal,showComb=FALSE,
                          showMoveRes=FALSE,showSO=FALSE,speciesWD=mainDir) {
  # Go through all files to get the review files
  moveDir <- file.path(speciesWD,birdFold,"movement_data")
  passFold <- paste0("Movement",passNum)
  stopFold <- paste0("stopover_",passNum,"analysis")
  
  # before 3rd analysis
  if (showCal) {
    soFileList <- list.files(file.path(moveDir,passFold,stopFold))
    soFileList <- soFileList[grep("calibration",soFileList)]
    soFileList <- soFileList[!grepl("Entire",soFileList)]
    combSo <- grepl("comb",soFileList)
    if (showComb) {
      soFileList <- soFileList[combSo]#get comb ones only
      sapply(file.path(moveDir,passFold,stopFold,soFileList),shell.exec)
    } else {
      soFileList <- soFileList[!combSo]#get uncombined stopovers only
      sapply(file.path(moveDir,passFold,stopFold,soFileList),shell.exec)
    }
  }
  
  
  # After 3rd Analysis
  # open files
  if (showMoveRes) {
    shell.exec(file.path(moveDir,passFold,paste0(birdFold,".",passNum,"map.pdf")))
    shell.exec(file.path(moveDir,passFold,paste0("lonlat",passNum,".pdf")))
  }
  
  if (showSO) {
    shell.exec(file.path(moveDir,passFold,stopFold,paste0("stopoverLocations",passNum,".pdf")))
    shell.exec(file.path(moveDir,passFold,stopFold,paste0("stopoverTiming",passNum,".pdf")))
  }
  Sys.sleep(2)
  wdr <- winDialog("ok","Done with this bird?")
  
  return(wdr)
}



#first FLightR run----------------
# Code to run through functions for first pass of analysis
#with nParticles=1e6 and running through all 8 birds it takes many hours - can run at 1e4 for faster, less precise run
birdList <- fileList[index.bird]
nParticles <- 1e6
threads <- 7
passNum <- 1
system.time(
  runError <- sapply(birdList,run.FLightR.movement,passNum=passNum, nParticles=nParticles,
                     runCalib=TRUE,runPF=TRUE,runMap=TRUE,Threads=threads)
)
ulError <- unlist(runError)
if(!is.null(ulError)) {
  write(ulError,file.path(mainDir,paste0(species,"_FLightR_runErrors_pass",passNum,"_",Sys.Date(),".txt")))
}
rm(passNum)
#first stopover review
passNum <- 1
system.time(
  runError <- sapply(birdList,run.stopover.review,passNum=passNum)
  # runerror <- run.stopover.review(birdList,passNum=passNum)
)
ulError <- unlist(runError)
if(!is.null(ulError)) {
  write(ulError,file.path(mainDir,paste0(species,"_stopover_runErrors_pass",passNum,"_",Sys.Date(),".txt")))
}
#first get excluded stopover periods
stopPeriodsExcluded <- lapply(birdList,getExcludedStopovers,passNum=passNum)
names(stopPeriodsExcluded) <- birdList
stopPeriodsExcluded
sapply(birdList, openGraphics, passNum=1,showCal=FALSE,showComb=FALSE,showMoveRes=FALSE,showSO=TRUE)
save(stopPeriodsExcluded,file=file.path(mainDir,paste0("Excluded_stopovers",passNum,".RData")))
rm(passNum)


#second FLightR run------------------
# Code to run through functions for second pass of analysis
#with nParticles=1e6 and running through all 8 birds it takes many hours - can run at 1e4 for faster, less precise run
nParticles <- 1e6
threads <- 7
passNum <- 2
if(!exists("stopPeriodsExcluded")) {
  load(file.path(mainDir,"Excluded_stopovers1.RData"))
  cat("stopPeriodsExcluded loaded from file")
}
system.time(
  runError <- sapply(birdList,run.FLightR.movement,passNum=passNum,nParticles=nParticles,
                     stopPeriodsExcluded=stopPeriodsExcluded,
                     runCalib=TRUE,runPF=TRUE,runMap=TRUE,Threads=threads)
)
ulError <- unlist(runError)
if(!is.null(ulError)) {
  write(ulError,file.path(mainDir,paste0(species,"_FLightR_runErrors_pass",passNum,"_",Sys.Date(),".txt")))
}
rm(passNum)
#second stopover review
passNum <- 2
system.time(
  runError <- sapply(birdList,run.stopover.review,passNum=passNum)
  # runError <- run.stopover.review(birdList,passNum=passNum)
)
ulError <- unlist(runError)
if(!is.null(ulError)) {
  write(ulError,file.path(mainDir,paste0(species,"_stopover_runErrors_pass",passNum,"_",Sys.Date(),".txt")))
}
# combine stopover periods - change function
stopoversComb <- lapply(birdList,getCombinedStopovers,passNum=passNum)
# stopoversComb <- getCombinedStopovers(birdList,passNum = passNum)
names(stopoversComb) <- birdList
sapply(birdList, openGraphics, passNum=2,showCal=FALSE,showComb=FALSE,showMoveRes=FALSE,showSO=TRUE)
## After a check if the combined stopovers need to be changed do so with the following 3 lines of code
# reviewSO <- lapply(stopoversComb, function(x) as.vector(x$stopNum))
# reviewSO$CANV_208739815[1] <- "e,f"
# stopoversComb <- mapply(a=birdList,b=passNum,c=reviewSO,
#        FUN=function(a,b,c) lapply(a,FUN = setCombinedStopovers,passNum=b,soCombination=c))

save(stopoversComb,file=file.path(mainDir,paste0("Combined_stopovers",passNum,".RData")))
# get coordinates
allOut <- mapply(a=birdList,b=passNum,c=stopoversComb,
                 FUN=function(a,b,c) lapply(a, run.combstopover.review,passNum=b,stopovers=c))
# allOut <- run.combstopover.review(birdList,passNum = passNum,stopovers = stopoversComb)

# #may not be necessary to run these lines
ulError <- unlist(lapply(allOut,function(x) x$rErr))
if(!is.null(ulError)) {
  write(ulError,file.path(mainDir,paste0(species,"_Combstopover_runErrors_pass",passNum,"_",Sys.Date(),".txt")))
}

sapply(birdList, openGraphics, passNum=2,showCal=TRUE,showComb=TRUE,showMoveRes=FALSE,showSO=FALSE)
# openGraphics(birdList,passNum = 2,showCal=TRUE,showComb=TRUE,showMoveRes=FALSE,showSO=TRUE)
rm(passNum)

#third FLightR run--------------
# Code to run through functions for third pass of analysis
#with nParticles=1e6 and running through all 8 birds it takes many hours - can run at 1e4 for faster, less precise run
nParticles <- 1e6
threads <- 7
passNum <- 3
system.time(
  runError <- sapply(birdList,run.FLightR.movement,passNum=passNum, nParticles=nParticles,
                     runCalib=TRUE,runPF=TRUE,runMap=TRUE,Threads=threads)
  # runError <- run.FLightR.movement(birdList,passNum = passNum,nParticles = nParticles)
)
ulError <- unlist(runError)
if(!is.null(ulError)) {
  write(ulError,file.path(mainDir,paste0(species,"_FLightR_runErrors_pass",passNum,"_",Sys.Date(),".txt")))
}
rm(passNum)
#third stopover review
passNum <- 3
system.time(
  runError <- sapply(birdList,run.stopover.review,passNum=passNum)
  # runError <- run.stopover.review(birdList,passNum=passNum)
)

ulError <- unlist(runError)
if(!is.null(ulError)) {
  write(ulError,file.path(mainDir,paste0(species,"_stopover_runErrors_pass",passNum,"_",Sys.Date(),".txt")))
}
# look over graphics outputs
sapply(birdList, openGraphics, passNum=3,showCal=TRUE,showComb=TRUE,showMoveRes=TRUE,showSO=TRUE)
# openGraphics(birdList,passNum = 2,showCal=TRUE,showComb=TRUE,showMoveRes=TRUE,showSO=TRUE)
rm(passNum)



# Mapping the stopover locations and timelines####