## FLightR analysis using functions that use multiple passes
## updated 5/8/20
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

## Source functions from script------------
source("Functions.R")



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