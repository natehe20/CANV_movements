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


## map color and stopover timing setup -------
# plot for presentations
## set colors for different stopovers
#spring stopovers
sColor <- "darkgreen" #dark green
#breeding area color
bColor <- "green3"
#molt/staging area color
mColor <- rgb(238,154,0,alpha=150, max=255)#tan
mColor <- "gold2"
#fall stopovers
fColor <- rgb(238,154,0,alpha=250, max=255)
fColor <- 'goldenrod3'
#wintering area color
wColor <- rgb(108,166,205,alpha=150, max=255)#light blue
wColor <- "skyblue"

##Visually inspect the stopover timeline, stopover locations and pStatPeriod from stopover analysis in order to 
#assign stopovers to Spring, breeding, molting, fall or winter timeframe

#Set the stopover types: Spring, breeding, molting, fall and wintering
sso <- list(CANV_208739806=c(1,2),
            CANV_208739809=c(1:4),
            CANV_208739815=c(2:5),
            CANV_208739817=c(1:2),
            CANV_210701439=c(2:3),
            CANV_210701300=c(2:6),
            CANV_210701364=c(5:6),
            CANV_210701385=c(2:4),
            CANV_210701300=c(13:19),
            CANV_210701364=c(14:15)
            ##CANV_210701385=c()
)
bso <- list(CANV_208739806=c(3),
            CANV_208739809=c(5:6),
            CANV_208739815=c(6:7),
            CANV_208739817=c(3:6),
            CANV_210701439=c(4),
            CANV_210701300=c(7:8),
            CANV_210701364=c(7),
            CANV_210701385=c(5),
            CANV_210701300=c(20),
            CANV_210701364=c(16)
            ##CANV_210701385=c()
)
mso <- list(CANV_208739806=c(4),
            CANV_208739809=c(8:9),
            CANV_208739815=c(8:9),
            CANV_208739817=c(7),
            CANV_210701439=c(5),
            CANV_210701300=c(9),
            CANV_210701364=c(8),
            CANV_210701385=c(7),
            CANV_210701300=c(21),
            CANV_210701364=c(16)
            ##CANV_210701385=c()
)
fso <- list(CANV_208739806=c(5),
            CANV_208739809=c(7,10:11),
            CANV_208739815=c(10:12),
            CANV_208739817=c(),
            CANV_210701439=c(6:9),
            CANV_210701300=c(10:11),
            CANV_210701364=c(9:12),
            CANV_210701385=c(6,8:9),
            CANV_210701300=c(22:24),
            CANV_210701364=c(17:18)
            ##CANV_210701385=c()
)
wso <- list(CANV_208739806=c(6),
            CANV_208739809=c(12:13),
            CANV_208739815=c(),#
            CANV_208739817=c(),#
            CANV_210701439=c(10),
            CANV_210701300=c(12),
            CANV_210701364=c(13),
            CANV_210701385=c(10),
            CANV_210701300=c(),#
            CANV_210701364=c(19:20)
            ##CANV_210701385=c()
)



##Maps
birdList <- fileList[index.bird]

# Map stopovers and errors -----------------
get.maps <- function(birdFold, passNum=3, sCol=sColor,bCol=bColor,mCol=mColor,fCol=fColor,wCol=wColor,
                     ol=FALSE,oso=NULL,ouso=NULL,omupnt=NULL,oerr=NULL,opath=NULL,oaerr=NULL,omloc=NULL,
                     springSO=sso,breedSO=bso,moltSO=mso,fallSO=fso,winterSO=wso,
                     projName="merc",projectionName=CRS("+init=epsg:3857"),birdList=fileList[index.bird],
                     speciesWD=mainDir, geolocatorDF=geo_df,
                     species=NULL,calp=cal.period, lcorr=l.correct) {
  require(shape)
  require(maps)
  require(sp)
  require(rgeos)
  require(grDevices)
  require(rgdal)
  require(raster)
  
  tryCatch({
    runError <- c()
    env <- environment()
    #sets overwrite authorization for shapefiles
    #sets null overwrite values to default
    if (is.null(oso)) {oso <- ol}
    if (is.null(ouso)) {ouso <- ol}
    if (is.null(omupnt)) {omupnt <- ol}
    if (is.null(oerr)) {oerr <- ol}
    if (is.null(opath)) {opath <- ol}
    if (is.null(oaerr)) {oaerr <- ol}
    if (is.null(omloc)) {omloc <- ol}
    
    if (is.null(species)) {
      species=unique(geo_df$Species)
      print("species pulled from geolocator data frame")
    }
    birdNum <- which(birdList %in% birdFold)
    
    m.folder <- paste(calp,"_",lcorr,sep="")
    passFold <- paste0("Movement",passNum)
    passDir <- file.path(speciesWD,birdFold,m.folder,passFold)
    stopoverFold <- paste0("stopover_",passNum,"analysis")
    # sets name and creates folder for GIS files if it doesn't exist
    gisDir <- file.path(speciesWD,birdFold,"GIS")
    if (!dir.exists(file.path(speciesWD,birdFold,"GIS"))) {
      dir.create(file.path(speciesWD,birdFold,"GIS"))
    }
    #load location values from Rdata output files
    load(file.path(passDir,stopoverFold,paste0("Stopover",passNum,"analysis.RData")))
    
    quants <- Result$Results$Quantiles
    
    #gets 1st quantile, mean and 3rd quantile of location values for each twilight
    l <- "FstQu."
    m <- "Mean"
    u <- "TrdQu."
    llon <- paste0(l,"lon")
    mlon <- paste0(m,"lon")
    ulon <- paste0(u,"lon")
    llat <- paste0(l,"lat")
    mlat <- paste0(m,"lat")
    ulat <- paste0(u,"lat")
    
    # get subset of twilights to be plotted - use twilights after geo deployment
    twlStart <- which(as.Date(quants$time) > as.Date(geo_df$Date.of.Initial.fitting[birdNum]))[1]
    subby <- twlStart:nrow(quants)
    
    #gets extent of lat and long for plotting
    lats <- range(c(quants[subby,ulat],quants[subby,llat]))
    lons <- range(c(quants[subby,ulon],quants[subby,llon]))
    
    #starts plot of locations
    plot(quants[subby,mlon],quants[subby,mlat],
         xlim=c(lons),
         ylim=c(lats),
         col="red",pch=nvn, main=paste0(birdFold,"_stopoverLocations"),type="n")
    ##    lines(statPeriod$Meanlon,statPeriod$Meanlat,col="blue")
    map("world",add=TRUE)
    map("state",add=TRUE)
    
    # get vertices coordinates for ellipse as estimated error in locations for each twilight
    get.error.ellipse <- function (subby,quants) {
      elliList <- mapply(a=list(abs(quants[subby,ulon]-quants[subby,mlon]),#radius TR x (lon)
                                abs(quants[subby,llon]-quants[subby,mlon]),#radius TL x (lon)
                                abs(quants[subby,llon]-quants[subby,mlon]),#radius BL x (lon)
                                abs(quants[subby,ulon]-quants[subby,mlon])),#radius BR x (lon)
                         b=list(abs(quants[subby,ulat]-quants[subby,mlat]),#radius TR y (lat)
                                abs(quants[subby,ulat]-quants[subby,mlat]),#radius TL y (lat)
                                abs(quants[subby,llat]-quants[subby,mlat]),#radius BL y (lat)
                                abs(quants[subby,llat]-quants[subby,mlat])), #radius BR y (lat)
                         c=list(matrix(c(quants[subby,mlon],quants[subby,mlat]),ncol=2)),
                         # c=list(c(quants[subby,mlon],quants[subby,mlat])),
                         d=0,
                         e=2*pi,
                         FUN=function(a,b,c,d,e) lapply(a,FUN=getellipse,ry=b,mid=c,from=d,to=e))
      
      elliList[[1]] <- elliList[[1]][-nrow(elliList[[1]]),]
      elliList[[4]] <- elliList[[4]][-1,]
      elliList[[1]] <- elliList[[1]][elliList[[1]][,1] >= quants[subby,mlon] & 
                                       elliList[[1]][,2] >= quants[subby,mlat],]
      elliList[[2]] <- elliList[[2]][elliList[[2]][,1] <= quants[subby,mlon] & 
                                       elliList[[2]][,2] >= quants[subby,mlat],]
      elliList[[3]] <- elliList[[3]][elliList[[3]][,1] <= quants[subby,mlon] & 
                                       elliList[[3]][,2] <= quants[subby,mlat],]
      elliList[[4]] <- elliList[[4]][elliList[[4]][,1] >= quants[subby,mlon] & 
                                       elliList[[4]][,2] <= quants[subby,mlat],]
      
      #get all vertices locations into one matrix for plotting
      uEL <- unlist(elliList)
      xel <- uEL[uEL < 0]
      yel <- uEL[uEL > 0]
      el <- cbind(xel,yel)
      return(el)
    }
    EL <- lapply(subby,get.error.ellipse,quants)
    
    #convert the vertices coordinates into Polygon, then Polygons class
    EL <- lapply(EL,Polygon,hole=FALSE)
    EL <- mapply(FUN=function(x,y) Polygons(list(x),y),EL,subby, SIMPLIFY = FALSE)
    
    #convert to Spatial Polygons class then create Unioned SpatialPolygons
    spEL <- SpatialPolygons(EL,1:length(EL),proj4string=projectionName)
    uSpEL <- gUnaryUnion(spEL)
    plot(uSpEL,add=TRUE)
    
    uSpELdf <- SpatialPolygonsDataFrame(uSpEL,data = data.frame(1))
    if (oerr) {
      writeOGR(uSpELdf,gisDir,layer = paste0(birdFold,"Error",projName),driver = "ESRI Shapefile",overwrite_layer = TRUE)
    }
    
    #create two data frames - one with all stopovers, the second with the combined breed, molt and winter stopovers
    psoTime <- pStatPeriod
    psoTime[1,"start"] <- 1
    psoTimerow <- which(psoTime$startDate < geo_df$Date.of.Initial.fitting[birdNum])
    psoTime$startDate[psoTimerow] <- geo_df$Date.of.Initial.fitting[psoTimerow]
    psoTime$start[psoTimerow] <- which(as.Date(quants$time) > as.Date(psoTime$startDate[psoTimerow]))[1]
    psoTime$startDate[psoTimerow] <- quants$time[psoTime$start[psoTimerow]]
    
    #spring
    soTime <- psoTime[springSO[[birdFold]],c("stopoverID","start","end","Duration","startDate","endDate","lengthGT14")]
    soTime[,"lengthGT14"] <- "spring"
    #breed
    if (!is.null(breedSO[[birdFold]])) {
      soTime <- rbind(soTime,psoTime[breedSO[[birdFold]][1],c("stopoverID","start","end","Duration","startDate","endDate","lengthGT14")])
      soTime[nrow(soTime),c("end","endDate","lengthGT14")] <- c(psoTime[breedSO[[birdFold]][length(breedSO[[birdFold]])],c("end","endDate")],"breed")
    }
    #molt
    if (!is.null(moltSO[[birdFold]])) {
      soTime <- rbind(soTime,psoTime[moltSO[[birdFold]][1],c("stopoverID","start","end","Duration","startDate","endDate","lengthGT14")])
      soTime[nrow(soTime),c("end","endDate","lengthGT14")] <- c(psoTime[moltSO[[birdFold]][length(moltSO[[birdFold]])],c("end","endDate")],"molt")
    }
    #fall
    if (!is.null(fallSO[[birdFold]])) {
      fsoTime <- psoTime[fallSO[[birdFold]],c("stopoverID","start","end","Duration","startDate","endDate","lengthGT14")]
      fsoTime[,"lengthGT14"] <- "fall"
      soTime <- rbind(soTime,fsoTime)
      rm(fsoTime)
    }
    #winter
    if (!is.null(winterSO[[birdFold]])) {
      soTime <- rbind(soTime,psoTime[winterSO[[birdFold]][1],c("stopoverID","start","end","Duration","startDate","endDate","lengthGT14")])
      soTime[nrow(soTime),c("end","endDate","lengthGT14")] <- c(psoTime[winterSO[[birdFold]][length(winterSO[[birdFold]])],c("end","endDate")],"winter")
    }
    # changes stopover twilght duration and set column name to "type" for typr of stopover (spring,fall,breed,etc)
    names(soTime)[7] <- "type"
    soTime$Duration <- soTime$end - soTime$start + 1
    
    #plot each stopover with different colors depending on time (spring,fall,breeding,molting,wintering)
    plot.so.color <- function (soTime,spEL,sC=sCol,bC=bCol,mC=mCol,fC=fCol,wC=wCol) {
      soDur <- soTime[2]:soTime[3]
      soErr <- spEL[as.character(soDur)]
      UsoErr <- gUnaryUnion(soErr)
      plot(UsoErr, add=TRUE, lwd=2)
      
      # require(scales)
      if (soTime[7] == "spring") {
        soCol <- sC
      } else if (soTime[7] == "breed") {
        soCol <- bC
      } else if (soTime[7] == "molt") {
        soCol <- mC
      } else if (soTime[7] == "fall") {
        soCol <- fC
      } else if (soTime[7] == "winter") {
        soCol <- wC
      }
      plot(soErr,col=adjustcolor(soCol,alpha.f = 0.05),border=FALSE,add=TRUE)
      
      coords <- c(mean(quants$Meanlon[soDur]),mean(quants$Meanlat[soDur]))
      return(list(soErr,UsoErr,
                  muCoord=c(mean(quants$Meanlon[soDur]),mean(quants$Meanlat[soDur]))))
    }
    
    #list for  non-stopover estimated error
    errList <- apply(soTime,1,FUN=function(x) plot.so.color(x,spEL))
    #list for stopover estimated error
    soErr <- lapply(errList,function(x) x[[1]])
    names(soErr) <- paste0("soErr",1:length(errList))
    list2env(soErr,env)
    #list for combined stopover estimated error - breeding, molting, wintering areas
    UsoErr <- lapply(errList,function(x) x[[2]])
    names(UsoErr) <- paste0("UsoErr",1:length(errList))
    list2env(UsoErr,env)
    
    soCoords <- sapply(errList,function(x) x[[3]])
    soCoordsdf <- SpatialPointsDataFrame(coords=matrix(c(soCoords[1,],soCoords[2,]),ncol=2),
                                         data = soTime, proj4string = projectionName)
    if (omupnt) {
      writeOGR(soCoordsdf,dsn = gisDir,layer = paste0(birdFold,"SoMeanPoints",projName),driver = "ESRI Shapefile",overwrite_layer = TRUE)
    }
    
    soErrT <- eval(parse(text=sprintf("bind(%s)",
                                      paste(paste0("soErr",1:nrow(soTime)),collapse=","))))
    soErrdf <- SpatialPolygonsDataFrame(soErrT,data = data.frame(twlID=unlist(mapply(seq,soTime$start,soTime$end)),
                                                                 Type=rep(soTime$type,soTime$Duration)))
    if (oso) {
      writeOGR(soErrdf,dsn = gisDir,layer = paste0(birdFold,"soError",projName),driver = "ESRI Shapefile",overwrite_layer = TRUE)
    }
    
    UsoErrT <- eval(parse(text=sprintf("bind(%s)",
                                       paste(paste0("UsoErr",1:nrow(soTime)),collapse=","))))
    UsoErrdf <- SpatialPolygonsDataFrame(UsoErrT,data = soTime, match.ID = FALSE)
    if (ouso) {
      writeOGR(UsoErrdf,dsn = gisDir,layer = paste0(birdFold,"UsoError",projName),driver = "ESRI Shapefile",overwrite_layer = TRUE)
    }
    
    mL <- Line(quants[subby,c("Meanlon","Meanlat")])
    mLs <- Lines(mL,"mean")
    msLs <- SpatialLines(list(mLs),proj4string = projectionName)    
    msLsdf <- SpatialLinesDataFrame(msLs, data = data.frame(1), match.ID = FALSE)
    
    if (opath) {
      writeOGR(msLsdf,dsn = gisDir,layer = paste0(birdFold,"meanPath",projName),driver = "ESRI Shapefile",overwrite_layer = TRUE)
    }
    
    twlMove <- unlist(apply(pMovePeriod,1,FUN=function(x) x[1]:x[2]))
    quantsMove <- quants[twlMove,]
    mLocdf <- SpatialPointsDataFrame(coords=quantsMove[,c("Meanlon","Meanlat")],
                                     data = quantsMove, proj4string = projectionName)
    if (omloc) {
      writeOGR(mLocdf,dsn = gisDir,layer = paste0(birdFold,"meanMVPoints",projName),driver = "ESRI Shapefile",overwrite_layer = TRUE)
    }
    
  }, error = function(err) {
    runError <- c(runError,paste("error on",birdFold,"(birdNum=",birdNum,")",err))
    print(paste("MY_ERROR: ",birdFold,err))
    if (!is.null(dev.list())) {dev.off()}
  })
  tryCatch({
    save(list=ls(),file=file.path(gisDir,paste0("StopoverMapping_values.RData")))
  }, error = function(err) {
    runError <- c(runError,paste("error on",birdFold,"(birdNum=",birdNum,")",err))
    print(paste("MY_ERROR: ",birdFold,err))
  })
  return(runError)
}

#run the mapping function for all birds in birdList
runError <- lapply(birdList,get.maps,ol=TRUE,
                   projName="wgs",projectionName=CRS("+init=epsg:4326"))

ulError <- unlist(runError)


# Timeline Plotting ####
## this will only plot the first year if there is more than one, so you have to run it again for subsequent year (if needed)
##it will also try to match the dates for each bird so they are the same x axes
plot.stopover.timeline <- function(birdList=fileList[index.bird], passNum=3,
                                   plotType=c("reg","cvert","choriz"),pNames=NULL,plotNesting=FALSE,
                                   sCol=sColor,bCol=bColor,mCol=mColor,fCol=fColor,wCol=wColor,
                                   springSO=sso,breedSO=bso,moltSO=mso,fallSO=fso,winterSO=wso,
                                   speciesWD=spDir, geolocatorDF=geo_df,
                                   species=NULL,calp=cal.period, lcorr=l.correct) {
  
  require(lubridate)
  runError <- c()
  tryCatch({
    if (is.null(species)) {
      species=unique(geolocatorDF$Species)
      print("species pulled from geolocator data frame")
    }
    
    plotType <- plotType[1]
    # set plotting values for vertical and horizontal legend for combined plots
    if (plotType != "reg") {
      density <- c(20,NA,0,20,30)
      color <- c("black","black",NA,"black","black")
      angle <- c(45,NA,NA,-45,90)
    }
    
    if (plotType == "reg") {
      par(mfrow=c(1,1),mar=c(5,4,4,2))
      #set plotting values
      density <- NA
      color <- c(sC,bC,mC,fC,wC)
      angle <- NA
    } else if (plotType == "cvert") {
      # layout(rbind(c(1,6),c(2,6),c(3,6),c(4,6),c(5,6)),
      layout(matrix(data=c(seq(1,length(birdList)),rep.int(4,length(birdList))),ncol=2),
             widths = c(2.5,1))
      par(mar=c(2,1,1,2))
    } else if (plotType == "choriz") {
      par(mfrow=c(6,1),mar=c(2,1,1,2))
      
      #Plot Horizontal Legend
      plot(quants$time,rep(0,length(quants$time)),
           xlab="", ylab="", yaxt="n",xaxt="n",
           xlim = c(0,30),
           ylim=c(0,8),
           type="n", main = "Legend",bty="n")
      legend(x=0,y=10,
             legend = c("Spring Stopovers","Breeding Area","Molting Area",
                        "Fall Stopovers","Winter Area"),
             density = density,
             angle = angle,
             col=color,
             cex=1.5,
             border = "black",bty="n",
             horiz=TRUE)
    } else {
      stop("Plotting type is not a recognized value")
    }
    
    #function for plotting timeline
    plot.timeline.sub <- function(birdFold,passNum=3,pty=plotType,pNames=birdList,pNest=plotNesting,
                                  col=color,den=density,ang=angle,
                                  ss=springSO,bs=breedSO,ms=moltSO,fs=fallSO,ws=winterSO,
                                  bL=birdList,spWD=speciesWD, geoDF=geolocatorDF,
                                  spec=species,cp=calp, lc=lcorr) {
      
      runError <- c()
      tryCatch({
        if (is.null(ss)) {ss <- NA}
        if (is.null(bs)) {bs <- NA}
        if (is.null(ms)) {ms <- NA}
        if (is.null(fs)) {fs <- NA}
        if (is.null(ws)) {ws <- NA}
        
        m.folder <- paste(cp,"_",lc,sep="")
        passFold <- paste0("Movement",passNum)
        passDir <- file.path(spWD,birdFold,m.folder,passFold)
        stopoverFold <- paste0("stopover_",passNum,"analysis")
        
        load(file.path(passDir,stopoverFold,paste0("Stopover",passNum,"analysis.RData")))
        birdNum <- which(bL %in% birdFold)
        quants <- Result$Results$Quantiles
        
        twlStart <- which(as.Date(quants$time) > as.Date(geo_df$Date.of.Initial.fitting[birdNum]))[1]
        pStatPeriod$startDate[pStatPeriod$start < twlStart] <- quants$time[twlStart]
        pStatPeriod$start[pStatPeriod$start < twlStart] <- twlStart
        
        if (nrow(pStatPeriod) < max(c(ss,bs,ms,fs,ws),na.rm=TRUE)) {
          plot(quants$time,rep(0,length(quants$time)), type="n",
               main = "Error in number of stopovers",bty="n")
          stop("Indicated stopovers types excede the number of stopovers - You say there are more than the data show")
        }
        
        ticks <- c(as.POSIXct(paste0(year(min(quants$time)),"-01-01"),tz="GMT"),
                   as.POSIXct(paste0(year(min(quants$time)),"-02-01"),tz="GMT"),
                   as.POSIXct(paste0(year(min(quants$time)),"-03-01"),tz="GMT"),
                   as.POSIXct(paste0(year(min(quants$time)),"-04-01"),tz="GMT"),
                   as.POSIXct(paste0(year(min(quants$time)),"-05-01"),tz="GMT"),
                   as.POSIXct(paste0(year(min(quants$time)),"-06-01"),tz="GMT"),
                   as.POSIXct(paste0(year(min(quants$time)),"-07-01"),tz="GMT"),
                   as.POSIXct(paste0(year(min(quants$time)),"-08-01"),tz="GMT"),
                   as.POSIXct(paste0(year(min(quants$time)),"-09-01"),tz="GMT"),
                   as.POSIXct(paste0(year(min(quants$time)),"-10-01"),tz="GMT"),
                   as.POSIXct(paste0(year(min(quants$time)),"-11-01"),tz="GMT"),
                   as.POSIXct(paste0(year(min(quants$time)),"-12-01"),tz="GMT"),
                   as.POSIXct(paste0(year(min(quants$time))+1,"-01-01"),tz="GMT"))
        
        ## Plot the timeline with median start and end times for stopover locations
        if (pty == "reg") {
          ylim <- c(min(statPeriod$Distance.Moved),
                    max(statPeriod$Distance.Moved))
          yb <- lylim - 1
          yt <- uylim + 1
          
        } else {
          ylim <- c(-2,8)
          yb <- 0
          yt <- 5
        }
        
        plot(quants$time,rep(0,length(quants$time)),
             xlab="Date", ylab="", yaxt="n", xaxt="n",
             xlim = c(as.POSIXct(paste0(year(min(quants$time)),"-01-15"),tz="GMT"),
                      as.POSIXct(paste0(year(min(quants$time)),"-12-31"),tz="GMT")),
             ylim=ylim,
             type="n", main = pNames[birdNum],bty="n")
        # plots x axis tickmarks
        if (birdNum == par()$mfrow[1] | par()$mfrow[1] == 1) {
          axis.POSIXct(1,at=ticks,lwd=0,lwd.ticks=1)
        }
        
        if (pty != "reg") {
          sDate <- which(quants$time > as.Date(geoDF$Date.of.Initial.fitting[grep(sub(paste0(spec,"_"),"",birdFold),geoDF$Metal)]))[1]
          lines(x=c(quants$time[sDate], quants$time[nrow(quants)]),y=c(0,0),lwd=3)
        }
        
        #spring stopover
        rect(xleft = pStatPeriod$startDate[ss], xright = pStatPeriod$endDate[ss],
             ybottom = yb,  ytop = yt, density = den[1], angle = ang[1], col=col[1])
        #breeding area
        rect(xleft = pStatPeriod$startDate[bs[1]], xright = pStatPeriod$endDate[bs[length(bs)]],
             ybottom = yb,  ytop = yt, density = den[2],angle = ang[2], col=col[2])
        #molt/stage stopover
        rect(xleft = pStatPeriod$startDate[ms[1]], xright = pStatPeriod$endDate[ms[length(ms)]],
             ybottom = yb,  ytop = yt, density = den[3],angle = ang[3], col=col[3])
        #fall stopover
        rect(xleft = pStatPeriod$startDate[fs], xright = pStatPeriod$endDate[fs],
             ybottom = yb,  ytop = yt,density = den[4], angle = ang[4], col=col[4])
        #winter area
        rect(xleft = pStatPeriod$startDate[ws[1]], xright = pStatPeriod$endDate[ws[length(ws)]],
             ybottom = yb,  ytop = yt, density = den[5], angle = ang[5], col=col[5],
             lty=3)
        rect(xleft = pStatPeriod$startDate[ws[1]], xright = pStatPeriod$endDate[ws[length(ws)]],
             ybottom = yb,  ytop = yt, border=col[5])
      }, error = function(err) {
        runError <- c(runError,paste("error on",birdFold,"(birdNum=",birdNum,")",err))
        print(paste("MY_ERROR: ",birdFold,err))
      })
      return(runError)
    }
    
    rError <- mapply(a=birdList,b=springSO,c=breedSO,d=moltSO,e=fallSO,f=winterSO,
                     g=passNum,h=list(pNames),
                     FUN=function(a,b,c,d,e,f,g,h) plot.timeline.sub(birdFold=a,ss=b,bs=c,ms=d,fs=e,ws=f,
                                                                     passNum=g,pNames=h))
    runError <- unlist(rError)
    if (plotType == "cvert") {
      #vertical legend
      plot(1,1,
           xlab="", ylab="", yaxt="n",xaxt="n",
           xlim = c(0,10),
           ylim=c(0,30),
           type="n", main = "Legend",bty="n")
      legend(x=0,y=30,
             legend = c("Spring Stopovers","Breeding Area","Molting Area",
                        "Fall Stopovers","Winter Area","Geolocator active"),
             density = c(density,0),
             angle = c(angle,0),
             col=c(color,"black"),
             lty = c(NA,NA,NA,NA,NA,1),
             lwd=c(NA,NA,NA,NA,NA,3),
             cex=1.5,
             # border = "black",
             bty="n",
             xjust=0,
             y.intersp=3)
    }
    
  }, error = function(err) {
    runError <- c(runError,paste("error on",species,err))
    print(paste("MY_ERROR: ",species,err))
  })
  return(runError)
}

runError <- plot.stopover.timeline(birdList=fileList[index.bird], passNum=3,plotType="cvert",
                                   # pNames=c("Canada Male #1","Canada Male #2","Alaska Male #1",
                                   #          "Canada Female #1","Canada Female #2"))
                                   pNames=c("Canada Male #3","Canada Female #3","Canada Male #4"))






