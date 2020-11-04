#author: Dominic Fawcett
#last modified: 12/10/2020
#description: script to extract NDVI time series per crown and fit logistic function  to retrieve phenological parameters
#             generates boxplots as shown in publication, can be modified to include evergreen crowns and use no brightest pixel masking
#NOTE: variables SOS, MOS and EOS refer to SOS, MOG and SOP as presented in the manuscript
#NOTE: adjust paths for Shapefiles, NDVI and green-band orthomosaics

library(rgdal)
library(lubridate)
library(RColorBrewer)
library(ggplot2)
library(raster)

setwd('D:/Microclimate/Sequoia/combNDVI')#change directory here to location of NDVI orthomosaics
filelist <- list.files(getwd())
filelistGRE <- list.files('D:/Microclimate/Sequoia/combGRE')#change directory here to location of green-band orthomosaics

refShapes <- readOGR(dsn = 'D:/Microclimate/GPS', layer = "TreeClassPolysreprojExtNoZeroNEW")#change directory here to location of tree shapefile

#tree species and reference surface class names
classnames <- c('Hawthorn','Oak','Goat Willow','Grey Willow','Holly','Elm','Pine','Gorse','Sycamore','Chestnut','Lime','Alder','Beech','Ash','Apple','Elder','Spruce','Concrete')

NDVIstack <- stack()
GREstack <- stack()
NDVIlist <- list() 
datelist <- list()

#get dates
for (i in 1:length(filelist)){
  stringpart <-  substr(filelist[[i]],11,20)
  datelist <- append(datelist,as.Date(stringpart,"%Y_%m_%d"))
}

#get doys
doylist <- yday(datelist)
doystart <- doylist[1]
doyend <- doylist[length(doylist)]

#resample datasets for consistent extents
for (i in 1:length(filelist)){
  setwd('D:/Microclimate/Sequoia/combNDVI')
  NDVIraster <- raster(filelist[[i]])
  setwd('D:/Microclimate/Sequoia/combGRE')
  GREraster <- raster(filelistGRE[[i]])
  if(i==1){
    GREraster <- resample(GREraster,NDVIraster)
  }
  if(i>1){
    NDVIraster <- resample(NDVIraster,NDVIstack[[1]])
    GREraster <- resample(GREraster,NDVIstack[[1]])
  }
  NDVIstack <- stack(NDVIstack,NDVIraster)
  GREstack <- stack(GREstack,GREraster)
}


refShapesReproj <- spTransform(refShapes,crs(NDVIraster))

#extract 80th percentile values per tree from green band
treesGRE80thpercentile <- raster::extract(GREstack,refShapesReproj,fun=function(X,na.rm=TRUE){quantile(X,0.8,na.rm=TRUE)}, df=TRUE)

treesNDVI <- raster::extract(NDVIstack,refShapesReproj, df=TRUE)
treesGRE <- raster::extract(GREstack,refShapesReproj, df=TRUE)

treesNDVITS <- matrix(nrow=nrow(treesGRE80thpercentile),ncol=21)
treesNDVISDTS <- matrix(nrow=nrow(treesGRE80thpercentile),ncol=21)

#perform brightest pixel masking of NDVI values for reduced background effect
for (i in 1:nrow(treesGRE80thpercentile)){
  
  treesNDVITS[i,1] <- i
  for(j in 1:20){
    GREmask <- treesGRE[treesGRE[,1]==i,j+1]>treesGRE80thpercentile[i,j+1]
    NDVIvals <- treesNDVI[treesNDVI[,1]==i,j+1]
    NDVIgtGRE80pc <- NDVIvals[GREmask]#replace with NDVIvals to disble brightest pixel masking
    treesNDVITS[i,j+1] <- mean(NDVIgtGRE80pc,na.rm=TRUE)
    treesNDVISDTS[i,j+1] <- sd(NDVIgtGRE80pc,na.rm=TRUE)
    
  }
}

#prepare empty lists and index vectors
treeclasses <- refShapes@data$class
treeclassnrs <- c(seq(1:17),99)
linecols <- colorRampPalette(brewer.pal(9,"Set1"))(18)
linecolDF <- data.frame(class=treeclassnrs,linecols=linecols)
ndviplotslist <- list()
SOSdoyplusSD <- list()
SOSdoyminusSD <- list()
slp <- list()
midpointlist<- list()
SOSdoylist <- matrix(NA,nrow=nrow(treesNDVITS),ncol=4)
residualSTDerror <- matrix(NA,nrow=nrow(treesNDVITS),ncol=2)
indcount <- 1

#perform model fitting and extract SOS, MOS and EOS values from fitted model

#iterate through classes
for(i in 1:18){
  #i=1#debug

  SOSdoy <- list()
  classNDVIvals <- treesNDVITS[treeclasses==treeclassnrs[i],2:21]
  
  if(is.null(nrow(classNDVIvals))){
  ntrees <- 1
  }else{
  ntrees <- nrow(classNDVIvals)
  }
  doylistcurrent <- doylist
  
  SOSdoymat <- matrix(NA,nrow=1,ncol=ntrees)
  
  if(i %in% c(5,7,8,17,18)){#disabled fitting for evergreen and reference surface
    SOSdoylist[indcount,1] <- NA
    SOSdoylist[indcount,2] <- NA
    SOSdoylist[indcount,3] <- NA
    SOSdoylist[indcount,4] <- i
    indcount <- indcount+1
  }else{
  
#iterate through individuals per class
  for(j in c(1:ntrees)){

     
    if(i==1){#remove hawthorn flowering observations for successful fitting
      doylistcurrent <- doylist[-10:-13]
      ndvidat <- unlist(classNDVIvals[j,-10:-13])
    }else if(ntrees==1){
      ndvidat <- unlist(classNDVIvals)
    }else{
      ndvidat <- unlist(classNDVIvals[j,])
    }
    if(i %in% c(13)){#different inital parameters for beech for successful fitting
    fitmodel <- nls(ndvidat ~ a/(1 + exp(-b * (doylistcurrent-c)) ) + d,control=nls.control(warnOnly=TRUE,maxiter = 50), start=list(a=0.25,b=0.2,c=100,d=0.6))
    }else{
    tt <- tryCatch(nls(ndvidat ~ a/(1 + exp(-b * (doylistcurrent-c)) ) + d,control=nls.control(warnOnly=TRUE,maxiter = 50), start=list(a=0.25,b=0.1,c=120,d=0.6)),error=function(e) e, warning=function(w) w)
    fitmodel <- nls(ndvidat ~ a/(1 + exp(-b * (doylistcurrent-c)) ) + d,control=nls.control(warnOnly=TRUE,maxiter = 50), start=list(a=0.25,b=0.1,c=120,d=0.6))
    }
    if(is(tt,'warning')){
      print(paste('fail at class ',i, ' individual ',j))
      SOSdoylist[indcount,1:3] <- NA
      SOSdoylist[indcount,4] <- i
    }else{
      
    #get model coefficients
    modsum <- summary(fitmodel)
    a <- modsum$parameters[1]
    b <- modsum$parameters[2]
    c <- modsum$parameters[3]
    d <- modsum$parameters[4]
    preddoylist <- seq(doystart,doyend,1)
    predNDVImod <- predict(fitmodel,data.frame(doylistcurrent=preddoylist))
    

    #slope
    firstDerivNDVImod <- (a*b*exp(-b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^2
    #curvature
    secondDerivNDVImod <- a*(((2*b^2*exp(-2*b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^3)-((b^2*exp(-b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^2))
    
    curvature <- secondDerivNDVImod/((1+(firstDerivNDVImod)^2)^1.5)
    
    ROCcurvature <- c(curvature[2:(doyend-doystart+1)],NA)-curvature
    
    if(length((which(diff(sign(diff(ROCcurvature)))==-2)+1))<2){ #if SOS can't be retrieved, ignore MOS and EOS
    print(paste('SOS fail at class', i,', individual ',j))
    SOS <- NA
    MOS <- NA
    EOS <- NA
    }else{ 
    SOS <- preddoylist[(which(diff(sign(diff(ROCcurvature)))==-2)+1)[1]]#get DOY of ROC local maximum 1
    MOS <- preddoylist[firstDerivNDVImod==max(firstDerivNDVImod,na.rm=TRUE)]
    EOS <- preddoylist[(which(diff(sign(diff(ROCcurvature)))==-2)+1)[2]]#get DOY of ROC local maximum 2
    }
    
    SOSdoylist[indcount,1] <- SOS

    SOSdoylist[indcount,2] <- MOS
   
    SOSdoylist[indcount,3] <- EOS
    SOSdoylist[indcount,4] <- i

    residualSTDerror[indcount,1] <- modsum$sigma
    residualSTDerror[indcount,2] <- i
    indcount <- indcount+1
  
    }
  }
  }
}

#make data frame with extracted dates per individual and class
SOSdf <- data.frame(SOS=SOSdoylist[,1],MOS=SOSdoylist[,2],EOS=SOSdoylist[,3],class=SOSdoylist[,4],RSE=residualSTDerror[,1])

#aggregate species / types
SOSclassMedians80perc <- aggregate(.~class,SOSdf,median,na.rm=T)
SOSclassMeans80perc <- aggregate(.~class,SOSdf,mean,na.rm=T)
SOSclassCount80perc <- aggregate(.~class,SOSdf,length)

#names of analysed classes 
boxplotnames <- classnames[c(-5,-7,-8,-17,-18)]

classnamesdf <- data.frame(class=seq(1,18,1),names = classnames,colours=linecols)




par(mar=c(6,5,2,3))
#SOS

new_order <- with(SOSdf, reorder(SOSdf$class, SOSdf$SOS, median , na.rm=T))

ordernrdf <- data.frame(class=as.numeric(levels(new_order)))

reorderedclassdf <- merge(ordernrdf,classnamesdf,by='class',sort=FALSE)

bpSOS <- boxplot(SOSdf$SOS~new_order,names=reorderedclassdf$names[1:13],las=2,col=as.vector(reorderedclassdf$colours[1:13]),drop=TRUE,ylab='SOS [DOY]',outline=T)

#MOS

new_order <- with(SOSdf, reorder(SOSdf$class, SOSdf$MOS, median , na.rm=T))

ordernrdf <- data.frame(class=as.numeric(levels(new_order)))

reorderedclassdf <- merge(ordernrdf,classnamesdf,by='class',sort=FALSE)

bpMOS <- boxplot(SOSdf$MOS~new_order,names=reorderedclassdf$names[1:13],las=2,col=as.vector(reorderedclassdf$colours[1:13]),drop=TRUE,ylab='MOG [DOY]',outline=T)

#EOS

new_order <- with(SOSdf, reorder(SOSdf$class, SOSdf$EOS, median , na.rm=T))

ordernrdf <- data.frame(class=as.numeric(levels(new_order)))

reorderedclassdf <- merge(ordernrdf,classnamesdf,by='class',sort=FALSE)

bpEOS <- boxplot(SOSdf$EOS~new_order,names=reorderedclassdf$names[1:13],las=2,col=as.vector(reorderedclassdf$colours[1:13]),drop=TRUE,ylab='SOP [DOY]',outline=T)

######
#mean standard errors of residuals per class

stderrdf <- data.frame(stderr=residualSTDerror[,1],class=residualSTDerror[,2])
errperclass <- aggregate(stderrdf$stderr,by=list(stderrdf$class),mean)


