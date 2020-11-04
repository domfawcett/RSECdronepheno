library(rgdal)
library(lubridate)
library(RColorBrewer)
library(ggplot2)
library(raster)
setwd('D:/Microclimate/Sequoia/combNDVI')
filelist <- list.files(getwd())
filelistGRE <- list.files('D:/Microclimate/Sequoia/combGRE')

refShapes <- readOGR(dsn = 'D:/Microclimate/GPS', layer = "TreeClassPolysreprojExtNoZeroNEW")

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
  #strdatelist <- append(strdatelist,substr(filelist[[i]],11,20))
}


refShapesReproj <- spTransform(refShapes,crs(NDVIraster))

#writeOGR(refShapesReproj,dsn= 'D:/Microclimate/GPS/treereprojtest.shp',layer='treereprojtest',driver='ESRI Shapefile')

#treesNDVITS <- extract(NDVIstack,refShapesReproj,fun=mean, df=TRUE)
#treesNDVISDTS<- extract(NDVIstack,refShapesReproj,fun=sd, df=TRUE)

treesGRE80thpercentile <- raster::extract(GREstack,refShapesReproj,fun=function(X,na.rm=TRUE){quantile(X,0.8,na.rm=TRUE)}, df=TRUE)

treesNDVI <- raster::extract(NDVIstack,refShapesReproj, df=TRUE)
treesGRE <- raster::extract(GREstack,refShapesReproj, df=TRUE)

treesNDVITS <- matrix(nrow=nrow(treesGRE80thpercentile),ncol=21)
treesNDVISDTS <- matrix(nrow=nrow(treesGRE80thpercentile),ncol=21)


for (i in 1:nrow(treesGRE80thpercentile)){
  
  treesNDVITS[i,1] <- i
  for(j in 1:20){
    GREmask <- treesGRE[treesGRE[,1]==i,j+1]>treesGRE80thpercentile[i,j+1]
    NDVIvals <- treesNDVI[treesNDVI[,1]==i,j+1]
    NDVIgtGRE80pc <- NDVIvals[GREmask]
    treesNDVITS[i,j+1] <- mean(NDVIgtGRE80pc,na.rm=TRUE)
    treesNDVISDTS[i,j+1] <- sd(NDVIgtGRE80pc,na.rm=TRUE)
    
  }
}


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
  
  if(i %in% c(5,7,8,17,18)){
    SOSdoylist[indcount,1] <- NA
    SOSdoylist[indcount,2] <- NA
    SOSdoylist[indcount,3] <- NA
    SOSdoylist[indcount,4] <- i
    indcount <- indcount+1
  }else{
  
  for(j in c(1:ntrees)){
    #debug
   # j=24
    #j=30
   #j=9
    #doylist2 <- doylist[1:11]
    #ndvidat2 <- ndvidat[1:11]
     
    if(i==1){
      doylistcurrent <- doylist[-10:-13]
      ndvidat <- unlist(classNDVIvals[j,-10:-13])
    }else if(ntrees==1){
      ndvidat <- unlist(classNDVIvals)
    }else{
      ndvidat <- unlist(classNDVIvals[j,])
    }
    if(i %in% c(13)){
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
    
    #DEBUG
    plot(doylistcurrent,ndvidat)
    points(preddoylist,predNDVImod)
    #slope
    firstDerivNDVImod <- (a*b*exp(-b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^2
    #curvature
    secondDerivNDVImod <- a*(((2*b^2*exp(-2*b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^3)-((b^2*exp(-b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^2))
    
    curvature <- secondDerivNDVImod/((1+(firstDerivNDVImod)^2)^1.5)
    
    ROCcurvature <- c(curvature[2:(doyend-doystart+1)],NA)-curvature
    
    if(length((which(diff(sign(diff(ROCcurvature)))==-2)+1))<2){
    print(paste('SOS fail at class', i,', individual ',j))
    SOS <- NA
    MOS <- NA#preddoylist[firstDerivNDVImod==max(firstDerivNDVImod,na.rm=TRUE)]#ROCcurvature ==min(ROCcurvature,na.rm=TRUE)][1]
    EOS <- NA#preddoylist[(which(diff(sign(diff(ROCcurvature)))==-2)+1)[1]]#preddoylist[secondDerivNDVImod==min(secondDerivNDVImod)]
    }else{ 
    SOS <- preddoylist[(which(diff(sign(diff(ROCcurvature)))==-2)+1)[1]]#preddoylist[secondDerivNDVImod==max(secondDerivNDVImod)]
    MOS <- preddoylist[firstDerivNDVImod==max(firstDerivNDVImod,na.rm=TRUE)]#ROCcurvature ==min(ROCcurvature,na.rm=TRUE)][1]
    EOS <- preddoylist[(which(diff(sign(diff(ROCcurvature)))==-2)+1)[2]]#preddoylist[secondDerivNDVImod==min(secondDerivNDVImod)]
    }
    
    SOSdoylist[indcount,1] <- SOS
    #if(MOS>150|MOS<90){
    #SOSdoylist[indcount,2] <- NA 
    #print(paste('MOS outlier removed at class',i,', individual ',j))
    #}else{
    SOSdoylist[indcount,2] <- MOS
    #}
    SOSdoylist[indcount,3] <- EOS
    SOSdoylist[indcount,4] <- i

    residualSTDerror[indcount,1] <- modsum$sigma
    residualSTDerror[indcount,2] <- i
    indcount <- indcount+1
    #slp <- append(slp,max(slopeNDVImod))
    }
  }
  }
}

SOSdf <- data.frame(SOS=SOSdoylist[,1],MOS=SOSdoylist[,2],EOS=SOSdoylist[,3],class=SOSdoylist[,4],RSE=residualSTDerror[,1])

#aggregate species / types
SOSclassMedians80perc <- aggregate(.~class,SOSdf,median,na.rm=T)
SOSclassMeans80perc <- aggregate(.~class,SOSdf,mean,na.rm=T)
SOSclassCount80perc <- aggregate(.~class,SOSdf,length)
SOSclassMeansOrig <- aggregate(.~class,SOSdf,mean,na.rm=T)


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

#SOSclassMeansOrig <- SOSclassMedians80perc

SOSdifmat <- SOSclassMeansOrig-SOSclassMeans80perc
#########
#midpoint, non-interpolated version
  
# 
#   test <- unlist(classNDVIvals)
# 
#   classmeanNDVIvals <- colMeans(classNDVIvals,na.rm=TRUE)
#   classSDNDVIvals <- apply(classNDVIvals,2,sd,na.rm=TRUE)
#   NDVIdata <- data.frame(doy=doylist,ndvi=classmeanNDVIvals,sderror=classSDNDVIvals)
# 
#   #if evergreen, don't do phenology (could assess evergreenness with descriptor, here done with known classes)
#   if(i %in% c(5,7,17,18)){
#   midpointlist[[i]] <- NA
#   SOSdoy[[i]] <- NA
#   SOSdoyplusSD[[i]] <- NA
#   SOSdoyminusSD[[i]] <- NA
#   slp[[i]] <- NA
#   }else{
#   midpoint <- min(NDVIdata$ndvi,na.rm=TRUE)+(max(NDVIdata$ndvi,na.rm=TRUE)-min(NDVIdata$ndvi,na.rm=TRUE))/2
#   midpointplusSD <- midpoint+mean(NDVIdata$sderror,na.rm=TRUE)
#   midpointminusSD <- midpoint-mean(NDVIdata$sderror,na.rm=TRUE)
#   midpointlist[[i]] <- midpoint
#   #approximate doy at midpoint NDVI
#   #constraint of using only first half of data for vegetation with flowering stage (hawthorne, Gorse)
#   if(i %in% c(1)){
#   interpMid<- approx(NDVIdata$ndvi[1:9],NDVIdata$doy[1:9],xout=midpoint)
#   SOSdoy[[i]] <- interpMid$y
#   SOSdoyplusSD[[i]] <- approx(NDVIdata$ndvi[1:9],NDVIdata$doy[1:9],xout=midpointplusSD)$y
#   SOSdoyminusSD[[i]] <- approx(NDVIdata$ndvi[1:9],NDVIdata$doy[1:9],xout=midpointminusSD)$y
#   }else if(i %in% c(8)){
#     interpMid<- approx(NDVIdata$ndvi[1:11],NDVIdata$doy[1:11],xout=midpoint)
#     SOSdoy[[i]] <- interpMid$y
#     SOSdoyplusSD[[i]] <- approx(NDVIdata$ndvi[1:11],NDVIdata$doy[1:11],xout=midpointplusSD)$y
#     SOSdoyminusSD[[i]] <- approx(NDVIdata$ndvi[1:11],NDVIdata$doy[1:11],xout=midpointminusSD)$y
#   }else{
#     interpMid<- approx(NDVIdata$ndvi,NDVIdata$doy,xout=midpoint)
#     SOSdoy[[i]] <- interpMid$y
#     SOSdoyplusSD[[i]] <- approx(NDVIdata$ndvi,NDVIdata$doy,xout=midpointplusSD)$y
#     SOSdoyminusSD[[i]] <- approx(NDVIdata$ndvi,NDVIdata$doy,xout=midpointminusSD)$y
#   } 
# 
#   #approximate slope (better method?)
#   delta<- approx(NDVIdata$doy,NDVIdata$ndvi,xout=interpMid$y+1)$y-approx(NDVIdata$doy,NDVIdata$ndvi,xout=interpMid$y-1)$y
#   slp[[i]] <- delta/2
#   }
#   lines(doylist,classmeanNDVIvals,col=linecols[i])
#   if(i<18){
#   p <- ggplot(NDVIdata, aes(x = doy, y = ndvi)) +
#     ggtitle(paste(classnames[i],', n=',nrow(classNDVIvals),sep=''))+
#     theme_classic() +
#     xlab('DOY') +
#     ylab('NDVI') +
#     ylim(c(0.5,1))+
#     geom_line(col=linecols[i]) +
#     geom_ribbon(aes(ymin = ndvi - sderror,
#                     ymax = ndvi + sderror), alpha = 0.2,fill=linecols[i])
#   }else{
#     p <- ggplot(NDVIdata, aes(x = doy, y = ndvi)) +
#       ggtitle(paste(classnames[i],', n=',nrow(classNDVIvals),sep=''))+
#       theme_classic() +
#       xlab('DOY') +
#       ylab('NDVI') +
#       ylim(c(0,0.5))+
#       geom_line(col=linecols[i]) +
#       geom_ribbon(aes(ymin = ndvi - sderror,
#                       ymax = ndvi + sderror), alpha = 0.2,fill=linecols[i])
#   }
#   ndviplotslist[[i]] <- p
# }

#plot the SOS vs slope of every class
plot(slp,SOSdoy,type='n',ylim=c(100,160),xlab='NDVI slope at midpoint',ylab='SOS DOY')
for(i in 1:18){
points(slp[[i]],SOSdoy[[i]],pch=16,col=linecols[i])
arrows(slp[[i]], SOSdoyminusSD[[i]], slp[[i]], SOSdoyplusSD[[i]], length=0.05, angle=90, code=3,col=linecols[i])
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

multiplot(plotlist = ndviplotslist, cols=6 )


classnamesNewLegend <- c('Hawthorne','Chestnut','Poplar','Alder','Beech','Ash','Apple','Elder','Fir2','Oak','Goat Willow','Grey Willow','Holly','Elm','Fir','Gorse','Sycamore','Concrete')

########
#Visualise classes
colRegions <- as.vector(linecolDF$linecols[match(levels(refShapesReproj$class), linecolDF$class)])
spplot(refShapesReproj, zcol = "class", col.regions = colRegions,
       colorkey = list(labels = list( labels = classnamesNewLegend,at=c(1:18),        
                                      width = 1, cex = 1)))

