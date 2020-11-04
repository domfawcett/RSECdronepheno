#author: Dominic Fawcett
#last modified: 12/10/2020
#description: script to extract NDVI time series per crown and make ribbon plots per species and concrete reference surface
#NOTE: modify paths to shapefiles, NDVI and green-band orthomosaics


library(raster)
library(rgdal)
library(lubridate)
library(RColorBrewer)
library(ggplot2)

setwd('D:/Microclimate/Sequoia/combNDVI')#modify path to NDVI orthomosaics
filelist <- list.files(getwd())
filelistGRE <- list.files('D:/Microclimate/Sequoia/combGRE')#modify path to green-band orthomosaics

#read digitised crown extents
refShapes <- readOGR(dsn = 'D:/Microclimate/GPS', layer = "TreeClassPolysreprojExtNoZeroNEW")#modify path to tree shapefile

#class names
classnames <- c('Hawthorn','Oak','Goat Willow','Grey Willow','Holly','Elm','Pine','Gorse','Sycamore','Chestnut','Lime','Alder','Beech','Ash','Apple','Elder','Spruce','Concrete')

NDVIstack <- stack()
GREstack <- stack()
NDVIlist <- list()
datelist <- list()

#get dates from filenames
for (i in 1:length(filelist)){
  stringpart <-  substr(filelist[[i]],11,20)
  datelist <- append(datelist,as.Date(stringpart,"%Y_%m_%d"))
}

#get doys
doylist <- yday(datelist)
doystart <- doylist[1]
doyend <- doylist[length(doylist)]

#stack all NDVI rasters
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

treesNDVISDTS<- extract(NDVIstack,refShapesReproj,fun=sd, df=TRUE)

#extract 80th percentile of green band for masking
treesGRE80thpercentile <- extract(GREstack,refShapesReproj,fun=function(X,na.rm=TRUE){quantile(X,0.8,na.rm=TRUE)}, df=TRUE)

treesNDVI <- extract(NDVIstack,refShapesReproj, df=TRUE)
treesGRE <- extract(GREstack,refShapesReproj, df=TRUE)

treesNDVITS <- matrix(nrow=nrow(treesGRE80thpercentile),ncol=21)

#extract NDVI values per tree, toggle brightest pixel masking
for (i in 1:nrow(treesGRE80thpercentile)){
  
  treesNDVITS[i,1] <- i
  for(j in 1:20){
    GREmask <- treesGRE[treesGRE[,1]==i,j+1]>treesGRE80thpercentile[i,j+1]
    NDVIvals <- treesNDVI[treesNDVI[,1]==i,j+1]
    NDVIgtGRE80pc <- NDVIvals[GREmask]
    treesNDVITS[i,j+1] <- mean(NDVIgtGRE80pc,na.rm=TRUE) #toggle to add brightest pixel masking
    #treesNDVITS[i,j+1] <- mean(NDVIvals,na.rm=TRUE) #toggle to remove brightest pixel masking
  }
}


#prepare vectors and lists for pheno extraction
treeclasses <- refShapes@data$class
treeclassnrs <- c(seq(1:17),99)
linecols <- colorRampPalette(brewer.pal(9,"Set1"))(18)
linecolDF <- data.frame(class=treeclassnrs,linecols=linecols)
ndviplotslist <- list()
SOSdoyplusSD <- list()
SOSdoyminusSD <- list()
slp <- list()
midpointlist<- list()
SOSdoylist <- matrix(NA,nrow=nrow(treesNDVITS),ncol=2)

ndvistart80perc <- matrix(NA,nrow=13,ncol=1)
ndviend80perc <- matrix(NA,nrow=13,ncol=1)
#########
declist <- c(1,2,3,4,6,9,10,11,12,13,14,15,16)#list of deciduous 

#extract NDVI means and SD per class

for(i in 1:18){
  classNDVIvals <- treesNDVITS[treeclasses==treeclassnrs[i],2:21]
  
  if(is.vector(classNDVIvals)){
    classmeanNDVIvals <- classNDVIvals
    classSDNDVIvals <- 0
    size <- 1
  }else{
  classmeanNDVIvals <- colMeans(classNDVIvals,na.rm=TRUE)
  classSDNDVIvals <- apply(classNDVIvals,2,sd,na.rm=TRUE)
  size <- nrow(classNDVIvals)
  }
  NDVIdata <- data.frame(doy=doylist,ndvi=classmeanNDVIvals,sderror=classSDNDVIvals)
  ndvistart80perc[i] <- mean(classmeanNDVIvals[1:3])#for leaf-off / leaf-on difference calculation
  ndviend80perc[i] <- mean(classmeanNDVIvals[17:20])#for leaf-off / leaf-on difference calculation
  
  if(i<18){#non-reference surface plots
  p <- ggplot(NDVIdata, aes(x = doy, y = ndvi)) + #vegetation crown plot
    ggtitle(paste(classnames[i],', n=',size,sep=''))+
    theme_classic() +
    xlab('DOY') +
    ylab('NDVI') +
    ylim(c(0.4,1))+
    geom_line(col=linecols[i]) +
    geom_ribbon(aes(ymin = ndvi - sderror,
                    ymax = ndvi + sderror), alpha = 0.2,fill=linecols[i])
  }else{
    p <- ggplot(NDVIdata, aes(x = doy, y = ndvi)) + #plot for concrete surface, different scaling
      ggtitle(paste(classnames[i],', n=',size,sep=''))+
      theme_classic() +
      xlab('DOY') +
      ylab('NDVI') +
      ylim(c(0,0.5))+
      geom_line(col=linecols[i]) +
      geom_ribbon(aes(ymin = ndvi - sderror,
                      ymax = ndvi + sderror), alpha = 0.2,fill=linecols[i])
  }
  ndviplotslist[[i]] <- p
}


#multiplot function, source: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  if (is.null(layout)) {

    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {

    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    

    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#display ribbon plots on one page
multiplot(plotlist = ndviplotslist, cols=3 )


