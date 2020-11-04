#author: Dominic Fawcett
#last modified: 12/10/2020
#description: script to compare NDVI orthomosaic values for direct and diffuse illumination conditions and leaf on and leaf off conditions
#NOTE: For use, adjust paths to datasets

library(raster)
library(rgdal)
library(RColorBrewer)
library(ggplot2)
library(stats)
library(reshape)

#import digitised crowns (version without concrete reference surface)
refShapes <- readOGR(dsn = 'D:/Microclimate/GPS', layer = "TreeClassPolysreprojExtNoVal") #adjust path to tree shapefiles

NDVIpath <- "D:/Microclimate/Sequoia/combNDVI/" #adjust path to NDVI orthomosaics
GREpath <- "D:/Microclimate/Sequoia/combGRE/" #adjust path to green-band orthomosaics

#read direct/diffuse NDVI and green band (GRE) image pairs

leafOffDif <- raster(paste0(NDVIpath,"NDVItrelus2019_03_22_dif.tif"))
leafOffDir <- raster(paste0(NDVIpath,"NDVItrelus2019_03_24_dir.tif"))  
leafOnDif <-  raster(paste0(NDVIpath,"NDVItrelus2019_07_22_dif.tif")) 
leafOnDir <-  raster(paste0(NDVIpath,"NDVItrelus2019_07_23_dir.tif")) 
leafOffDir2 <- raster(paste0(NDVIpath,"NDVItrelus2019_03_31_dir.tif")) 
leafOffDif2 <- raster(paste0(NDVIpath,"NDVItrelus2019_04_01_dif.tif")) 

leafOffDifGRE <- raster(paste0(GREpath,"GREtrelus2019_03_22_dif.tif")) 
leafOffDirGRE <- raster(paste0(GREpath,"GREtrelus2019_03_24_dir.tif"))   
leafOnDifGRE <- raster(paste0(GREpath,"GREtrelus2019_07_22_dif.tif")) 
leafOnDirGRE <- raster(paste0(GREpath,"GREtrelus2019_07_23_dir.tif")) 
leafOffDirGRE2 <- raster(paste0(GREpath,"GREtrelus2019_03_31_dir.tif")) 
leafOffDifGRE2 <- raster(paste0(GREpath,"GREtrelus2019_04_01_dif.tif")) 


NDVIstack <- stack(leafOffDif,resample(leafOffDir,leafOffDif),resample(leafOnDif,leafOffDif),resample(leafOnDir,leafOffDif),
                   resample(leafOffDif2,leafOffDif),resample(leafOffDir2,leafOffDif))#,resample(leafOnDif2,leafOffDif),resample(leafOnDir2,leafOffDif))
names(NDVIstack) <- c('leafOffDif','leafOffDir','leafOnDif','leafOnDir','leafOffDif2','leafOffDir2')#,'leafOnDif2','leafOnDir2')

GREstack <- stack(leafOffDifGRE,resample(leafOffDirGRE,leafOffDifGRE),resample(leafOnDifGRE,leafOffDifGRE),resample(leafOnDirGRE,leafOffDifGRE),
                  resample(leafOffDifGRE2,leafOffDifGRE),resample(leafOffDirGRE2,leafOffDifGRE))#,resample(leafOnDifGRE2,leafOffDifGRE),resample(leafOnDirGRE2,leafOffDifGRE))
names(GREstack) <- c('leafOffDifGRE','leafOffDirGRE','leafOnDifGRE','leafOnDirGRE','leafOffDifGRE2','leafOffDirGRE2')#,'leafOnDifGRE2','leafOnDirGRE2')
refShapesReproj <- spTransform(refShapes,crs(leafOffDif))

#extract NDVI and GRE values per tree
treesNDVI <- extract(NDVIstack,refShapesReproj, df=TRUE)
treesGRE <- extract(GREstack,refShapesReproj, df=TRUE)

#get the 80th percentile for masking
treesGRE80thpercentile <- extract(GREstack,refShapesReproj,fun=function(X,na.rm){quantile(X,0.8)}, df=TRUE)

#prepare dataframes for density plots
treesNDVIdf <- melt(treesNDVI,id=c('ID'))
names(treesNDVIdf) <- c('ID','illumination','NDVI')

treesNDVIdf$illumination <- gsub('leafOffDif','leaf-off dif.',treesNDVIdf$illumination)
treesNDVIdf$illumination <- gsub('leafOffDir','leaf-off dir.',treesNDVIdf$illumination)
treesNDVIdf$illumination <- gsub('leafOnDif','leaf-on dif.',treesNDVIdf$illumination)
treesNDVIdf$illumination<- gsub('leafOnDir','leaf-on dir.',treesNDVIdf$illumination)

#create density plots
treesNDVI4plot <- treesNDVI[,c(1,2,3)]
treesNDVIdf <- melt(treesNDVI4plot,id=c('ID'))
names(treesNDVIdf) <- c('ID','illumination','NDVI')

treesNDVIdf$illumination <- gsub('leafOffDif','diffuse',treesNDVIdf$illumination)
treesNDVIdf$illumination <- gsub('leafOffDir','direct',treesNDVIdf$illumination)

gg1 <- ggplot(treesNDVIdf) + geom_density(aes(x = NDVI,fill=illumination,size=illumination,linetype=illumination,color=illumination), alpha = 0.2)+
  scale_fill_manual(values = c("blue","red"))+
  scale_color_manual(values = c( "blue","red"))+
  scale_linetype_manual(values=c(2, 1))+
  scale_size_manual(values=c(0.5, 0.5))+
  xlim(0.25, 1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(subtitle = "a)")

treesNDVI4plot <- treesNDVI[,c(1,6,7)]
treesNDVIdf <- melt(treesNDVI4plot,id=c('ID'))
names(treesNDVIdf) <- c('ID','illumination','NDVI')

treesNDVIdf$illumination <- gsub('leafOffDif2','diffuse',treesNDVIdf$illumination)
treesNDVIdf$illumination <- gsub('leafOffDir2','direct',treesNDVIdf$illumination)

gg2 <- ggplot(treesNDVIdf) + geom_density(aes(x = NDVI,fill=illumination,size=illumination,linetype=illumination,color=illumination), alpha = 0.2)+
  scale_fill_manual(values = c("blue","red"))+
  scale_color_manual(values = c( "blue","red"))+
  scale_linetype_manual(values=c(2, 1))+
  scale_size_manual(values=c(0.5, 0.5))+
  xlim(0.25, 1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(subtitle = "b)")

treesNDVI4plot <- treesNDVI[,c(1,4,5)]
treesNDVIdf <- melt(treesNDVI4plot,id=c('ID'))
names(treesNDVIdf) <- c('ID','illumination','NDVI')

treesNDVIdf$illumination <- gsub('leafOnDif','diffuse',treesNDVIdf$illumination)
treesNDVIdf$illumination <- gsub('leafOnDir','direct',treesNDVIdf$illumination)

gg3 <- ggplot(treesNDVIdf) + geom_density(aes(x = NDVI,fill=illumination,size=illumination,linetype=illumination,color=illumination), alpha = 0.2)+
  scale_fill_manual(values = c("blue","red"))+
  scale_color_manual(values = c( "blue","red"))+
  scale_linetype_manual(values=c(2, 1))+
  scale_size_manual(values=c(0.5, 0.5))+
  xlim(0.25, 1)+
  #scale_color_manual(values = c("darkgreen", "darkred","green","red"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(subtitle = "c)")

grobplot <- arrangeGrob(gg1,gg2,gg3,nrow=3,ncol=1)

#modify plot output directory
ggsave("OUTDIR",grobplot)

#get NDVI metrics per tree for statistics

treeNDVIagVals <- matrix(nrow=nrow(treesGRE80thpercentile),ncol=7)
treeNDVIVals <- matrix(nrow=nrow(treesGRE80thpercentile),ncol=7)
treeNDVImedianVals <- matrix(nrow=nrow(treesGRE80thpercentile),ncol=7)

for (i in 1:nrow(treesGRE80thpercentile)){
  
  treeNDVIagVals[i,1] <- i
  treeNDVIVals[i,1] <- i
  for(j in 1:6){
  GREmask <- treesGRE[treesGRE[,1]==i,j+1]>treesGRE80thpercentile[i,j+1]
  NDVIvals <- treesNDVI[treesNDVI[,1]==i,j+1]
  NDVIgtGRE80pc <- NDVIvals[GREmask]
  treeNDVIagVals[i,j+1] <- mean(NDVIgtGRE80pc,na.rm=TRUE)
  treeNDVIVals[i,j+1] <- mean(NDVIvals,na.rm=TRUE)  
  treeNDVImedianVals[i,j+1] <- median(NDVIvals,na.rm=TRUE)  
  }
}

treeNDVIagVals <- data.frame(treeNDVIagVals)
names(treeNDVIagVals) <- c('class','leafOffDif','leafOffDir','leafOnDif','leafOnDir','leafOffDif2','leafOffDir2')#,'leafOnDif2','leafOnDir2')

treeNDVIVals <- data.frame(treeNDVIVals)
names(treeNDVIVals) <-  c('class','leafOffDif','leafOffDir','leafOnDif','leafOnDir','leafOffDif2','leafOffDir2')#,'leafOnDif2','leafOnDir2')#c('class','leaf-off dif.','leaf-off dir.','leaf-on dif.','leaf-on dir.')

treeNDVImedianVals <- data.frame(treeNDVImedianVals)
names(treeNDVImedianVals) <-  c('class','leafOffDif','leafOffDir','leafOnDif','leafOnDir','leafOffDif2','leafOffDir2')#,'leafOnDif2','leafOnDir2')#c('class','leaf-off dif.','leaf-off dir.','leaf-on dif.','leaf-on dir.')


#extract statistics

#leafOffcomp 80perc

MADOff <- mean(abs(treeNDVIagVals$leafOffDir-treeNDVIagVals$leafOffDif))
percMADOff <- MADOff/mean(treeNDVIagVals$leafOffDir)*100
biasoff <- mean(treeNDVIagVals$leafOffDir-treeNDVIagVals$leafOffDif)
ttestOff <- t.test(treeNDVIagVals$leafOffDir,treeNDVIagVals$leafOffDif)
ttestpvalOff <- ttestOff$p.value
offModel <- lm(treeNDVIagVals$leafOffDir~treeNDVIagVals$leafOffDif)
rsqOff <- summary(offModel)$adj.r.squared

#leafOffcomp mean

MADOff <- mean(abs(treeNDVIVals$leafOffDir-treeNDVIVals$leafOffDif))
percMADOff <- MADOff/mean(treeNDVIVals$leafOffDir)*100
biasoff <- mean(treeNDVIVals$leafOffDir-treeNDVIVals$leafOffDif)
ttestOff <- t.test(treeNDVIVals$leafOffDir,treeNDVIVals$leafOffDif)
ttestpvalOff <- ttestOff$p.value
offModel <- lm(treeNDVIVals$leafOffDir~treeNDVIVals$leafOffDif)
rsqOff <- summary(offModel)$adj.r.squared

#leafOffcomp P2 80perc

MADOff <- mean(abs(treeNDVIagVals$leafOffDir2-treeNDVIagVals$leafOffDif2))
percMADOff <- MADOff/mean(treeNDVIagVals$leafOffDir2)*100
biasoff <- mean(treeNDVIagVals$leafOffDir2-treeNDVIagVals$leafOffDif2)
ttestOff <- t.test(treeNDVIagVals$leafOffDir2,treeNDVIagVals$leafOffDif2)
ttestpvalOff <- ttestOff$p.value
offModel <- lm(treeNDVIagVals$leafOffDir2~treeNDVIagVals$leafOffDif2)
rsqOff <- summary(offModel)$adj.r.squared

#leafOffcomp P2 mean

MADOff <- mean(abs(treeNDVIVals$leafOffDir2-treeNDVIVals$leafOffDif2))
percMADOff <- MADOff/mean(treeNDVIVals$leafOffDir2)*100
biasoff <- mean(treeNDVIVals$leafOffDir2-treeNDVIVals$leafOffDif2)
ttestOff <- t.test(treeNDVIVals$leafOffDir2,treeNDVIVals$leafOffDif2)
ttestpvalOff <- ttestOff$p.value
offModel <- lm(treeNDVIVals$leafOffDir2~treeNDVIVals$leafOffDif2)
rsqOff <- summary(offModel)$adj.r.squared

#leafOffcomp P2 median

MADOff <- mean(abs(treeNDVImedianVals$leafOffDir2-treeNDVImedianVals$leafOffDif2))
percMADOff <- MADOff/mean(treeNDVImedianVals$leafOffDir2)*100
biasoff <- mean(treeNDVImedianVals$leafOffDir2-treeNDVImedianVals$leafOffDif2)
ttestOff <- t.test(treeNDVImedianVals$leafOffDir2,treeNDVImedianVals$leafOffDif2)
ttestpvalOff <- ttestOff$p.value
offModel <- lm(treeNDVImedianVals$leafOffDir2~treeNDVImedianVals$leafOffDif2)
rsqOff <- summary(offModel)$adj.r.squared


#analyse difference between pre green-up NDVI mean and fully greened up

treeclasses <- refShapes@data$class
treeclassnrs <- c(seq(0:17))

#without 80 percentile
meanDifOff <-mean(treeNDVIVals$leafOffDif[treeclasses %in% c(0,1,2,3,4,6,9,10,11,12,13,14,15,16)],na.rm=TRUE) 
meanDifOn <-mean(treeNDVIVals$leafOnDif[treeclasses %in% c(0,1,2,3,4,6,9,10,11,12,13,14,15,16)],na.rm=TRUE) 
onOffDifference <- meanDifOn-meanDifOff

#with 80th percentile
meanDifOff <-mean(treeNDVIagVals$leafOffDif[treeclasses %in% c(0,1,2,3,4,6,9,10,11,12,13,14,15,16)],na.rm=TRUE) 
meanDifOn <-mean(treeNDVIagVals$leafOnDif[treeclasses %in% c(0,1,2,3,4,6,9,10,11,12,13,14,15,16)],na.rm=TRUE) 
onOffDifference <- meanDifOn-meanDifOff




