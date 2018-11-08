################################################################
## R code for **Browning et al., 2017. Breaks in MODIS time series 
## portend vegetation change – verification using long-term data 
## in an arid grassland ecosystem**
## Coding and commenting by Jonathan Maynard (USDA-ARS)
## Date: 03-13-2017
## Questions? jonathan.maynard@ars.usda.gov
###############################################################

###############################################################
## Required libraries
###############################################################
library(MODIS)
library(maptools)
library(rgdal)
library(raster)
library(maps)
library(doMPI)
library(foreach)
library(bfast)
library(strucchange)
library(Cairo)

#set up parrallel backend
cl <- startMPIcluster(count=24)
registerDoMPI(cl)
getDoParWorkers()
###############################################################
#Set repository for MODIS imagery
###############################################################
MODIS:::checkTools()
MODISoptions(localArcPath="~/MODIS/MODIS_ARC/",
             outDirPath="~/MODIS/MODIS_ARC/PROCESSES")
 

# #Vegetation Indicies (250 M)
#MODIS tile for USDA Jornada Experimental Range
getHdf(product="MOD13Q1",begin="2015321", tileH=9,tileV=5, collection='005') 
 

## Read in the Jornada boundary shapefile (specify source directory)
JRN <- readOGR(dsn="~/JRN_boundary",layer="JRN_boundary")
JRN.point <- readOGR(dsn="~/Jornada_1M_plot_centroid",layer="Jornada_1m_plot_centroid")



## Set up a raster "template" to use in rasterize()
##expand extent 20000 m on every side to prevent clipping across study area
ext <-  extent (JRN)+40000
xy <- abs(apply(as.matrix(bbox(ext)), 1, diff))  #this calculated the absolute distance in the x and y directions
r <- raster(ext, ncol=xy[1]/250, nrow=xy[2]/250) #this creates a raster box of the extended study area

## Rasterize the shapefile
JRN_raster <-rasterize(JRN, r) #we then create a raster of the JRN boundary within the extended study area
res(JRN_raster) <- 250 #assign resolution
plot(JRN_raster)

#find out which SDS numbers correspond to which bands by examining one downloaded image
#MOD13Q1
getSds(HdfName="MOD13Q1.A2000049.h08v04.005.2006269063040.hdf")

# $SDSnames
# [1] "250m 16 days NDVI"                       [2] "250m 16 days EVI"                      
# [3] "250m 16 days VI Quality"                 [4] "250m 16 days red reflectance"          
# [5] "250m 16 days NIR reflectance"            [6] "250m 16 days blue reflectance"         
# [7] "250m 16 days MIR reflectance"            [8] "250m 16 days view zenith angle"        
# [9] "250m 16 days sun zenith angle"           [10] "250m 16 days relative azimuth angle"   
# [11] "250m 16 days composite day of the year" [12]"250m 16 days pixel reliability" 

#run MODIS function runGdal to extract MODIS Data
#MOD13Q1
runGdal(job="JRN_Modis_VI", product="MOD13Q1", collection='005', extent=JRN_raster, begin="2000001", SDSstring="111111")

#set the path to get the correct data on the source directory: 
#path to 250M VI's
path <- paste(options("MODIS_outDirPath"),"/JRN_Modis_VI",sep="") 
ndvi <- preStack(path=path, pattern="*_NDVI.tif$")
JRN_NDVI_brick<- brick(stack(ndvi),  filename='/output_directory/JRN_NDVI_brick', overwrite=TRUE)


#Raster layers need to be in a vector list to parallelize accross a cluster.
stack2list <- function(stack){
  list<-vector("list", nlayers(stack))
  for(i in 1:nlayers(stack)) {
  list[[i]] <- stack[[i]]
  }
  return(list)
}

JRN_NDVI_list<- stack2list(stack(ndvi))
####################################################################################################################################

######################################################################################################
#MODIS QA Processing

## QA for MOD13 products
## MODIS13Q1 has 16 bit data so we want the fist 16 integers.
## MODIS QA data is parsed from right to left, and the individual bit groups (bit words) are read from left to right.
## Example
## as.integer(intToBits(7425)[16:1])
## produces 0 0 0 1 1 1 0 1 0 0 0 0 0 0 0 1 , the bit words or combinations of integers going from right to left mean different things
## for different MODIS products

###############MODIS13Q1 QA lookup table

QC_Data <- data.frame(Integer_Value = 0:65535,
                      Bit15 = NA,
                      Bit14 = NA,
                      Bit13 = NA,
                      Bit12 = NA,
                      Bit11 = NA,
                      Bit10 = NA,
                      Bit9 = NA,
                      Bit8 = NA,
                      Bit7 = NA,
                      Bit6 = NA,
                      Bit5 = NA,
                      Bit4 = NA,
                      Bit3 = NA,
                      Bit2 = NA,
                      Bit1 = NA,
                      Bit0 = NA,
                      QA_word1 = NA,
                      QA_word2 = NA,
                      QA_word3 = NA,
                      QA_word4 = NA,
                      QA_word5 = NA,
                      QA_word6 = NA,
                      QA_word7 = NA,
                      QA_word8 = NA,
                      QA_word9 = NA
                      
)

for(i in QC_Data$Integer_Value){
  AsInt <- as.integer(intToBits(i)[1:16])
  QC_Data[i+1,2:17]<- AsInt[16:1]
}

#VI Quality
QC_Data$QA_word1[QC_Data$Bit1 == 0 & QC_Data$Bit0==0] <- "VI GOOD"
QC_Data$QA_word1[QC_Data$Bit1 == 0 & QC_Data$Bit0==1] <- "VI Produced,Other Quality"
QC_Data$QA_word1[QC_Data$Bit1 == 1 & QC_Data$Bit0==0] <- "No Pixel,clouds"
QC_Data$QA_word1[QC_Data$Bit1 == 1 & QC_Data$Bit0==1] <- "No Pixel, Other QA"

#VI Usefullness
QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==0 & QC_Data$Bit3 == 0 & QC_Data$Bit2==0] <- "Highest Quality"
QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==0 & QC_Data$Bit3 == 0 & QC_Data$Bit2==1] <- "High Quality"
QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==0 & QC_Data$Bit3 == 1 & QC_Data$Bit2==0] <- "Good Quality"
QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==0 & QC_Data$Bit3 == 1 & QC_Data$Bit2==1] <- "Acceptable Quality"
QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==1 & QC_Data$Bit3 == 0 & QC_Data$Bit2==0] <- "Fair Quality"
QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==1 & QC_Data$Bit3 == 0 & QC_Data$Bit2==1] <- "Intermediate Quality"
QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==1 & QC_Data$Bit3 == 1 & QC_Data$Bit2==0] <- "Below Intermediate Quality"
QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==1 & QC_Data$Bit3 == 1 & QC_Data$Bit2==1] <- "Average Quality"
QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==0 & QC_Data$Bit3 == 0 & QC_Data$Bit2==0] <- "Below Average Quality"
QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==0 & QC_Data$Bit3 == 0 & QC_Data$Bit2==1] <- "Questionable Quality"
QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==0 & QC_Data$Bit3 == 1 & QC_Data$Bit2==0] <- "Above Marginal Quality"
QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==0 & QC_Data$Bit3 == 1 & QC_Data$Bit2==1] <- "Marginal Quality"
QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==1 & QC_Data$Bit3 == 0 & QC_Data$Bit2==0] <- "Low Quality"
QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==1 & QC_Data$Bit3 == 0 & QC_Data$Bit2==1] <- "No Atmopheric Correction"
QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==1 & QC_Data$Bit3 == 1 & QC_Data$Bit2==0] <- "Quality too low to be useful"
QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==1 & QC_Data$Bit3 == 1 & QC_Data$Bit2==1] <- "Not useful for other reasons"

#Aerosol quantity
QC_Data$QA_word3[QC_Data$Bit7 == 0 & QC_Data$Bit6==0] <- "Climatology"
QC_Data$QA_word3[QC_Data$Bit7 == 0 & QC_Data$Bit6==1] <- "Low"
QC_Data$QA_word3[QC_Data$Bit7 == 1 & QC_Data$Bit6==0] <- "Intermediate"
QC_Data$QA_word3[QC_Data$Bit7 == 1 & QC_Data$Bit6==1] <- "High"

#Adjacent Clouds
QC_Data$QA_word4[QC_Data$Bit8 == 0] <- "No"
QC_Data$QA_word4[QC_Data$Bit8 == 1] <- "Yes"

#Atomospheric BRDF correction
QC_Data$QA_word5[QC_Data$Bit9 == 0] <- "No"
QC_Data$QA_word5[QC_Data$Bit9 == 1] <- "Yes"

#Mixed Clouds
QC_Data$QA_word6[QC_Data$Bit10 == 0] <- "No"
QC_Data$QA_word6[QC_Data$Bit10 == 1] <- "Yes"

#Land/Water Mask
QC_Data$QA_word7[QC_Data$Bit13==0 & QC_Data$Bit12 == 0 & QC_Data$Bit11==0] <- "Shallow ocean"
QC_Data$QA_word7[QC_Data$Bit13==0 & QC_Data$Bit12 == 0 & QC_Data$Bit11==1] <- "Land"
QC_Data$QA_word7[QC_Data$Bit13==0 & QC_Data$Bit12 == 1 & QC_Data$Bit11==0] <- "Coastline/Shoreline"
QC_Data$QA_word7[QC_Data$Bit13==0 & QC_Data$Bit12 == 1 & QC_Data$Bit11==1] <- "Shallow inland water"
QC_Data$QA_word7[QC_Data$Bit13==1 & QC_Data$Bit12 == 0 & QC_Data$Bit11==0] <- "Ephemeral water"
QC_Data$QA_word7[QC_Data$Bit13==1 & QC_Data$Bit12 == 0 & QC_Data$Bit11==1] <- "Deep inland water"
QC_Data$QA_word7[QC_Data$Bit13==1 & QC_Data$Bit12 == 1 & QC_Data$Bit11==0] <- "Moderate or continental ocean"
QC_Data$QA_word7[QC_Data$Bit13==1 & QC_Data$Bit12 == 1 & QC_Data$Bit11==1] <- "Deep ocean"

#Possible snow/ice
QC_Data$QA_word8[QC_Data$Bit14 == 0] <- "No"
QC_Data$QA_word8[QC_Data$Bit14 == 1] <- "Yes"

#Possible shadow
QC_Data$QA_word9[QC_Data$Bit15 == 0] <- "No"
QC_Data$QA_word9[QC_Data$Bit15 == 1] <- "Yes"
##############################################################################################
#Run QA on MODIS 250m NDVI
#extract unique QA values from raster brick
JRN_QA_unique<- foreach(i=1:length(JRN_QA_list), .combine='c', .packages = c("raster", "rgdal")) %dopar% {   
  unique(JRN_QA_list[[i]])
}

#final list of unique QA integers from image brick
JRN_QA_unique<-unique(JRN_QA_unique)

#now subset QA_Data to extract only relevant QA integers
JRN_QC_Data <- subset(QC_Data, Integer_Value %in% JRN_QA_unique)

##Now subset based on the most important QA variables
#select all with VI quality of 00
VI.good <- subset(JRN_QC_Data, QA_word1=="VI GOOD")

#select subset of VI Usefullness 
VI.okay <- subset(subset(JRN_QC_Data, QA_word1=="VI Produced,Other Quality"), QA_word2=="Highest Quality" | QA_word2=="High Quality" | QA_word2=="Good Quality"
                  | QA_word2=="Fair Quality" | QA_word2=="Intermediate Quality" | QA_word2=="Below Intermediate Quality"| QA_word2=="Average Quality")

#combine good and okay integer values and then assign to a vector.
VI.QA <- rbind(VI.good, VI.okay)
QA.int<-VI.QA$Integer_Value

#Create JRN_QA raster list object to modify
JRN_QA <- JRN_QA_list

#First loop assigns all good pixels a value of 1 based on previous selection
for(i in 1:length(JRN_QA)){
  
  for(j in 1:length(QA.int)){
    JRN_QA[[i]][JRN_QA[[i]]==QA.int[j]] <-1
    
  }}

#Second loop assigns all values greater than 1 to NA
for(i in 1:length(JRN_QA)){
  JRN_QA[[i]][JRN_QA[[i]]>1] <-NA
}

## Final loop creates mask to mask out bad pixels
## This loop uses the foreach fucntion in a cluster configuration and thus requires 
## that the raster data be in a list format so that it can subset the data

#################Mask NDVI
JRN_NDVI_mask<-JRN_NDVI_list
JRN_NDVI_mask<- foreach(i=1:length(JRN_NDVI_mask), .packages = c("raster", "rgdal"))  %dopar% { 
  JRN_NDVI_mask[[i]]<-mask(JRN_NDVI_mask[[i]], (JRN_QA[[i]]))
}

#Final NDVI product
JRN_NDVI_mask <- brick(stack(JRN_NDVI_mask),  filename='/output_directory/JRN_NDVI_mask', overwrite=TRUE)

#Test it to see results
plot(JRN_NDVI_mask[[1]])
medianNDVI.mask <- calc(JRN_NDVI_mask, fun=function(x) median(x, na.rm = TRUE)) 
plot(medianNDVI.mask)
###################################################################################################
########################################################
####select sites based on following names: BASN CALI COLL EAST GRAV IBPE NORT RABB SAND SMAL SUMM TAYL TOBO WELL WEST
#NPP Sites
BASN <- JRN.point[JRN.point@data$SITE %in% c("BASN"), ]
CALI <- JRN.point[JRN.point@data$SITE %in% c("CALI"), ]
COLL <- JRN.point[JRN.point@data$SITE %in% c("COLL"), ]
EAST <- JRN.point[JRN.point@data$SITE %in% c("EAST"), ]
GRAV <- JRN.point[JRN.point@data$SITE %in% c("GRAV"), ]
IBPE <- JRN.point[JRN.point@data$SITE %in% c("IBPE"), ]
NORT <- JRN.point[JRN.point@data$SITE %in% c("NORT"), ]
RABB <- JRN.point[JRN.point@data$SITE %in% c("RABB"), ]
SAND <- JRN.point[JRN.point@data$SITE %in% c("SAND"), ]
SMAL <- JRN.point[JRN.point@data$SITE %in% c("SMAL"), ]
SUMM <- JRN.point[JRN.point@data$SITE %in% c("SUMM"), ]
TAYL <- JRN.point[JRN.point@data$SITE %in% c("TAYL"), ]
TOBO <- JRN.point[JRN.point@data$SITE %in% c("TOBO"), ]
WELL <- JRN.point[JRN.point@data$SITE %in% c("WELL"), ]
WEST <- JRN.point[JRN.point@data$SITE %in% c("WEST"), ]

#######extract plot points from raster brick: BASN CALI COLL EAST GRAV IBPE NORT RABB SAND SMAL SUMM TAYL TOBO WELL WEST
###############NDVI_mask
BASN.ndvi.points<-extract(JRN_NDVI_mask, BASN) #49 points from BASN plot
CALI.ndvi.points<-extract(JRN_NDVI_mask, CALI)
COLL.ndvi.points<-extract(JRN_NDVI_mask, COLL)
EAST.ndvi.points<-extract(JRN_NDVI_mask, EAST)
GRAV.ndvi.points<-extract(JRN_NDVI_mask, GRAV)
IBPE.ndvi.points<-extract(JRN_NDVI_mask, IBPE)
NORT.ndvi.points<-extract(JRN_NDVI_mask, NORT)
RABB.ndvi.points<-extract(JRN_NDVI_mask, RABB)
SAND.ndvi.points<-extract(JRN_NDVI_mask, SAND)
SMAL.ndvi.points<-extract(JRN_NDVI_mask, SMAL)
SUMM.ndvi.points<-extract(JRN_NDVI_mask, SUMM)
TAYL.ndvi.points<-extract(JRN_NDVI_mask, TAYL)
TOBO.ndvi.points<-extract(JRN_NDVI_mask, TOBO)
WELL.ndvi.points<-extract(JRN_NDVI_mask, WELL)
WEST.ndvi.points<-extract(JRN_NDVI_mask, WEST)

#mean across 49 points
BASN.ndvi.m <- colMeans(BASN.ndvi.points)
CALI.ndvi.m <- colMeans(CALI.ndvi.points)
COLL.ndvi.m <- colMeans(COLL.ndvi.points)
EAST.ndvi.m <- colMeans(EAST.ndvi.points)
GRAV.ndvi.m <- colMeans(GRAV.ndvi.points)
IBPE.ndvi.m <- colMeans(IBPE.ndvi.points)
NORT.ndvi.m <- colMeans(NORT.ndvi.points)
RABB.ndvi.m <- colMeans(RABB.ndvi.points)
SAND.ndvi.m <- colMeans(SAND.ndvi.points)
SMAL.ndvi.m <- colMeans(SMAL.ndvi.points)
SUMM.ndvi.m <- colMeans(SUMM.ndvi.points)
TAYL.ndvi.m <- colMeans(TAYL.ndvi.points)
TOBO.ndvi.m <- colMeans(TOBO.ndvi.points)
WELL.ndvi.m <- colMeans(WELL.ndvi.points)
WEST.ndvi.m <- colMeans(WEST.ndvi.points)

############################################################################################################
########## The above code (lns 24-310) automates the downloading and processing of MODIS imagery for the USDA-ARS 
########## Jornada experimental range and the extraction of spatially weighted NDVI data for each NPP site. This
########## paper has include a .csv file with the processed NDIV data for each NPP site on a github repository.

#Load in JRN_NPP_NDVI.csv file from github
JRN.NPP.NDVI <- read.csv("JRN_NPP_NDVI.csv")
BASN.ndvi.m <- JRN.NPP.NDVI$BASN.ndvi.m
CALI.ndvi.m <- JRN.NPP.NDVI$CALI.ndvi.m
COLL.ndvi.m <- JRN.NPP.NDVI$COLL.ndvi.m
EAST.ndvi.m <- JRN.NPP.NDVI$EAST.ndvi.m
GRAV.ndvi.m <- JRN.NPP.NDVI$GRAV.ndvi.m
IBPE.ndvi.m <- JRN.NPP.NDVI$IBPE.ndvi.m
NORT.ndvi.m <- JRN.NPP.NDVI$NORT.ndvi.m
RABB.ndvi.m <- JRN.NPP.NDVI$RABB.ndvi.m
SAND.ndvi.m <- JRN.NPP.NDVI$SAND.ndvi.m
SMAL.ndvi.m <- JRN.NPP.NDVI$SMAL.ndvi.m
SUMM.ndvi.m <- JRN.NPP.NDVI$SUMM.ndvi.m
TAYL.ndvi.m <- JRN.NPP.NDVI$TAYL.ndvi.m
TOBO.ndvi.m <- JRN.NPP.NDVI$TOBO.ndvi.m
WELL.ndvi.m <- JRN.NPP.NDVI$WELL.ndvi.m
WEST.ndvi.m <- JRN.NPP.NDVI$WEST.ndvi.m
############################################################################################
######################################################
####Now we process the MODIS data with the bfast package
## Processing raster bricks (satellite image time series of 16-day NDVI images) 

#modified code for function bfastts to account for leap year
bfastts.na <- function (data, dates, type = c("irregular", "16-day", "quarterly")) 
{
  yday365 <- function(x) {
    x <- as.POSIXlt(x)
    mdays <- c(31L, 28L, 31L, 30L, 31L, 30L, 31L, 31L, 30L, 
               31L, 30L, 31L)
    cumsum(c(0L, mdays))[1L + x$mon] + x$mday
  }
  if (type == "irregular") {
    zz <- zoo(data, 1900 + as.POSIXlt(dates)$year + (yday365(dates) - 
                                                       1)/365, frequency = 365)
  }
  if (type == "16-day") {
    z <- zoo(data, dates)
    yr <- as.numeric(format(time(z), "%Y"))
    jul <- as.numeric(format(time(z), "%j"))
    delta <- min(unlist(tapply(jul, yr, diff)))
    zz <- aggregate(z, yr + (jul - 1)/delta/23)
  }
  
  if (type == "quarterly") {
    zz <- zoo(data, dates)
  }
  #tso <- na.approx(as.ts(zz))
  tso <- na.approx(ts(coredata(zz),frequency=(365.25/16),start=start(zz))) ####modified part of bfastts code to account for leap year
  return(tso)
}

##############################################################################################################################
#Old code for bfast function, verions 1.4.4. This version uses the 'strucchange' package for breakpoint
# analysis whereas newer versions, e.g., version 1.5.7, do not and produce a different result.
bfast_old <- function(Yt, h=0.15, season =c("dummy","harmonic","none"), max.iter = NULL, breaks = NULL, hpc = "none")
{
  season <- match.arg(season)
  ti <- time(Yt)
  f <- frequency(Yt)      # on cycle every f time points (seasonal cycle)
  if(class(Yt)!="ts")
    stop ("Not a time series object")
  ## return value
  output <- list()
  
  # Start the iterative procedure and for first iteration St=decompose result
  St <- stl(Yt, "periodic")$time.series[, "seasonal"]
  Tt <- 0
  
  # seasonal model setup
  if (season=="harmonic") {
    w <- 1/f # f = 23 when freq=23 :-)
    tl <- 1:length(Yt)
    co <- cos(2*pi*tl*w); si <- sin(2*pi*tl*w)
    co2 <- cos(2*pi*tl*w*2);si2 <- sin(2*pi*tl*w*2)
    co3 <- cos(2*pi*tl*w*3);si3 <- sin(2*pi*tl*w*3)
    smod <- Wt ~ co+si+co2+si2+co3+si3
  } else if (season=="dummy") {
    D <- seasonaldummy(Yt)
    D[rowSums(D)==0,] <- -1
    smod <- Wt ~ -1+D
  } else if (season=="none") {
    #print("No seasonal model will be fitted!")
  } else stop("Not a correct seasonal model is selected ('harmonic' or 'dummy') ")
  
  # number/timing of structural breaks in the trend/seasonal component
  Vt.bp <- 0
  Wt.bp <- 0 
  CheckTimeTt <- 1
  CheckTimeSt <- 1
  
  i <- 0
  while ( (!identical(CheckTimeTt,Vt.bp) | !identical(CheckTimeSt,Wt.bp)) & i < max.iter)
  {
    CheckTimeTt <- Vt.bp
    CheckTimeSt <- Wt.bp
    # TREND
    Vt <- Yt-St
    #         p.Vt <- sctest(efp(Vt ~ ti, h=h, type= "OLS-MOSUM"))
    #         if (p.Vt$p.value <=0.05) 
    #         {
    bp.Vt <- breakpoints(Vt ~ ti, h=h,breaks=breaks, hpc = hpc)
    nobp.Vt <- is.na(breakpoints (bp.Vt)[1])
    #         } 
    #         else 
    #         {
    #           nobp.Vt <- TRUE
    #           bp.Vt <- NA       
    #         }
    if (nobp.Vt)
    {
      fm0 <- lm(Vt ~  ti)
      Vt.bp <- 0      # no breaks times
      Tt <- ts(fitted(fm0))     # Data minus trend
      tsp(Tt) <- tsp(Yt)
      ci.Vt <- NA
    } 
    else
    {
      fm1 <- lm(Vt ~ breakfactor(bp.Vt)/ti)
      ci.Vt <- confint(bp.Vt, het.err = FALSE)
      Vt.bp <- ci.Vt$confint[,2]
      Tt <- ts(fitted(fm1))     # Data minus trend
      tsp(Tt) <- tsp(Yt)
    }
    
    # SEASONAL COMPONENT
    if (season=="none") {
      Wt <- 0
      St <- 0
      bp.Wt <- NA; ci.Wt <- NA; nobp.Wt<- TRUE
    } else
    {
      Wt <- Yt-Tt
      #        p.Wt <- sctest(efp(smod, h=h, type= "OLS-MOSUM"))      # preliminary test 
      #        if (p.Wt$p.value <=0.05) # OR statement 
      #        {
      bp.Wt <- breakpoints(smod, h=h,breaks=breaks, hpc = hpc) # Breakpoints in the seasonal component
      nobp.Wt <- is.na(breakpoints (bp.Wt)[1])
      #        } 
      #        else 
      #        {
      #            nobp.Wt <- TRUE
      #            bp.Wt <- NA       
      #        }
      if (nobp.Wt)
      {
        sm0 <- lm(smod)
        St <- ts(fitted(sm0))  #  The fitted seasonal component
        tsp(St) <- tsp(Yt)
        Wt.bp <- 0             # no seasonal breaks
        ci.Wt <- NA
      } 
      else
      {
        if(season=="dummy") sm1 <-lm(Wt ~ -1+D %in% breakfactor(bp.Wt))
        if(season=="harmonic") sm1 <- lm(Wt ~ (co+si+co2+si2+co3+si3) %in% breakfactor(bp.Wt)) 
        St <- ts(fitted(sm1))  #  The fitted seasonal component
        tsp(St) <- tsp(Yt)
        ci.Wt <- confint(bp.Wt, het.err = FALSE)
        Wt.bp <- ci.Wt$confint[,2] 
      }
    }
    i <- i+1
    output[[i]] <- list(Tt=Tt,St=St,Nt=Yt-Tt-St,
                        Vt=Vt, bp.Vt=bp.Vt, Vt.bp=Vt.bp, ci.Vt=ci.Vt,
                        Wt=Wt, bp.Wt=bp.Wt, Wt.bp=Wt.bp, ci.Wt=ci.Wt)
  }
  if (!nobp.Vt) # probably only works well for dummy model!
  {
    Vt.nrbp <- length(bp.Vt$breakpoints)
    co <- coef(fm1) # final fitted trend model
    Mag <- matrix(NA,Vt.nrbp,3)
    for (r in 1:Vt.nrbp) 
    {
      if (r==1) 
        y1 <- co[1]+co[r+Vt.nrbp+1]*ti[Vt.bp[r]]
      else 
        y1 <- co[1]+co[r]+co[r+Vt.nrbp+1]*ti[Vt.bp[r]]
      y2 <- (co[1]+co[r+1])+co[r+Vt.nrbp+2]*ti[Vt.bp[r]+1]
      Mag[r,1] <- y1
      Mag[r,2] <- y2
      Mag[r,3] <- y2-y1   
    }
    index <- which.max(abs(Mag[,3]))
    m.x <- rep(Vt.bp[index],2)
    m.y <- c(Mag[index,1],Mag[index,2]) #Magnitude position
    Magnitude <- Mag[index,3] # Magnitude of biggest change
    Time <- Vt.bp[index]
  } 
  else 
  {
    m.x <- NA; m.y <- NA
    Magnitude <- 0  # if we do not detect a break then the magnitude is zero
    Time <- NA # if we do not detect a break then we have no timing of the break
    Mag <- 0
  }
  return(structure(list(Yt=Yt,output=output,nobp=list(Vt=nobp.Vt,Wt=nobp.Wt),Magnitude=Magnitude,Mags=Mag,
                        Time=Time,jump=list(x=ti[m.x],y=m.y)),class="bfast"))  
}

###############################################################################################################

### BASN
############################################################################################
## apply on MODIS composite point
ndvi.date <- (orgTime(ndvi,nDays=16,begin="2000049", end="2014177", pillow=0))$inputLayerDates  #, pillow=90  #original model run until end="2014177" but later cropped for figure
basn.ndvi <- bfastts.na(as.numeric(BASN.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(basn.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
basn.ndvi.fit<-stl(basn.ndvi, s.window=8)
colnames(basn.ndvi.fit$time.series)[2]<-"abrupt"
plot(basn.ndvi.fit)

basn.fit <- bfast_old(basn.ndvi, season="harmonic",h=0.15, max.iter=10)
plot(basn.fit) 

############################################################################################
### CALI
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
cali.ndvi <- bfastts.na(as.numeric(CALI.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(cali.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
cali.ndvi.fit<-stl(cali.ndvi, s.window=8)
colnames(cali.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
cali.fit <- bfast_old(cali.ndvi, season="harmonic", max.iter=10) 

############################################################################################
### COLL
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
coll.ndvi <- bfastts.na(as.numeric(COLL.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(coll.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
coll.ndvi.fit<-stl(coll.ndvi, s.window=8)
colnames(coll.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
coll.fit <- bfast_old(coll.ndvi, season="harmonic", max.iter=10) 

############################################################################################
### EAST
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
east.ndvi <- bfastts.na(as.numeric(EAST.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(east.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
east.ndvi.fit<-stl(east.ndvi, s.window=18)
colnames(east.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
east.fit <- bfast_old(east.ndvi, season="harmonic",h=.17, max.iter=10) 

############################################################################################
### GRAV
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
grav.ndvi <- bfastts.na(as.numeric(GRAV.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(grav.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
grav.ndvi.fit<-stl(grav.ndvi, s.window=8)
colnames(grav.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
grav.fit <- bfast_old(grav.ndvi, season="harmonic", max.iter=10) 

############################################################################################
### IBPE
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
ibpe.ndvi <- bfastts.na(as.numeric(IBPE.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(ibpe.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
ibpe.ndvi.fit<-stl(ibpe.ndvi, s.window=8)
colnames(ibpe.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
ibpe.fit <- bfast_old(ibpe.ndvi, season="harmonic", max.iter=10) 

############################################################################################
### NORT
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
nort.ndvi <- bfastts.na(as.numeric(NORT.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(nort.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
nort.ndvi.fit<-stl(nort.ndvi, s.window=8)
colnames(nort.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
nort.fit <- bfast_old(nort.ndvi, season="harmonic", max.iter=10) 

############################################################################################
### RABB
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
rabb.ndvi <- bfastts.na(as.numeric(RABB.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(rabb.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
rabb.ndvi.fit<-stl(rabb.ndvi, s.window=8)
colnames(rabb.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
rabb.fit <- bfast_old(rabb.ndvi, season="harmonic", max.iter=10) 

############################################################################################
### SAND
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
sand.ndvi <- bfastts.na(as.numeric(SAND.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(sand.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
sand.ndvi.fit<-stl(sand.ndvi, s.window=8)
colnames(sand.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
sand.fit <- bfast_old(sand.ndvi, season="harmonic", max.iter=10) 

############################################################################################
### SMAL
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
smal.ndvi <- bfastts.na(as.numeric(SMAL.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(smal.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
smal.ndvi.fit<-stl(smal.ndvi, s.window=8)
colnames(smal.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
smal.fit <- bfast_old(smal.ndvi, season="harmonic", max.iter=10) 

############################################################################################
### SUMM
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
summ.ndvi <- bfastts.na(as.numeric(SUMM.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(summ.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
summ.ndvi.fit<-stl(summ.ndvi, s.window=8)
colnames(summ.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
summ.fit <- bfast_old(summ.ndvi, season="harmonic", max.iter=10) 

############################################################################################
### TAYL
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
tayl.ndvi <- bfastts.na(as.numeric(TAYL.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(tayl.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
tayl.ndvi.fit<-stl(tayl.ndvi, s.window=8)
colnames(tayl.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
tayl.fit <- bfast_old(tayl.ndvi, season="harmonic", max.iter=10) 

############################################################################################
### TOBO
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
tobo.ndvi <- bfastts.na(as.numeric(TOBO.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(tobo.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
tobo.ndvi.fit<-stl(tobo.ndvi, s.window=8)
colnames(tobo.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
tobo.fit <- bfast_old(tobo.ndvi, season="harmonic", max.iter=10) 

############################################################################################
### WELL
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
well.ndvi <- bfastts.na(as.numeric(WELL.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(well.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
well.ndvi.fit<-stl(well.ndvi, s.window=8)
colnames(well.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
well.fit <- bfast_old(well.ndvi, season="harmonic", max.iter=10) 

############################################################################################
### WEST
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
west.ndvi <- bfastts.na(as.numeric(WEST.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(west.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
west.ndvi.fit<-stl(west.ndvi, s.window=8)
colnames(west.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
west.fit <- bfast_old(west.ndvi, season="harmonic", max.iter=10) 

####################################################################################################
#Function to extract Trend and Seasonal breaks
bfast.breaks <- function(bfast,date){
  
  trend <-list(list())
  if(is.na(bfast$output[[length(bfast$output)]]$bp.Vt$breakpoints)[1]==TRUE){
    trend <-NA
  } else {
    for(i in 1:length(bfast$output[[length(bfast$output)]]$bp.Vt$breakpoints)){
      trend[[i]] <- date[bfast$output[[length(bfast$output)]]$bp.Vt$breakpoints[i]]
      
    }
  }



  trend<-data.frame(as.character(as.Date(unlist(trend))))
  names(trend)<- c("Trend")
  season <-list(list())
  if(is.na(bfast$output[[length(bfast$output)]]$bp.Wt$breakpoints)[1]==TRUE){
    season <-NA
  }  else {
    for(i in 1:length(bfast$output[[length(bfast$output)]]$bp.Wt$breakpoints)){
      season[[i]] <- date[bfast$output[[length(bfast$output)]]$bp.Wt$breakpoints[i]]
      
    }
  }


  season<-data.frame(as.character(as.Date(unlist(season))))
  names(season)<- c("Season")
  breaks <- list(trend, season)
  names(breaks) <- c("Trend", "Season")
  return(breaks)
}


basn.breaks <- bfast.breaks(basn.fit, ndvi.date)
coll.breaks <- bfast.breaks(coll.fit, ndvi.date)
east.breaks <- bfast.breaks(east.fit, ndvi.date)
grav.breaks <- bfast.breaks(grav.fit, ndvi.date)
ibpe.breaks <- bfast.breaks(ibpe.fit, ndvi.date)
nort.breaks <- bfast.breaks(nort.fit, ndvi.date)
rabb.breaks <- bfast.breaks(rabb.fit, ndvi.date)
sand.breaks <- bfast.breaks(sand.fit, ndvi.date)
smal.breaks <- bfast.breaks(smal.fit, ndvi.date)
summ.breaks <- bfast.breaks(summ.fit, ndvi.date)
tayl.breaks <- bfast.breaks(tayl.fit, ndvi.date)
tobo.breaks <- bfast.breaks(tobo.fit, ndvi.date)
well.breaks <- bfast.breaks(well.fit, ndvi.date)
west.breaks <- bfast.breaks(west.fit, ndvi.date)

basn.breaks 
coll.breaks 
east.breaks 
grav.breaks 
ibpe.breaks 
nort.breaks
rabb.breaks
sand.breaks
smal.breaks
summ.breaks
tayl.breaks
tobo.breaks
well.breaks
west.breaks


###########################################################################
#function to extract CI interval dates
CI.interval.extract <- function(bfast,date, site){
      CI.trend <- data.frame(date[bfast$output[[length(bfast$output)]]$ci.Vt$confint[,1]],date[bfast$output[[length(bfast$output)]]$ci.Vt$confint[,2]],date[bfast$output[[length(bfast$output)]]$ci.Vt$confint[,3]])
      colnames(CI.trend) <- c("CI.start","Break", "CI.end")
      rownames(CI.trend) <- c(paste0(site,rep("T",nrow(CI.trend)),seq(1,nrow(CI.trend),1)))
      
      if(is.na(bfast$output[[length(bfast$output)]]$ci.Wt)[1]==TRUE){
      CI.final <- CI.trend
      } else {
      CI.season <- data.frame(date[bfast$output[[length(bfast$output)]]$ci.Wt$confint[,1]],date[bfast$output[[length(bfast$output)]]$ci.Wt$confint[,2]],date[bfast$output[[length(bfast$output)]]$ci.Wt$confint[,3]])
      colnames(CI.season) <- c("CI.start","Break", "CI.end")
      rownames(CI.season) <- c(paste0(site,rep("S",nrow(CI.season)),seq(1,nrow(CI.season),1)))
      CI.final <- rbind(CI.trend, CI.season)
      }
      return(CI.final)
  }


basn.CI.bounds <- CI.interval.extract( basn.fit,ndvi.date, "basn") 
cali.CI.bounds <- CI.interval.extract( cali.fit,ndvi.date, "cali")
coll.CI.bounds <- CI.interval.extract( coll.fit,ndvi.date, "coll") 
east.CI.bounds <- CI.interval.extract( east.fit,ndvi.date, "east") 
grav.CI.bounds <- CI.interval.extract( grav.fit,ndvi.date, "grav") 
ibpe.CI.bounds <- CI.interval.extract( ibpe.fit,ndvi.date, "ibpe") 
nort.CI.bounds <- CI.interval.extract( nort.fit,ndvi.date, "nort")
rabb.CI.bounds <- CI.interval.extract( rabb.fit,ndvi.date, "rabb")
sand.CI.bounds <- CI.interval.extract( sand.fit,ndvi.date, "sand")
smal.CI.bounds <- CI.interval.extract( smal.fit,ndvi.date, "smal")
summ.CI.bounds <- CI.interval.extract( summ.fit,ndvi.date, "summ")
tayl.CI.bounds <- CI.interval.extract( tayl.fit,ndvi.date, "tayl")
tobo.CI.bounds <- CI.interval.extract( tobo.fit,ndvi.date, "tobo")
well.CI.bounds <- CI.interval.extract( well.fit,ndvi.date, "well")
west.CI.bounds <- CI.interval.extract( west.fit,ndvi.date, "west")

CI.bounds <- rbind(basn.CI.bounds,  
                    coll.CI.bounds,  
                    cali.CI.bounds, 
                    east.CI.bounds,  
                    grav.CI.bounds,  
                    ibpe.CI.bounds,  
                    nort.CI.bounds, 
                    rabb.CI.bounds, 
                    sand.CI.bounds, 
                    smal.CI.bounds, 
                    summ.CI.bounds, 
                    tayl.CI.bounds, 
                    tobo.CI.bounds, 
                    well.CI.bounds, 
                    west.CI.bounds)

write.csv(CI.bounds, "JRN_CI_bounds.csv")
