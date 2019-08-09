library(MODIS)
library(rgdal)
library(doMPI)

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
getHdf(product="MOD13Q1",begin="2002017", end="2003001", tileH=9,tileV=5, collection='006')


## Read in the Jornada boundary shapefile (specify source directory)
JRN <- readOGR(dsn="Infrastructure2.gdb",layer="JER_Bdry_25Aug15")
JRN.point <- readOGR(dsn="Infrastructure2.gdb",layer="NPP_Sites_study011_Centroid")



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
getSds(HdfName="MOD13Q1.A2000065.h09v05.006.2015136022550.hdf")

# $SDSnames
# [1] "250m 16 days NDVI"                       [2] "250m 16 days EVI"                      
# [3] "250m 16 days VI Quality"                 [4] "250m 16 days red reflectance"          
# [5] "250m 16 days NIR reflectance"            [6] "250m 16 days blue reflectance"         
# [7] "250m 16 days MIR reflectance"            [8] "250m 16 days view zenith angle"        
# [9] "250m 16 days sun zenith angle"           [10] "250m 16 days relative azimuth angle"   
# [11] "250m 16 days composite day of the year" [12]"250m 16 days pixel reliability" 

#run MODIS function runGdal to extract MODIS Data
#MOD13Q1
runGdal(job="JRN_Modis_VI", product="MOD13Q1", collection='006', extent=JRN_raster, begin="2000065", SDSstring="111111")

#set the path to get the correct data on the source directory: 
#path to 250M VI's
path <- paste(options("MODIS_outDirPath"),"/JRN_Modis_VI",sep="") 
ndvi <- preStack(path=path, pattern="*_NDVI.tif$")
JRN_NDVI_brick<- brick(stack(ndvi),  filename='JRN_NDVI_brick', overwrite=TRUE)


#Raster layers need to be in a vector list to parallelize accross a cluster.
stack2list <- function(stack){
  list<-vector("list", nlayers(stack))
  for(i in 1:nlayers(stack)) {
    list[[i]] <- stack[[i]]
  }
  return(list)
}

JRN_NDVI_list<- stack2list(stack(ndvi))


##########################################################################################################
# I downloaded all the tifs from appEEARS, start here
path = 'C:/Users/echriste/Documents/Projects/Active/bfast/MODIS_appeears/NDVI'
ndvi <- preStack(path=path, pattern='*_NDVI_*')
JRN_NDVI_brick <- brick(stack(ndvi), filename='JRN_NDVI_brick', overwrite=T)

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
JRN_QA_list <- JRN_NDVI_list
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
  print(i)
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
JRN_NDVI_mask <- brick(stack(JRN_NDVI_mask),  filename='JRN_NDVI_mask', overwrite=TRUE)

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


# Create one big data frame
JRN_NDVI <- cbind(BASN.ndvi.m,
                  CALI.ndvi.m,
                  COLL.ndvi.m,
                  EAST.ndvi.m,
                  GRAV.ndvi.m,
                  IBPE.ndvi.m,
                  NORT.ndvi.m,
                  RABB.ndvi.m,
                  SAND.ndvi.m,
                  SMAL.ndvi.m,
                  SUMM.ndvi.m,
                  TAYL.ndvi.m,
                  TOBO.ndvi.m,
                  WELL.ndvi.m,
                  WEST.ndvi.m)
JRN_NDVI <- data.frame(JRN_NDVI)

dates = str_extract(rownames(JRN_NDVI),'[0-9]{7}')
JRN_NDVI$date = dates


write.csv(JRN_NDVI,'C:/Users/echriste/Documents/Projects/Active/bfast/JRN_NDVI_collection6.csv')
