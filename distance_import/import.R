# importer for the GIS files from Jason Roberts into Distance
# compatible format

# David L Miller 2017

library(rgdal)

# base path
base_path <- "../rawdata/"

## line transects/effort table

# effort data as csv
segs <- readOGR(paste0(base_path, "Analysis.gdb"),"Segment_Centroids")
segs <- as.data.frame(segs)
segs$x <- segs$POINT_X
segs$y <- segs$POINT_Y
segs$Effort <- segs$Length
segs$Sample.Label <- segs$SegmentID

# get rid of nuisance columns
segs$Length <- segs$coords.x1 <- segs$coords.x2 <-
  segs$POINT_X <- segs$POINT_Y <- NULL
write.csv(segs, file="segments.csv")

# segments as shapefile
segs_shp <- readOGR(paste0(base_path, "Analysis.gdb"),"Segments")
# this should be all we need...
writeOGR(segs_shp, "segments.shp", "data", "ESRI Shapefile" )


## observation and distance data

# obs table
obs <- readOGR(paste0(base_path, "Analysis.gdb"), "Sightings")
obs <- as.data.frame(obs)
obs$distance <- obs$Distance
obs$object <- obs$SightingID
obs$Sample.Label <- obs$SegmentID
obs$size <- obs$GroupSize

# get rid of nuisance columns
obs$Distance <- obs$SightingID <- obs$SegmentID <- obs$GroupSize  <- segs$coords.x1 <- segs$coords.x2 <- NULL
write.csv(obs, file="obs.csv")

# distance data
write.csv(obs, file="dist.csv")

## study area
study_area <- readOGR(paste0(base_path, "Analysis.gdb"), "Study_Area")
writeOGR(study_area, "studyarea.shp", "data", "ESRI Shapefile" )


## prediction grid
# from the Distance manual this requires a **point** not polygon set
library(raster)

# build a predictor stack
predictorStack <- stack(paste0(base_path,
                        c("Covariates_for_Study_Area/Depth.img",
                          "Covariates_for_Study_Area/GLOB/CMC/CMC0.2deg/analysed_sst/2004/20040602-CMC-L4SSTfnd-GLOB-v02-fv02.0-CMC0.2deg-analysed_sst.img",
                          "Covariates_for_Study_Area/VGPM/Rasters/vgpm.2004153.hdf.gz.img",
                          "Covariates_for_Study_Area/DistToCanyonsAndSeamounts.img",
                          "Covariates_for_Study_Area/Global/DT\ all\ sat/MSLA_ke/2004/MSLA_ke_2004153.img"
                          )))
# rename the layers in our stack to match those in the model 
names(predictorStack) <- c("Depth","SST","NPP", "DistToCAS", "EKE")

predgrid <- as.data.frame(predictorStack, xy=TRUE)
predgrid$off.set <- (10*1000)^2

predgrid <- predgrid[!is.na(predgrid$Depth),]

coords <- cbind(predgrid$x, predgrid$y)
sp <- SpatialPoints(coords)
# make spatial data frame
spdf <- SpatialPointsDataFrame(coords, predgrid)

writeOGR(spdf, "predgrid.shp", "data", "ESRI Shapefile" )



