# importer for the GIS files from Jason Roberts into Distance and R
# compatible formats

# David L Miller 2017
#   bugs added by Rexstad

library(rgdal)
library(lubridate)
options(stringsAsFactors = FALSE)

# base path
base_path <- "rawdata/"
# export path
export_path <- "distance_import/"

## line transects/effort table

# effort data as csv
segs <- readOGR(paste0(base_path, "Analysis.gdb"),"Segment_Centroids")
segs <- as.data.frame(segs)
segs$x <- segs$POINT_X
segs$y <- segs$POINT_Y
segs$Effort <- segs$Length
segs$Sample.Label <- segs$SegmentID
segs$CenterTime <- ymd_hms(segs$CenterTime)

# get rid of nuisance columns
segs$Length <- segs$coords.x1 <- segs$coords.x2 <-
  segs$POINT_X <- segs$POINT_Y <- NULL

# segments as shapefile
segs_shp <- readOGR(paste0(base_path, "Analysis.gdb"),"Segments")
# this should be all we need...
writeOGR(segs_shp, paste0(export_path, "segments.shp"), "data",
         "ESRI Shapefile" )

# add in the survey ID from the shapefile
dd <- unique(segs_shp@data[,c("SegmentID", "FIRST_Survey")])
names(dd) <- c("Sample.Label", "Survey")
segs <- merge(segs, dd, by="Sample.Label")

# add in the Beaufort

endat <- as.data.frame(readOGR(paste0(base_path, "Analysis.gdb"),"EN_Trackline2"))
# make a time window
endat$en04_effort_csv_begdatetime <- ymd_hms(endat$en04_effort_csv_begdatetime)
endat$en04_effort_csv_enddatetime <- ymd_hms(endat$en04_effort_csv_enddatetime)
endat$window <- endat$en04_effort_csv_begdatetime %--% endat$en04_effort_csv_enddatetime

# get the EN survey
segsEN <- subset(segs, Survey=="en04395")
# make a beaufort column
segsEN$Beaufort <- NA

# loop over the segments
for(j in 1:nrow(segsEN)){
  # which interval does this segment lie in?
  ind <- which(segsEN$CenterTime[j] %within% endat$window)
  if(length(ind)==1){
    segsEN$Beaufort[j] <- endat$en04_effort_csv_beaufort[ind]
  }
}

# some were missing, interpolate
naind <- is.na(segsEN$Beaufort)
napred <- approx(as.numeric(segsEN$CenterTime)[!naind], segsEN$Beaufort[!naind],
                 as.numeric(segsEN$CenterTime)[naind])

# check plot
#plot(segsEN$CenterTime[!naind], segsEN$Beaufort[!naind], pch=19)
#points(as.numeric(segsEN$CenterTime)[naind], napred$y, pch=21)

#now for GU
gudat <- as.data.frame(readOGR(paste0(base_path, "Analysis.gdb"),"GU_Trackline2"))
# make a time window
gudat$DateTime1 <- ymd_hms(gudat$DateTime1)
gudat$DateTime2 <- ymd_hms(gudat$DateTime2)
gudat$window <- gudat$DateTime1 %--% gudat$DateTime2
# get the EN survey
segsGU <- subset(segs, Survey=="GU0403")
# make a beaufort column
segsGU$Beaufort <- NA

# loop over the segments
for(j in 1:nrow(segsGU)){
  # which interval does this segment lie in?
  ind <- which(segsGU$CenterTime[j] %within% gudat$window)
  if(length(ind)==1){
    segsGU$Beaufort[j] <- gudat$SeaState[ind]
  }
}

# some were missing, interpolate
naind <- is.na(segsGU$Beaufort)
napred <- approx(as.numeric(segsGU$CenterTime)[!naind], segsGU$Beaufort[!naind],
                 as.numeric(segsGU$CenterTime)[naind], method="constant")

# check plot
#plot(segsGU$CenterTime[!naind], segsGU$Beaufort[!naind], pch=19)
#points(as.numeric(segsGU$CenterTime)[naind], napred$y, pch=21)


# smoosh them back together again
segs <- rbind(segsEN, segsGU)

# write the csv
write.csv(segs, file=paste0(export_path, "segments.csv"), row.names=FALSE)

## transects
transects <- readOGR(paste0(base_path, "Analysis.gdb"), "Tracklines")
writeOGR(transects, paste0(export_path, "transects.shp"),
         "data", "ESRI Shapefile" )


## observation and distance data

# obs table
obs <- readOGR(paste0(base_path, "Analysis.gdb"), "Sightings")
obs <- as.data.frame(obs)
obs$distance <- obs$Distance
obs$object <- obs$SightingID
obs$Sample.Label <- obs$SegmentID
obs$size <- obs$GroupSize
obs$observer <- obs$detected <- 1
obs$Beaufort <- obs$SeaState

# get rid of nuisance columns
obs$Distance <- obs$SightingID <- obs$SegmentID <- obs$GroupSize <- obs$SeaState<- NULL
write.csv(obs, file=paste0(export_path, "obs.csv"), row.names=FALSE)

# distance data
write.csv(obs, file=paste0(export_path, "dist.csv"), row.names=FALSE)

## study area
study_area <- readOGR(paste0(base_path, "Analysis.gdb"), "Study_Area")
writeOGR(study_area, paste0(export_path, "studyarea.shp"),
         "data", "ESRI Shapefile")


## prediction grid
# from the Distance manual this requires a **point** not polygon set
library(raster)

# build a predictor stack
predictorStack <- stack(paste0(base_path,
                        c("Covariates_for_Study_Area/Depth.img",
                          "Covariates_for_Study_Area/GLOB/CMC/CMC0.2deg/analysed_sst/2004/20040624-CMC-L4SSTfnd-GLOB-v02-fv02.0-CMC0.2deg-analysed_sst.img",
                          "Covariates_for_Study_Area/VGPM/Rasters/vgpm.2004153.hdf.gz.img",
                          "Covariates_for_Study_Area/DistToCanyonsAndSeamounts.img",
                          "Covariates_for_Study_Area/Global/DT\ all\ sat/MSLA_ke/2004/MSLA_ke_2004176.img"
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

writeOGR(spdf, paste0(export_path, "predgrid.shp"), "data", "ESRI Shapefile" )

# as csv too
predgrid <- round(as.data.frame(spdf@data),4)
predgrid$LinkID <- seq(from=1, to=dim(spdf@data)[1])
write.csv(predgrid, file=paste0(export_path, "predgrid.csv"),
          row.names=FALSE, quote = FALSE)

# now save some things for R
dist <- obs
save(segs, obs, dist, study_area, predgrid, file="R_import/spermwhale.RData")

