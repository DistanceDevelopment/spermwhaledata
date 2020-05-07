# import sperm whale data into R

# Data from NOAA North East and South East Fisheries Science Centers

# re-write with sf/stars May 2020

library(sf)
library(stars)
library(lubridate)
options(stringsAsFactors = FALSE)

# path to geodatabase
gdb_path <- "rawdata/Analysis.gdb"
# export path
export_path <- "distance_import/"

## line transects/effort table

# effort data as csv
segs <- read_sf(gdb_path, "Segment_Centroids")
segs <- st_drop_geometry(segs)
segs$x <- segs$POINT_X
segs$y <- segs$POINT_Y
segs$Effort <- segs$Length
segs$Sample.Label <- segs$SegmentID
segs$CenterTime <- ymd_hms(segs$CenterTime)

# get rid of nuisance columns
segs$Length <- segs$coords.x1 <- segs$coords.x2 <-
  segs$POINT_X <- segs$POINT_Y <- NULL


# also save segments as shapefile
segs_shp <- read_sf(gdb_path, "Segments")
# this should be all we need...
write_sf(segs_shp, dsn=paste0(export_path, "segments.shp"),
         driver="ESRI Shapefile" )

# add the survey ID from the shapefile to segs
dd <- unique(st_drop_geometry(segs_shp)[,c("SegmentID", "FIRST_Survey")])
names(dd) <- c("Sample.Label", "Survey")
segs <- merge(segs, dd, by="Sample.Label")




## add in the Beaufort to the segments
endat <- as.data.frame(read_sf(gdb_path, "EN_Trackline2"))
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

# now do the same for the GU segments
gudat <- as.data.frame(read_sf(gdb_path, "GU_Trackline2"))
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
transects <- read_sf(gdb_path, "Tracklines")
write_sf(transects, dsn=paste0(export_path, "transects.shp"),
         driver="ESRI Shapefile" )


## observation and distance data

# obs table
obs <- read_sf(gdb_path, "Sightings")
obs <- st_drop_geometry(obs)
obs$distance <- obs$Distance
obs$object <- obs$SightingID
obs$Sample.Label <- obs$SegmentID
obs$size <- obs$GroupSize
obs$observer <- obs$detected <- 1
obs$Beaufort <- obs$SeaState

# get rid of nuisance columns
obs$Distance <- obs$SightingID <- obs$SegmentID <- obs$GroupSize <- obs$SeaState <- NULL
write.csv(obs, file=paste0(export_path, "obs.csv"), row.names=FALSE)

# distance data
write.csv(obs, file=paste0(export_path, "dist.csv"), row.names=FALSE)

## study area
study_area <- read_sf(gdb_path, "Study_Area")
write_sf(study_area, dsn=paste0(export_path, "studyarea.shp"),
         driver="ESRI Shapefile")


## prediction grid

# list of files to load
raster_files <- paste0("rawdata/",
                       c("Covariates_for_Study_Area/Depth.img",
                         "Covariates_for_Study_Area/GLOB/CMC/CMC0.2deg/analysed_sst/2004/20040624-CMC-L4SSTfnd-GLOB-v02-fv02.0-CMC0.2deg-analysed_sst.img",
                         "Covariates_for_Study_Area/VGPM/Rasters/vgpm.2004153.hdf.gz.img",
                         "Covariates_for_Study_Area/DistToCanyonsAndSeamounts.img",
                         "Covariates_for_Study_Area/Global/DT\ all\ sat/MSLA_ke/2004/MSLA_ke_2004176.img"
                       ))

# read all into one object
predictorStack <- read_stars(raster_files)

# rename the layers in our stack to match those in the model
names(predictorStack) <- c("Depth","SST","NPP", "DistToCAS", "EKE")

# coerce to data.frame
predgrid <- as.data.frame(predictorStack, xy=TRUE)
# size of grid cells
predgrid$off.set <- (10*1000)^2



# now create a shapefile for distance
# from the Distance manual this requires a set of points not polygons
pred_shp <- st_cast(st_as_sf(predictorStack), "POINT")
write_sf(pred_shp, dsn=paste0(export_path, "predgrid.shp"),
         driver="ESRI Shapefile" )

# as csv too
predgrid_csv <- predgrid
predgrid_csv$LinkID <- seq(from=1, to=nrow(predgrid_csv))
write.csv(predgrid, file=paste0(export_path, "predgrid.csv"),
          row.names=FALSE, quote = FALSE)


# raster is rectangular, prediction grid is not
# so remove the NA values were we won't predict
predgrid <- predgrid[!is.na(predgrid$Depth), ]


# now save some things for R
dist <- obs
save(segs, obs, dist, study_area, predgrid, file="R_import/spermwhale.RData")

