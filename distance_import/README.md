Data to import into Distance
============================

Data in this folder can be imported straight into Distance.

`import.R` does the work and looks for files in `../rawdata` to convert. It will produce the following:

- `dist.csv` distances to be used for detection function analysis
- `obs.csv` observation table link file that links segments to detected animals
- `segments.csv` effort data on the segments (locations are segment centroids)
- `segments.shp` shapefile corresponding to the exact locations of the segment lines
- `predgrid.shp` prediction grid shape file (points) with data giving the covariate values (dynamic covariates fix at values from ~2004-06-02).
- `studyarea.shp` polygon shapefile with the study area (note: not the US EEZ)

