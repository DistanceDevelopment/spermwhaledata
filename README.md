# spermwhaledata

Sperm whale line transect survey data in R and Distance for Windows format.

Data provided by Debi Palka (NOAA North East Fisheries Science Center) and Lance Garrison (NOAA South East Fisheries Science Center). Initial data processing by Jason Roberts (Marine Geospatial Ecology Lab, Duke University).

More information available in the `background/` directory.

# Data notes

- Data provided as an example dataset for teaching and should not be used for abundance estimation.
- Beaufort is interpolated on segments with no detections using the `approx()` function in R.
