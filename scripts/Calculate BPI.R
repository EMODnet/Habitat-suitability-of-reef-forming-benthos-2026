setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(terra)
library(MultiscaleDTM)

### ------------------------------------------------------------
### Load and reproject bathymetry
### ------------------------------------------------------------

bathy <- rast("./data/raw_data/bathymetry_emodnet_80m.tif")

# Reproject to ETRS89 / LAEA Europe (units = metres)
bathy_m <- project(
  bathy,
  "EPSG:3035",
  method = "bilinear"
)

### ------------------------------------------------------------
### Helper functions
### ------------------------------------------------------------

# km → pixels
km_to_pixels <- function(r, km) {
  round((km * 1000) / mean(res(r)))
}

# Harmonise DEM to remove hard coastal / island boundaries
# (nearest-neighbour propagation, no smoothing)
harmonise_dem_for_bpi <- function(dem) {
  
  dem_work <- dem
  
  # Iteratively fill NA using focal mean (built-in, terra-safe)
  while (any(is.na(values(dem_work)))) {
    
    dem_work <- cover(
      dem_work,
      focal(
        dem_work,
        w = 3,
        fun = mean,
        na.rm = TRUE,
        expand = TRUE
      )
    )
  }
  
  return(dem_work)
}

### ============================================================
### SCALE 1: 5–10 km BPI  (~650 m resolution)
### ============================================================

# Aggregate bathymetry
bathy_650 <- aggregate(bathy_m, fact = 8, fun = mean)

# Preserve original marine mask
marine_650 <- bathy_650 < 0

# Harmonise DEM for BPI neighbourhoods
bathy_650_h <- harmonise_dem_for_bpi(bathy_650)

# Radii in pixels
rad_5km  <- km_to_pixels(bathy_650_h, 5)
rad_10km <- km_to_pixels(bathy_650_h, 10)

# --- BPI 5 km ---
bpi5 <- BPI(bathy_650_h, w = c(1, rad_5km))
bpi5 <- mask(bpi5, marine_650)
writeRaster(bpi5, "./data/derived_data/BPI5.tif", overwrite = TRUE)

# --- BPI 10 km ---
bpi10 <- BPI(bathy_650_h, w = c(1, rad_10km))
bpi10 <- mask(bpi10, marine_650)
writeRaster(bpi10, "./data/derived_data/BPI10.tif", overwrite = TRUE)

rm(bathy_650, bathy_650_h, marine_650, bpi5, bpi10); gc()

### ============================================================
### SCALE 2: 30–50 km BPI  (~1.3 km resolution)
### ============================================================

bathy_1300 <- aggregate(bathy_m, fact = 16, fun = mean)
marine_1300 <- bathy_1300 < 0
bathy_1300_h <- harmonise_dem_for_bpi(bathy_1300)

rad_30km <- km_to_pixels(bathy_1300_h, 30)
rad_50km <- km_to_pixels(bathy_1300_h, 50)

# --- BPI 30 km ---
bpi30 <- BPI(bathy_1300_h, w = c(1, rad_30km))
bpi30 <- mask(bpi30, marine_1300)
writeRaster(bpi30, "./data/derived_data/BPI30.tif", overwrite = TRUE)

# --- BPI 50 km ---
bpi50 <- BPI(bathy_1300_h, w = c(1, rad_50km))
bpi50 <- mask(bpi50, marine_1300)
writeRaster(bpi50, "./data/derived_data/BPI50.tif", overwrite = TRUE)

rm(bathy_1300, bathy_1300_h, marine_1300, bpi30, bpi50); gc()

### ============================================================
### SCALE 3: 75–100 km BPI  (~2.6 km resolution)
### ============================================================

bathy_2600 <- aggregate(bathy_m, fact = 32, fun = mean)
marine_2600 <- bathy_2600 < 0
bathy_2600_h <- harmonise_dem_for_bpi(bathy_2600)

rad_75km  <- km_to_pixels(bathy_2600_h, 75)
rad_100km <- km_to_pixels(bathy_2600_h, 100)

# --- BPI 75 km ---
bpi75 <- BPI(bathy_2600_h, w = c(1, rad_75km))
bpi75 <- mask(bpi75, marine_2600)
writeRaster(bpi75, "./data/derived_data/BPI75.tif", overwrite = TRUE)

# --- BPI 100 km ---
bpi100 <- BPI(bathy_2600_h, w = c(1, rad_100km))
bpi100 <- mask(bpi100, marine_2600)
writeRaster(bpi100, "./data/derived_data/BPI100.tif", overwrite = TRUE)

rm(bathy_2600, bathy_2600_h, marine_2600, bpi75, bpi100); gc()
