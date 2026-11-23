setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(copernicusR)
library(terra)

setup_copernicus(
  username = Sys.getenv("COPERNICUS_USERNAME"), 
  password = Sys.getenv("COPERNICUS_PASSWORD")
)

# --- Mean salinity ---

# Download
copernicus_download(
  dataset_id = "cmems_mod_nws_phy-sal_anfc_1.5km-3D_P1D-m",
  dataset_version = "202511",
  variables = "so",
  start_date = "2025-01-01",
  end_date = "2026-01-01",
  bbox = c(-12, 12, 48, 62),
  output_file = "salinity_3d.nc"
)

# Move to my own folder
file.rename(
  from = "salinity_3d.nc",
  to = file.path(getwd(), "Environment", "salinity_3d.nc")
)

# Calculate mean
nc <- rast("./Environment/salinity_3d.nc")
mean_salinity <- mean(nc, na.rm = TRUE)

# Save as GeoTIFF
writeRaster(
  mean_salinity,
  filename = "./Environment/mean_salinity.tif",
  overwrite = TRUE
)

# --- Mean bottom temperature ---

# Download
copernicus_download(
  dataset_id = "cmems_mod_nws_phy-bottomt_anfc_1.5km-2D_P1D-m",
  dataset_version = "202511",
  variables = "bottomT",
  start_date = "2025-01-01",
  end_date = "2026-01-01",
  bbox = c(-12, 12, 48, 62),
  output_file = "bottom_temperature.nc"
)

# Move to my own folder
file.rename(
  from = "bottom_temperature.nc",
  to = file.path(getwd(), "Environment", "bottom_temperature.nc")
)

# Calculate mean
nc <- rast("./Environment/bottom_temperature.nc")
mean_temp <- mean(nc, na.rm = TRUE)

# Save as GeoTIFF
writeRaster(
  mean_temp,
  filename = "./Environment/mean_bottom_temperature.tif",
  overwrite = TRUE
)

# --- Stratification ---

# Download
copernicus_download(
  dataset_id = "cmems_mod_nws_phy-mld_anfc_1.5km-2D_P1D-m",
  dataset_version = "202511",
  variables = "mlotst",
  start_date = "2025-01-01",
  end_date = "2026-01-01",
  bbox = c(-12, 12, 48, 62),
  output_file = "mld.nc"
)

# Move to my own folder
file.rename(
  from = "mld.nc",
  to = file.path(getwd(), "Environment", "mld.nc")
)

# Calculate mean
nc <- rast("./Environment/mld.nc")
mixed_layer_depth <- mean(nc, na.rm = TRUE)

# Save as GeoTIFF
writeRaster(
  mixed_layer_depth,
  filename = "./Environment/mixed_layer_depth.tif",
  overwrite = TRUE
)

# --- Bathymetry (from EMODnet) ---

files <- c(
  "./Environment/Bathymetry/D3_2024.tif",
  "./Environment/Bathymetry/D4_2024.tif",
  "./Environment/Bathymetry/D5_2024.tif",
  "./Environment/Bathymetry/E3_2024.tif",
  "./Environment/Bathymetry/E4_2024.tif",
  "./Environment/Bathymetry/E5_2024.tif"
)

# Create a raster collection
tiles <- sprc(files)

# Mosaic all tiles together
merged <- mosaic(tiles)

# Crop to your area
cropped <- crop(merged, ext(-12, 12, 48, 62))

writeRaster(cropped,
            "./Environment/bathymetry_emodnet_115m.tif",
            overwrite = TRUE)
