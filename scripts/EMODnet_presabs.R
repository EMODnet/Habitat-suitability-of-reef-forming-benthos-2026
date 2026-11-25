# --- EMODnet macrobenthos presence/absence ---

# Dependencies: tidyverse, sf, readr
required <- c("tidyverse", "sf", "readr")
new_pkgs <- required[!required %in% installed.packages()[,"Package"]]
if(length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)
invisible(lapply(required, library, character.only = TRUE))

# --- 0) Species of interest (WoRMS AphiaIDs) ---
sp_list <- tibble::tribble(
  ~aphia_id, ~species_name,
  131495,    "Lanice conchilega"    # WoRMS; for Sabellaria spinulosa: 130867
)

# --- 1) Read EMODnet derived presence/absence tables (lightweight)
# These tables are published by EMODnet's benthos-trends product and derived from the official PA NetCDF.
events_url <- "https://raw.githubusercontent.com/EMODnet/EMODnet-Biology-benthos-trends/master/data/derived_data/sample_events.csv"
pres_url   <- "https://raw.githubusercontent.com/EMODnet/EMODnet-Biology-benthos-trends/master/data/derived_data/pres_df.csv"

message("Downloading derived event table...")
sample_events <- readr::read_csv(events_url, show_col_types = FALSE)
message("Downloading derived presence table...")
pres_df       <- readr::read_csv(pres_url,   show_col_types = FALSE)

# Expect columns in sample_events: sample_id, lat, lon, date, year, includes_na
# Expect columns in pres_df: sample_id, aphia_id
stopifnot(all(c("sample_id","lat","lon","date","year","includes_na") %in% names(sample_events)))
stopifnot(all(c("sample_id","aphia_id") %in% names(pres_df)))

# --- 2) Filter events to complete macrobenthos inventories (valid PA inference)
events_pa <- sample_events %>%
  dplyr::filter(includes_na == FALSE) %>%  # critical for valid absences
  dplyr::mutate(event_id = as.character(sample_id)) %>%
  dplyr::filter(is.finite(lat), is.finite(lon))

pres_pa <- pres_df %>%
  dplyr::mutate(event_id = as.character(sample_id))

# --- 3) Sanity check: global presence counts (no spatial clip yet)
pres_counts_global <- pres_pa %>%
  dplyr::filter(aphia_id %in% sp_list$aphia_id) %>%
  dplyr::count(aphia_id, name = "n_present_global")
print(pres_counts_global)

if(any(pres_counts_global$n_present_global == 0)) {
  warning("One or more target AphiaIDs have zero presences globally in pres_df. Check IDs and columns.")
}

# --- 4) Build global presence/absence (denominator = all PA-suitable events)
make_pa_global <- function(one_aphia, one_name) {
  pres_sp <- pres_pa %>% dplyr::filter(aphia_id == one_aphia) %>% dplyr::distinct(event_id)
  events_pa %>%
    dplyr::mutate(
      presence_absence = as.integer(event_id %in% pres_sp$event_id),
      species_aphia = one_aphia,
      species_name  = one_name
    )
}
pa_global <- purrr::pmap_df(sp_list, ~ make_pa_global(..1, ..2))
message("\nGlobal PA counts (no spatial clip):")
print(pa_global %>% dplyr::count(species_name, presence_absence))

# --- 5) Spatial clip: Greater North Sea bounding box
# Official ICES bounds: Min Lat 48, Max Lat 62; Min Lon -5, Max Lon 13.0629 (Marine Regions).
# We expand slightly to avoid trimming coastal/western events (e.g., -6.x near UK west coast).
north_sea_wkt <- "POLYGON((-11 48, 13 48, 13 62, -11 62, -11 48))"
north_sea_poly <- sf::st_sfc(sf::st_as_sfc(north_sea_wkt), crs = 4326)

events_sf <- sf::st_as_sf(events_pa, coords = c("lon","lat"), crs = 4326, remove = FALSE)
events_ns <- events_sf[sf::st_within(events_sf, north_sea_poly, sparse = FALSE), ] %>%
  sf::st_drop_geometry()

# --- 6) Recompute PA for North Sea only
make_pa_ns <- function(one_aphia, one_name) {
  pres_sp <- pres_pa %>% dplyr::filter(aphia_id == one_aphia) %>% dplyr::distinct(event_id)
  events_ns %>%
    dplyr::mutate(
      presence_absence = as.integer(event_id %in% pres_sp$event_id),
      species_aphia = one_aphia,
      species_name  = one_name
    )
}
pa_ns <- purrr::pmap_df(sp_list, ~ make_pa_ns(..1, ..2))

message("\nNorth Sea PA counts (widened box):")
print(pa_ns %>% dplyr::count(species_name, presence_absence))

# --- 7) Export per-species CSVs
out_dir <- file.path(getwd(), "emodnet_pa_out")
if(!dir.exists(out_dir)) dir.create(out_dir)

write_species_csv <- function(df, sp_name) {
  fn <- file.path(out_dir, paste0(gsub(" ", "_", tolower(sp_name)), "_greater_north_sea_presence_absence.csv"))
  df %>%
    dplyr::filter(species_name == sp_name) %>%
    dplyr::select(event_id, date, lon, lat, species_aphia, species_name, presence_absence) %>%
    readr::write_csv(fn)
  message("Wrote: ", fn)
}

unique(pa_ns$species_name) %>% purrr::walk(~ write_species_csv(pa_ns, .x))

# --- 8) Summary table
summary_tbl <- pa_ns %>%
  dplyr::group_by(species_name) %>%
  dplyr::summarise(
    n_events_ns = dplyr::n(),
    n_present_ns = sum(presence_absence),
    n_absent_ns = n_events_ns - n_present_ns,
    .groups = "drop"
  )
message("\nGreater North Sea prevalence (widened box):")
print(summary_tbl)

# Visualise

# Packages
library(leaflet)
library(leaflet.extras)
library(sf)
library(dplyr)

leaflet(subset(pa_ns, species_name == "Lanice conchilega")) %>%
  addTiles() %>%
  addCircleMarkers(
    color = ~ifelse(presence_absence == 1, "green", "red"),
    radius = 2, stroke = FALSE, fillOpacity = 0.8
  ) %>%
  addMeasure(
    position = "topleft",
    primaryLengthUnit = "meters",
    secondaryLengthUnit = "kilometers",
    primaryAreaUnit   = "sqmeters",
    secondaryAreaUnit = "hectares"
  )

# Add additional OBIS points

# --- OBIS download ---
library(robis)
library(sf)
library(dplyr)

aphia <- 131495  # Lanice conchilega

obis_raw <- occurrence(taxonid = aphia)

# Keep only valid coordinates
obis_sf <- obis_raw %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# --- Limit OBIS to the EMODnet bounding box ---
obis_sf <- obis_sf[st_within(obis_sf, north_sea_poly, sparse = FALSE), ]

# --- Convert EMODnet PA events to sf ---
emod_sf <- pa_ns %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# --- Identify OBIS points overlapping EMODnet events (1 km buffer) ---

library(nngeo)

buffer_m <- 1000  # 1 km

# Nearest EMODnet event for each OBIS point
obis_sf_m  <- st_transform(obis_sf,  3857)
emod_sf_m  <- st_transform(emod_sf,  3857)

nn <- st_nn(obis_sf_m, emod_sf_m, k = 1, returnDist = TRUE)

nearest_event <- sapply(nn$nn, `[`, 1)
nearest_dist  <- sapply(nn$dist, `[`, 1)

# OBIS points within 1 km → flip absences
obis_overlapping <- obis_sf[nearest_dist <= buffer_m, ]
events_to_flip <- emod_sf$event_id[nearest_event[nearest_dist <= buffer_m]]

pa_ns_updated <- pa_ns %>%
  mutate(presence_absence = ifelse(event_id %in% events_to_flip, 1, presence_absence))

# OBIS points outside 1 km → add as new presence-only events
obis_new <- obis_sf[nearest_dist > buffer_m, ]

obis_new_df <- obis_new %>%
  st_drop_geometry() %>%
  transmute(
    event_id = paste0("OBIS_", row_number()),
    date = eventDate %>%
      stringr::str_extract("^\\d{4}-\\d{2}-\\d{2}") %>%  # extract first YYYY-MM-DD
      as.Date(),                                         # convert to Date
    lon = st_coordinates(obis_new)[,1],
    lat = st_coordinates(obis_new)[,2],
    species_aphia = aphia,
    species_name = "Lanice conchilega",
    presence_absence = 1
  )

# Standardise EMODnet PA columns
pa_ns_clean <- pa_ns_updated %>%
  dplyr::select(
    event_id,
    date,
    lon,
    lat,
    species_aphia,
    species_name,
    presence_absence
  )

# Combine EMODnet PA + new OBIS presences
pa_combined <- bind_rows(
  pa_ns_clean,
  obis_new_df
)

# Add pseudo-absences in the German EEZ (EEZ shapefile downloaded from MarineRegions)

eez_path <- "eez_v12.shp"
eez <- sf::st_read(eez_path)

# Inspect metadata
names(eez)

# Filter to Germany
eez_de <- eez %>% filter(grepl("Germany", SOVEREIGN1, ignore.case = TRUE))

# Reproject to WGS84 if needed
eez_de <- st_transform(eez_de, 4326)
eez_de <- st_make_valid(eez_de)
eez_de <- st_union(eez_de)

# Add pseudo-absences for Germany
library(purrr)
library(stringr)

# Parameters that can be tuned
buffer_km_default <- 2      # avoid labeling absence close to known presences
ratio_default     <- 1       # 1 pseudo-absence per presence
max_n_default     <- 500    # safety cap per species
seed_default      <- 42

#Helper: sample pseudo-absences for one species
make_pseudo_absences_for_species <- function(sp_row,
                                             pa_df,
                                             mask_poly,
                                             buffer_km = buffer_km_default,
                                             ratio = ratio_default,
                                             max_n = max_n_default,
                                             seed = seed_default) {
  set.seed(seed)
  
  sp_aphia <- sp_row$aphia_id
  sp_name  <- sp_row$species_name
  
  # 1) PRESENCES = EMODnet presences + OBIS/GBIF presences for this species
  #    (We assume pa_df already contains all sources and date parsed)
  pres_sp <- pa_df %>%
    filter(species_aphia == sp_aphia, presence_absence == 1) %>%
    filter(!is.na(lon), !is.na(lat)) %>%
    distinct(lon, lat, date, .keep_all = TRUE)
  
  n_pres <- nrow(pres_sp)
  if (n_pres == 0) {
    message("No presences for ", sp_name, "; skipping pseudo-absences.")
    return(NULL)
  }
  
  # 2) Decide how many pseudo-absences to draw
  n_pa <- min(max_n, ceiling(ratio * n_pres))
  
  # 3) Build a buffer around presences (meters) in a projected CRS
  #    Use ETRS89 / LAEA Europe (EPSG:3035) for robust meter buffering
  pres_pts <- st_as_sf(pres_sp, coords = c("lon","lat"), crs = 4326)
  pres_pts_3035 <- st_transform(pres_pts, 3035)
  mask_3035 <- st_transform(mask_poly, 3035)
  
  # We adaptively back off the buffer if area becomes empty
  buffer_seq_km <- c(buffer_km, 7.5, 5, 2.5, 1, 0.5, 0)  # km
  pa_pts <- NULL
  
  for (bk in buffer_seq_km) {
    buf_m <- bk * 1000
    if (bk > 0) {
      pres_buf <- st_union(st_buffer(pres_pts_3035, dist = buf_m))
      available <- suppressWarnings(st_difference(mask_3035, pres_buf))
    } else {
      # No buffer, just use the mask
      available <- mask_3035
    }
    
    # If available area is empty, try next (smaller) buffer
    if (is.null(available) || length(available) == 0 || all(st_is_empty(available))) {
      next
    }
    
    # 4) Sample pseudo-absence points
    # Note: st_sample on MULTIPOLYGON; 'random' is fine here
    # To avoid failures with very small areas, we try a slightly higher ask then thin if needed
    n_try <- max(n_pa, ceiling(1.2 * n_pa))
    samp  <- try(suppressWarnings(st_sample(available, size = n_try, type = "random")), silent = TRUE)
    if (inherits(samp, "try-error") || length(samp) == 0) next
    
    # Keep exactly n_pa
    samp <- samp[seq_len(min(length(samp), n_pa))]
    # Back-transform to WGS84
    samp_wgs84 <- st_transform(st_as_sf(samp), 4326)
    coords <- st_coordinates(samp_wgs84)
    
    # 5) Assign dates: sample from this species' presence dates (>=1990 in filtered dataset)
    sp_dates <- pres_sp %>% filter(!is.na(date)) %>% pull(date)
    if (length(sp_dates) == 0) {
      # Fallback: set to a fixed valid date in 1990+
      sampled_dates <- rep(as.Date("2000-06-15"), nrow(samp_wgs84))
    } else {
      sampled_dates <- sample(sp_dates, size = nrow(samp_wgs84), replace = TRUE)
    }
    
    # 6) Build pseudo-absence tibble
    pa_tbl <- tibble::tibble(
      event_id          = paste0("PSEUDO_DE_", str_replace_all(sp_name, " ", "_"), "_", seq_len(nrow(samp_wgs84))),
      date              = sampled_dates,
      lon               = coords[,1],
      lat               = coords[,2],
      species_aphia     = sp_aphia,
      species_name      = sp_name,
      presence_absence  = 0L,
      source            = paste0("PseudoAbsence (DE EEZ, buffer=", bk, "km)")
    )
    
    pa_pts <- pa_tbl
    message("Pseudo-absences for ", sp_name, ": n=", nrow(pa_pts), " (buffer=", bk, " km)")
    break
  }
  
  if (is.null(pa_pts)) {
    warning("Could not place pseudo-absences for ", sp_name, " — EEZ area too small for chosen buffer. Try reducing buffer_km.")
  }
  
  pa_pts
}

# -------------------------------
# Run pseudo-absence generation
# -------------------------------

pa_for_pa <- pa_combined

pseudo_abs_list <- purrr::pmap(sp_list, ~ make_pseudo_absences_for_species(
  sp_row = tibble::tibble(aphia_id = ..1, species_name = ..2),
  pa_df  = pa_for_pa,
  mask_poly = eez_de,
  buffer_km = 2,     # <-- tweak if needed
  ratio     = 1,      # <-- tweak (e.g., 2 or 3)
  max_n     = 500,
  seed      = 42
))

pseudo_abs_all <- pseudo_abs_list %>% dplyr::bind_rows()

# Combine with current table (which already contains EMODnet PA + OBIS presences)
pa_combined <- pa_combined %>%
  dplyr::bind_rows(pseudo_abs_all) %>%
  filter(!is.na(date) & date >= as.Date("1990-01-01"))

leaflet(subset(pa_combined, species_name == "Lanice conchilega")) %>%
  addTiles() %>%
  addCircleMarkers(
    color = ~ifelse(presence_absence == 1, "green", "red"),
    radius = 2, stroke = FALSE, fillOpacity = 0.8
  ) %>%
  addMeasure(
    position = "topleft",
    primaryLengthUnit = "meters",
    secondaryLengthUnit = "kilometers",
    primaryAreaUnit   = "sqmeters",
    secondaryAreaUnit = "hectares"
  )

pa_combined_sf <- st_as_sf(
  pa_combined,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)

# Save
st_write(
  pa_combined_sf,
  "./data/derived_data/pa_combined.shp",
  delete_dsn = TRUE
)
