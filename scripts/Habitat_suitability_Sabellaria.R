# ============================================================
# Sabellaria spinulosa – Habitat suitability modelling (RF)
# Spatial block cross-validation + final map
# ============================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# --- Libraries ------------------------------------------------
packages <- c(
  "sf","terra","ggplot2","patchwork","dplyr", "tidyr",
  "randomForest","pROC","ecospat"
)

new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(new_packages)) install.packages(new_packages)
invisible(lapply(packages, library, character.only = TRUE))

# --- Paths ----------------------------------------------------
output_dir       <- "Output/"
coastline_shp    <- "kustwgs.shp"
rasters_dir      <- "Data/Environment/"
pa_path          <- "pa_combined.shp"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- CRS ------------------------------------------------------
laea_crs <- "EPSG:3035"

# --- Coastline (for plotting only) ----------------------------
coastline <- st_read(coastline_shp, quiet = TRUE) |>
  st_transform(4258) |>
  st_make_valid() |>
  st_cast("MULTILINESTRING") |>
  st_simplify(dTolerance = 0.01)

# --- Presence / absence data ----------------------------------
pa_sf <- st_read(pa_path, quiet = TRUE) |>
  st_transform(4258) |>
  filter(!st_is_empty(geometry))

# Separate geometry and response
presabs <- pa_sf |>
  select(prsnc_b, geometry)

# --- Remove absences close to presences -----------------------
buffer_m <- 2000

presabs_laea <- st_transform(presabs, laea_crs)

pres_pts <- filter(presabs_laea, prsnc_b == 1L)
abs_pts  <- filter(presabs_laea, prsnc_b == 0L)

presence_buffer <- st_union(st_geometry(pres_pts)) |> st_buffer(buffer_m)

abs_pts <- abs_pts[!st_intersects(abs_pts, presence_buffer, sparse = FALSE)[,1],]

presabs <- bind_rows(pres_pts, abs_pts) |>
  st_transform(4258) |>
  mutate(presabs = as.integer(prsnc_b))

# --- Environmental rasters ------------------------------------
raster_stack_path <- file.path(output_dir, "raster_stack.tif")
tif_files <- list.files(rasters_dir, pattern = "\\.tif$", full.names = TRUE)

if (!file.exists(raster_stack_path)) {
  
  rasters <- lapply(tif_files, \(f) project(rast(f), laea_crs))
  template <- rasters[[1]]
  
  rasters <- lapply(rasters, \(r) resample(extend(r, template), template))
  rasters <- rast(rasters)
  
  writeRaster(rasters, raster_stack_path, overwrite = TRUE)
  
} else {
  rasters <- rast(raster_stack_path)
}

# --- Extract raster values at observation locations ------------
presabs <- st_transform(presabs, crs(rasters))
coords  <- st_coordinates(presabs)

data_RF <- cbind(
  presabs[, c("presabs", "geometry")],
  terra::extract(rasters, coords)
) |>
  na.omit()

# Add coordinates for spatial blocking
coords_laea <- st_coordinates(st_transform(data_RF, laea_crs))
data_RF$X <- coords_laea[,1]
data_RF$Y <- coords_laea[,2]

df <- st_drop_geometry(data_RF)

# --- Correlation analysis of environmental predictors ----------
selected_vars <- names(rasters)

# Subset predictor matrix only (no response, no coordinates)
pred_matrix <- df[, selected_vars]

# Compute correlation matrix
# Pearson is standard; switch to method = "spearman" if preferred
cor_mat <- cor(
  pred_matrix,
  use = "pairwise.complete.obs",
  method = "pearson"
)

# Inspect numeric correlation values
round(cor_mat, 2)

# Visualise correlations (very helpful for ecological judgement)
ggcorrplot::ggcorrplot(
  cor_mat,
  type = "lower",
  lab = TRUE,
  lab_size = 2.5,
  show.diag = TRUE,
  title = "Correlation matrix of environmental predictors"
)

# Identify highly correlated variable pairs
cor_threshold <- 0.9

high_cor_pairs <- which(
  abs(cor_mat) >= cor_threshold & abs(cor_mat) < 1,
  arr.ind = TRUE
)

high_cor_pairs <- data.frame(
  variable_1 = rownames(cor_mat)[high_cor_pairs[,1]],
  variable_2 = colnames(cor_mat)[high_cor_pairs[,2]],
  r = cor_mat[high_cor_pairs]
)

# Remove mirrored duplicates (A–B vs B–A)
high_cor_pairs <- high_cor_pairs[
  high_cor_pairs$variable_1 < high_cor_pairs$variable_2, ]

high_cor_pairs

# --- Predictor selection: allow up to 10% missing at presences --

# Number of presence records
n_pres <- sum(df$presabs == 1L)

vars_ok <- df |>
  filter(presabs == 1L) |>
  summarise(
    across(
      all_of(selected_vars),
      ~ sum(is.na(.)) / n_pres
    )
  ) |>
  pivot_longer(
    everything(),
    names_to  = "variable",
    values_to = "na_frac"
  ) |>
  filter(na_frac <= 0.10) |>
  pull(variable)

vars_ok

# --- Remove predictors dropped due to high collinearity ----------

vars_to_drop <- c(
  "bottom_chlfa_max",
  "bpi5",
  "bpi30",
  "bpi75",
  "DCSM_mean_stress"
)

vars_ok <- setdiff(vars_ok, vars_to_drop)

vars_ok

# --- Spatial block cross-validation ---------------------------
rf_train <- df |>
  select(presabs, all_of(vars_ok), X, Y) |>
  drop_na()

block_size <- 500000   # 500 km blocks

rf_train <- rf_train |>
  mutate(
    block_x  = floor(X / block_size),
    block_y  = floor(Y / block_size),
    block_id = interaction(block_x, block_y, drop = TRUE),
    presabs  = factor(presabs, levels = c(0,1))
  )

set.seed(42)

cv_results <- list()

for (b in levels(rf_train$block_id)) {
  
  test_id  <- rf_train$block_id == b
  if (sum(test_id) < 20) next
  
  train_dat <- rf_train[!test_id,]
  test_dat  <- rf_train[test_id,]
  
  rf_cv <- randomForest(
    presabs ~ .,
    data = train_dat |> select(presabs, all_of(vars_ok)),
    ntree = 500
  )
  
  pred <- predict(
    rf_cv,
    test_dat |> select(all_of(vars_ok)),
    type = "prob"
  )[,"1"]
  
  cv_results[[as.character(b)]] <- data.frame(
    obs  = test_dat$presabs,
    pred = pred
  )
}

cv_df <- do.call(rbind, cv_results)

# --- CV metrics -----------------------------------------------
auc_spatial <- roc(cv_df$obs, cv_df$pred)$auc
auc_spatial

roc_obj <- roc(cv_df$obs, cv_df$pred)
thresh  <- as.numeric(coords(roc_obj, "best",
                             ret = "threshold",
                             best.method = "youden")$threshold)

pred_bin <- ifelse(cv_df$pred >= thresh, 1, 0)

sens <- sum(pred_bin == 1 & cv_df$obs == 1) /
  sum(cv_df$obs == 1)

spec <- sum(pred_bin == 0 & cv_df$obs == 0) /
  sum(cv_df$obs == 0)

tss_spatial <- sens + spec - 1
tss_spatial

# --- Final model ----------------------------------------------
rf_final <- rf_train |>
  select(presabs, all_of(vars_ok))

rf_model <- randomForest(
  presabs ~ .,
  data = rf_final,
  ntree = 500,
  importance = TRUE
)

# --- Check variable importance --------------------------------
importance(rf_model)
# --- Prediction map -------------------------------------------
rf_prob_fun <- function(model, data) {
  predict(model, data, type = "prob")[,"1"]
}

pred_raster_rf <- predict(
  rasters,
  rf_model,
  fun = rf_prob_fun,
  na.rm = TRUE
)

writeRaster(
  pred_raster_rf,
  file.path(output_dir, "Sabellaria_spinulosa_probability_RF.tif"),
  overwrite = TRUE
)
