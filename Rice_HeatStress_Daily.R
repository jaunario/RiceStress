# Title: Seasonal Rice Heat Stress Assessment Engine (Daily Temperature Rasters)
# Author: Jorrel Khalil Aunario, Avy
# Description: Dedicated engine to assess seasonal rice heat stress incidence (daytime heat / nighttime heat)
#              using daily weather rasters (AgERA5). Aligned with the published methodology
#              for reproductive and flowering phases.
# Environment & Package Versions:
#   - R     : 4.6.1 (2026-06-24)
#   - terra : 1.9.34

suppressPackageStartupMessages({
  library(terra)
})

if (!exists("dateFromDoy", mode = "function")) {
  dateFromDoy <- function(doy, year) as.Date(paste(year, doy), "%Y %j")
}
if (!exists("monthFromDate", mode = "function")) {
  monthFromDate <- function(dt) as.integer(format(as.Date(dt), "%m"))
}

# ==============================================================================
# 1. USER CONFIGURATION
# ==============================================================================

# Select Heat Stress Type: "daytime_heat" (DTH) or "nighttime_heat" (NTH)
HEAT_STRESS_TYPE <- "daytime_heat"

# Target Crop Stage/Phase: "reproductive", "flowering", or "custom"
CROP_PHASE <- "reproductive"

# Custom stage offsets relative to harvest date (used only when CROP_PHASE == "custom")
CUSTOM_STAGE <- list(
  start_offset = -56,
  end_offset   = -26
)

# Daily Weather Rasters Source Directory (output of Download_Standardize_AgERA5.R)
SRC_DIRECTORY <- "./agera5_daily"
SRC_FILETYPE  <- "tif"

# Time Series & Seasons to Analyze
ANALYSIS_YEARS   <- 1993:2023
ANALYSIS_SEASONS <- 1:3
SPATIAL_RES      <- 0.1 # 6 arcmin / 0.1 degree

# Input Spatial Files (provide path to Rice Crop Calendar shapefile and optional Rice Mask raster)
# If you do not have a local Rice Crop Calendar shapefile, download RiceAtlas:
#   Dataset (Harvard Dataverse): https://doi.org/10.7910/DVN/JE6R2R
#   Reference: Laborte et al. (2017), Scientific Data 4, 170074 (https://doi.org/10.1038/sdata.2017.74)
CROPCAL_FILE  <- "path/to/rice_crop_calendar.shp"
RICEMASK_FILE <- NULL # Optional: "path/to/rice_mask.tif"

# Optional AOI Filter (set to NULL for full Crop Calendar extent)
# Example: AOI_FILTER <- list(field = "COUNTRY", value = "Philippines")
AOI_FILTER <- NULL

OUTPUT_DIR <- paste0("./out/", ifelse(grepl("night", HEAT_STRESS_TYPE, ignore.case = TRUE), "NTH", "DTH"), "_AgERA5_Daily/")

# ==============================================================================
# 2. HEAT STRESS PRESET & STAGE DEFINITION
# ==============================================================================

# Heat stress definition presets:
# Daytime heat:   tmax > 35°C, daycount threshold = 10
# Nighttime heat: tmin > 25°C, daycount threshold = 15
get_heat_preset <- function(stress_type) {
  if (grepl("night", stress_type, ignore.case = TRUE)) {
    return(list(name = "NTH", variable = "tmin", threshold = 25, comparator = ">", daycount = 15))
  } else {
    return(list(name = "DTH", variable = "tmax", threshold = 35, comparator = ">", daycount = 10))
  }
}

resolve_stage_dates <- function(date_plant, date_harv, phase_type, custom_cfg = NULL) {
  phase_key <- tolower(phase_type)
  if (phase_key %in% c("reproductive", "repr")) {
    # Panicle initiation to end of flowering: -65 + 14 - 5 to -30 + 4
    st <- date_harv - 56
    en <- date_harv - 26
  } else if (phase_key %in% c("flowering", "flwr", "anthesis")) {
    # Anthesis peak flowering window: -30 - 4 to -30 + 4
    st <- date_harv - 34
    en <- date_harv - 26
  } else if (phase_key == "custom" && !is.null(custom_cfg)) {
    st <- date_harv + custom_cfg$start_offset
    en <- date_harv + custom_cfg$end_offset
  } else {
    stop("Unknown CROP_PHASE: '", phase_type, "'. Must be 'reproductive', 'flowering', or 'custom'.")
  }
  return(list(st = st, en = en))
}

RUN_SETTING <- get_heat_preset(HEAT_STRESS_TYPE)

message("========================================================")
message("  Seasonal Rice Heat Stress Assessment Engine")
message("  Stress Type   : ", RUN_SETTING$name, " (", RUN_SETTING$variable, " ", RUN_SETTING$comparator, " ", RUN_SETTING$threshold, " C)")
message("  Day Threshold : >= ", RUN_SETTING$daycount, " days")
message("  Crop Stage    : ", toupper(CROP_PHASE))
message("  Years         : ", min(ANALYSIS_YEARS), " - ", max(ANALYSIS_YEARS))
message("  Output Dir    : ", OUTPUT_DIR)
message("========================================================")

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if (!dir.exists(file.path(OUTPUT_DIR, "Masked"))) dir.create(file.path(OUTPUT_DIR, "Masked"), recursive = TRUE)

# ==============================================================================
# 3. BASE RASTERS & DATA INVENTORY
# ==============================================================================

rst.base <- rast()
res(rst.base) <- SPATIAL_RES

if (!file.exists(CROPCAL_FILE)) {
  stop(
    "Crop calendar not found at: ", CROPCAL_FILE, "\n",
    "Download RiceAtlas (Laborte et al., 2017) from Harvard Dataverse: https://doi.org/10.7910/DVN/JE6R2R"
  )
}
vec.cropcal <- vect(CROPCAL_FILE)

if (!is.null(AOI_FILTER) && AOI_FILTER$field %in% names(vec.cropcal)) {
  vec.cropcal <- vec.cropcal[vec.cropcal[[AOI_FILTER$field]][, 1] == AOI_FILTER$value, ]
  rst.base <- crop(rst.base, ext(vec.cropcal), snap = "out")
  message("Filtered crop calendar to AOI: ", AOI_FILTER$field, " == ", AOI_FILTER$value, " (", length(vec.cropcal), " units)")
}
dat.cropcal <- as.data.frame(vec.cropcal)

# Optional rice mask for final products
rst.rice <- NULL
if (!is.null(RICEMASK_FILE) && nzchar(RICEMASK_FILE) && file.exists(RICEMASK_FILE)) {
  rst.rice <- rast(RICEMASK_FILE)
  if (!compareGeom(rst.rice, rst.base, stopOnError = FALSE)) {
    rst.rice <- resample(rst.rice, rst.base, method = "near")
  }
}

# Source Daily Weather Files Inventory
files.src <- dir(SRC_DIRECTORY, pattern = paste0(SRC_FILETYPE, "$"), recursive = TRUE, full.names = TRUE)
if (length(files.src) == 0) {
  warning("No weather raster files found in: ", SRC_DIRECTORY, ". Verify path or drive mount.")
}

PARAM_DELIM <- "\\."
PARAM_NAMES <- c("setname", "spat_res", "variable", "year", "month", "dom", "fileext")
NUM_FIELDS  <- 4:6

info.fromnames <- sapply(basename(files.src), FUN = strsplit, split = PARAM_DELIM, USE.NAMES = FALSE)
info.fromnames <- try(do.call(rbind, info.fromnames), silent = TRUE)
if (!inherits(info.fromnames, "try-error") && length(info.fromnames) > 0) {
  colnames(info.fromnames) <- PARAM_NAMES
  dat.inventory <- as.data.frame(info.fromnames)
  dat.inventory[, NUM_FIELDS] <- as.data.frame(lapply(dat.inventory[, NUM_FIELDS], as.numeric, SIMPLIFY = FALSE))
  dat.inventory$date <- as.Date(paste(dat.inventory$year, dat.inventory$month, dat.inventory$dom), "%Y %m %d")
  dat.inventory$filename <- files.src
} else {
  dat.inventory <- data.frame(filename = character(0), variable = character(0), date = as.Date(character(0)))
}

# ==============================================================================
# 4. ANNUAL HEAT STRESS COMPUTATION
# ==============================================================================

for (season in ANALYSIS_SEASONS) {
  txt.pl <- paste0("PLANT_PK", season)
  txt.hv <- paste0("HARV_PK", season)

  rice.units <- which(dat.cropcal[, txt.pl] > 0 & dat.cropcal[, txt.hv] > 0)
  if (length(rice.units) == 0) {
    message("Season ", season, ": No rice calendar units found. Skipping.")
    next
  }

  vec.thisseason <- vec.cropcal[rice.units, ]
  dat.thisseason <- dat.cropcal[rice.units, ]
  if (exists("rst.yearly")) rm(rst.yearly)

  for (y in ANALYSIS_YEARS) {
    message(y, "-S", season, ": deriving target period (", CROP_PHASE, ")")

    dat.thisseason$date.plant <- dateFromDoy(dat.thisseason[, txt.pl], y)
    dat.thisseason$date.harv <- dateFromDoy(
      dat.thisseason[, txt.hv],
      ifelse(dat.thisseason[, txt.hv] < dat.thisseason[, txt.pl], y + 1, y)
    )

    stage_dates <- resolve_stage_dates(
      date_plant = dat.thisseason$date.plant,
      date_harv  = dat.thisseason$date.harv,
      phase_type = CROP_PHASE,
      custom_cfg = CUSTOM_STAGE
    )

    dat.thisseason$date.targetperiod.st <- stage_dates$st
    dat.thisseason$date.targetperiod.en <- stage_dates$en
    dat.thisseason$month.plant.st       <- monthFromDate(dat.thisseason$date.plant)
    dat.thisseason$month.target.st      <- monthFromDate(dat.thisseason$date.targetperiod.st)

    if (exists("rst.monthly")) rm(rst.monthly)

    for (m in 1:12) {
      for (mm in 1:12) {
        pl.thismonth <- which(dat.thisseason$month.plant.st == m & dat.thisseason$month.target.st == mm)
        if (length(pl.thismonth) == 0) next

        message(y, "-S", season, ": AOIs planting ", month.abb[m], " target phase ", month.abb[mm])
        vec.thismonth <- vec.thisseason[pl.thismonth, ]
        vec.thismonth$reprst <- as.numeric(dat.thisseason$date.targetperiod.st[pl.thismonth])
        vec.thismonth$repren <- as.numeric(dat.thisseason$date.targetperiod.en[pl.thismonth])
        reqdate.st <- min(dat.thisseason$date.targetperiod.st[pl.thismonth])
        reqdate.en <- max(dat.thisseason$date.targetperiod.en[pl.thismonth])

        rst.reprst <- rasterize(vec.thismonth, rst.base, field = "reprst", touches = TRUE)
        rst.reprst <- crop(rst.reprst, vec.thismonth)
        rst.repren <- rasterize(vec.thismonth, rst.base, field = "repren", touches = TRUE)
        rst.repren <- crop(rst.repren, vec.thismonth)

        files.weather <- subset(
          dat.inventory,
          date >= reqdate.st & date <= reqdate.en & variable == RUN_SETTING$variable,
          select = c("date", "filename")
        )
        if (nrow(files.weather) == 0) next

        message(y, "-S", season, ": Loading ", nrow(files.weather), " ", RUN_SETTING$variable, " files")
        rst.wth <- rast(files.weather$filename)
        rst.wth <- crop(rst.wth, vec.thismonth)
        rst.wth <- mask(rst.wth, rst.reprst, maskvalue = NA)
        rst.wth <- rst.wth / 100 # convert to deg C

        # Temporal masking to pixel-specific window
        rst.date <- rast(rst.wth)
        values(rst.date) <- sapply(as.numeric(files.weather$date), rep, ncell(rst.date))
        rst.wth <- ifel(rst.date < rst.reprst | rst.date > rst.repren, NA, rst.wth)

        # Threshold checking and day count accumulation
        rst.vthres <- app(rst.wth, fun = RUN_SETTING$comparator, RUN_SETTING$threshold)
        rst.vthres_count <- sum(rst.vthres, na.rm = TRUE)

        if (!exists("rst.monthly")) {
          rst.monthly <- rst.vthres_count
        } else {
          rst.monthly <- merge(rst.monthly, rst.vthres_count)
        }

        rm(rst.vthres_count, rst.vthres, rst.wth, rst.date, rst.repren, rst.reprst, vec.thismonth)
        tmpFiles(orphan = TRUE, remove = TRUE)
        gc(reset = TRUE, verbose = FALSE)
      }
    }

    if (!exists("rst.monthly")) next
    names(rst.monthly) <- paste0("COUNT.VTHRES.Y", y)
    rst.monthly.resampled <- resample(rst.monthly, rst.base, method = "near")

    if (!exists("rst.yearly")) {
      rst.yearly <- rst.monthly.resampled
    } else {
      rst.yearly <- c(rst.yearly, rst.monthly.resampled)
    }
  }

  if (!exists("rst.yearly")) next

  # Save Annual Exceedance Day Count Rasters (INT2S)
  actual_years <- ANALYSIS_YEARS[seq_len(nlyr(rst.yearly))]
  names(rst.yearly) <- paste0("XTHRES_COUNT_Y", actual_years)
  fn_yearly <- paste0(OUTPUT_DIR, RUN_SETTING$name, "_SEASON", season, "_", names(rst.yearly), ".tif")
  writeRaster(rst.yearly, filename = fn_yearly, datatype = "INT2S", overwrite = TRUE)

  # Save Cumulative Stressed Years Count (1993-YYYY)
  rst.stressyrs <- sum(rst.yearly >= RUN_SETTING$daycount, na.rm = TRUE)
  fn_stressyrs <- paste0(OUTPUT_DIR, RUN_SETTING$name, "_SEASON", season, "_AgERA5_STRESSYRS_COUNT_", min(actual_years), "-", max(actual_years), ".tif")
  writeRaster(rst.stressyrs, filename = fn_stressyrs, datatype = "INT2S", overwrite = TRUE)

  # Save Masked Versions
  if (!is.null(rst.rice)) {
    rst.stressyrs.masked <- mask(rst.stressyrs, rst.rice, maskvalue = NA)
    fn_masked <- paste0(OUTPUT_DIR, "Masked/", RUN_SETTING$name, "_SEASON", season, "_AgERA5_STRESSYRS_COUNT_", min(actual_years), "-", max(actual_years), "_masked.tif")
    writeRaster(rst.stressyrs.masked, filename = fn_masked, datatype = "INT2S", overwrite = TRUE)
  }

  message("Completed Season ", season, ". Outputs written to: ", OUTPUT_DIR)
}

message("Rice Heat Stress Daily Processing Complete.")
