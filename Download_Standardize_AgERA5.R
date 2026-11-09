# Title: Generic AgERA5 Daily Weather Downloader and Standardizer (CDS API)
# Author: Jorrel Khalil Aunario, Avy
# Description: Generic script for downloading AgERA5 daily weather archives from the
#              Copernicus Climate Data Store (CDS) using a user-supplied CDS Personal
#              Access Token, supporting Area of Interest (AOI) selection via a
#              terra::SpatExtent, a terra::SpatVector, or a loaded/filtered Rice Crop
#              Calendar, with time-period validation and CDS server limitation warnings.
# Environment & Package Versions:
#   - R         : 4.6.1 (2026-06-24)
#   - terra     : 1.9.34
#   - Python    : 3.14.8
#   - cdsapi    : 0.7.7

suppressPackageStartupMessages({
  library(terra)
})

# ==============================================================================
# 1. USER CONFIGURATION
# ==============================================================================

# --- A. Copernicus CDS Authentication ---
# Leave CDS_API_TOKEN as NULL to be prompted interactively for your personal token.
# Obtain your Personal Access Token from: https://cds.climate.copernicus.eu/profile
CDS_API_TOKEN <- NULL
CDS_API_URL   <- "https://cds.climate.copernicus.eu/api"

# --- B. Area of Interest (AOI) Selection ---
# Choose ONE of the following three AOI methods:
#   1. SpatExtent object (or numeric c(xmin, xmax, ymin, ymax)):
#      e.g., AOI_EXTENT <- ext(120.8, 121.8, 13.8, 14.8)
#   2. SpatVector object (or path to a vector file .shp/.gpkg/.geojson):
#      e.g., AOI_VECTOR <- vect("path/to/study_area.shp")
#   3. Rice Crop Calendar (loaded as SpatVector or path to shapefile, optionally filtered):
#      e.g., CROPCAL_DATA   <- vect("path/to/rice_crop_calendar.shp")
#            CROPCAL_FILTER <- list(COUNTRY = "Philippines", SUB_REGION = "Laguna")
#      If you do not have a local Rice Crop Calendar shapefile, download RiceAtlas:
#        Dataset (Harvard Dataverse): https://doi.org/10.7910/DVN/JE6R2R
#        Reference: Laborte et al. (2017), Scientific Data 4, 170074 (https://doi.org/10.1038/sdata.2017.74)
AOI_EXTENT     <- NULL
AOI_VECTOR     <- NULL
CROPCAL_DATA   <- NULL
CROPCAL_FILTER <- NULL # Optional: named list of attribute filters, e.g. list(SUB_REGION = "Laguna")
ACTIVE_SEASON  <- NULL # Optional: 1, 2, or 3 to require PLANT_PKn > 0 & HARV_PKn > 0
AOI_BUFFER_DEG <- 0.2  # Padding in degrees around resolved AOI extent

# --- C. Time Period & Batching ---
DOWNLOAD_YEARS  <- 2024
DOWNLOAD_MONTHS <- 1:12
# "monthly" submits 1 request/month (default to stay within CDS per-request cost limits);
# "yearly" submits 1 request/year (only for single-statistic or short month subsets).
BATCH_MODE      <- "monthly"

# Thresholds for triggering CDS server limitation warnings
WARN_MAX_YEARS        <- 2
WARN_MAX_REQUESTS     <- 12
WARN_LARGE_AREA_SQDEG <- 2500
CONFIRM_ON_WARNING    <- TRUE

# --- D. Target Variables & Output Directories ---
AGERA5_VARIABLE   <- "2m_temperature"
AGERA5_STATISTICS <- c("24_hour_minimum", "24_hour_maximum")
AGERA5_VERSION    <- "2_0"

RAW_ZIP_DIR   <- "./agera5_raw"
STANDARD_DIR  <- "./agera5_daily"
SPATIAL_RES   <- 0.1 # 0.1 deg (6 arcmin)
PYTHON_BIN    <- "python3"
DOWNLOADER_PY <- "./download_agera5.py"

# ==============================================================================
# 2. HELPER FUNCTIONS: AUTHENTICATION, AOI RESOLUTION & LIMIT CHECKS
# ==============================================================================

prompt_cds_token <- function(token = NULL) {
  if (!is.null(token) && nzchar(trimws(token))) {
    return(trimws(token))
  }
  env_token <- Sys.getenv("CDS_API_TOKEN", unset = "")
  if (nzchar(trimws(env_token))) {
    return(trimws(env_token))
  }
  default_rc <- path.expand("~/.cdsapirc")
  if (interactive()) {
    message("\n=== Copernicus CDS Authentication ===")
    message("Obtain your Personal Access Token at: https://cds.climate.copernicus.eu/profile")
    prompt_str <- if (file.exists(default_rc)) {
      "Enter your CDS Personal Access Token (or press Enter to use ~/.cdsapirc): "
    } else {
      "Enter your CDS Personal Access Token: "
    }
    entered <- trimws(readline(prompt = prompt_str))
    if (nzchar(entered)) {
      return(entered)
    }
  }
  if (file.exists(default_rc)) {
    return("")
  }
  stop(
    "CDS Personal Access Token is required. Set CDS_API_TOKEN in the script, ",
    "pass it via environment variable CDS_API_TOKEN, or run interactively."
  )
}

resolve_aoi_extent <- function(aoi_extent = NULL,
                               aoi_vector = NULL,
                               cropcal = NULL,
                               cropcal_filter = NULL,
                               active_season = NULL,
                               buffer_deg = 0.2) {
  raw_ext <- NULL
  aoi_label <- NULL
  region_tag <- "aoi"

  # 1. SpatExtent or numeric c(xmin, xmax, ymin, ymax)
  if (!is.null(aoi_extent)) {
    if (inherits(aoi_extent, "SpatExtent")) {
      raw_ext <- aoi_extent
    } else if (is.numeric(aoi_extent) && length(aoi_extent) == 4) {
      raw_ext <- ext(aoi_extent[1], aoi_extent[2], aoi_extent[3], aoi_extent[4])
    } else {
      stop("AOI_EXTENT must be a terra::SpatExtent object or numeric c(xmin, xmax, ymin, ymax).")
    }
    aoi_label <- sprintf(
      "SpatExtent (xmin=%.4f, xmax=%.4f, ymin=%.4f, ymax=%.4f)",
      xmin(raw_ext), xmax(raw_ext), ymin(raw_ext), ymax(raw_ext)
    )
    region_tag <- sprintf("bbox_%.1f_%.1f_%.1f_%.1f", xmin(raw_ext), xmax(raw_ext), ymin(raw_ext), ymax(raw_ext))

  # 2. SpatVector object or vector file path
  } else if (!is.null(aoi_vector)) {
    vec_obj <- if (inherits(aoi_vector, "SpatVector")) {
      aoi_vector
    } else if (is.character(aoi_vector) && file.exists(aoi_vector)) {
      vect(aoi_vector)
    } else {
      stop("AOI_VECTOR must be a terra::SpatVector object or valid vector file path.")
    }
    if (nrow(vec_obj) == 0) stop("Provided SpatVector has 0 features.")
    raw_ext <- ext(vec_obj)
    aoi_label <- sprintf("SpatVector (%d feature(s))", nrow(vec_obj))
    region_tag <- if (is.character(aoi_vector)) tools::file_path_sans_ext(basename(aoi_vector)) else "spatvector"

  # 3. Rice Crop Calendar (SpatVector or shapefile path) with optional filtering
  } else if (!is.null(cropcal)) {
    vec_cal <- if (inherits(cropcal, "SpatVector")) {
      cropcal
    } else if (is.character(cropcal) && file.exists(cropcal)) {
      vect(cropcal)
    } else {
      stop(
        "Crop calendar not found: ", as.character(cropcal), "\n",
        "Download RiceAtlas (Laborte et al., 2017) from Harvard Dataverse: https://doi.org/10.7910/DVN/JE6R2R"
      )
    }

    if (!is.null(cropcal_filter) && length(cropcal_filter) > 0) {
      # Support both list(field = "SUB_REGION", value = "Laguna") and named list(SUB_REGION = "Laguna")
      if (all(c("field", "value") %in% names(cropcal_filter)) && length(cropcal_filter) == 2) {
        fld <- cropcal_filter$field
        val <- cropcal_filter$value
        if (!fld %in% names(vec_cal)) {
          stop("Filter field '", fld, "' not found in crop calendar. Available: ", paste(names(vec_cal), collapse = ", "))
        }
        vec_cal <- vec_cal[tolower(trimws(vec_cal[[fld]][, 1])) %in% tolower(trimws(val)), ]
        region_tag <- paste(val, collapse = "_")
      } else {
        for (fld in names(cropcal_filter)) {
          if (!fld %in% names(vec_cal)) {
            stop("Filter field '", fld, "' not found in crop calendar. Available: ", paste(names(vec_cal), collapse = ", "))
          }
          val <- cropcal_filter[[fld]]
          vec_cal <- vec_cal[tolower(trimws(vec_cal[[fld]][, 1])) %in% tolower(trimws(val)), ]
        }
        region_tag <- paste(unlist(cropcal_filter), collapse = "_")
      }
    }

    if (!is.null(active_season)) {
      pl_col <- paste0("PLANT_PK", active_season)
      hv_col <- paste0("HARV_PK", active_season)
      if (all(c(pl_col, hv_col) %in% names(vec_cal))) {
        keep_idx <- which(vec_cal[[pl_col]][, 1] > 0 & vec_cal[[hv_col]][, 1] > 0)
        vec_cal <- vec_cal[keep_idx, ]
      }
    }

    if (nrow(vec_cal) == 0) {
      stop("No Rice Crop Calendar features matched the specified filter criteria.")
    }

    raw_ext <- ext(vec_cal)
    aoi_label <- sprintf("Rice Crop Calendar (%d matched feature(s))", nrow(vec_cal))
  } else {
    stop("No Area of Interest configured. Set AOI_EXTENT, AOI_VECTOR, or CROPCAL_DATA.")
  }

  buffered_ext <- ext(
    max(-180, round(xmin(raw_ext) - buffer_deg, 2)),
    min(180,  round(xmax(raw_ext) + buffer_deg, 2)),
    max(-90,  round(ymin(raw_ext) - buffer_deg, 2)),
    min(90,   round(ymax(raw_ext) + buffer_deg, 2))
  )

  # CDS format: [North, West, South, East]
  cds_area <- c(
    North = ymax(buffered_ext),
    West  = xmin(buffered_ext),
    South = ymin(buffered_ext),
    East  = xmax(buffered_ext)
  )

  return(list(
    raw_extent      = raw_ext,
    buffered_extent = buffered_ext,
    cds_area        = cds_area,
    label           = aoi_label,
    tag             = gsub("[^A-Za-z0-9_-]", "_", tolower(region_tag))
  ))
}

check_cds_request_limits <- function(years,
                                     months,
                                     statistics,
                                     aoi_extent,
                                     batch_mode = "yearly",
                                     warn_max_years = 5,
                                     warn_max_requests = 10,
                                     warn_large_area = 2500,
                                     confirm = TRUE) {
  n_years  <- length(years)
  n_months <- length(months)
  total_months <- n_years * n_months
  approx_days  <- round(total_months * 30.4375)
  n_requests   <- if (batch_mode == "yearly") n_years else total_months
  area_sq_deg  <- (xmax(aoi_extent) - xmin(aoi_extent)) * (ymax(aoi_extent) - ymin(aoi_extent))
  fields_per_req <- if (batch_mode == "yearly") round(n_months * 30.4375) * length(statistics) else 31 * length(statistics)

  msgs <- character(0)

  if (n_years > warn_max_years) {
    msgs <- c(msgs, sprintf(
      paste0(
        "Long time period selected (%d years: %d-%d, ~%d days). ",
        "Copernicus CDS enforces strict queue concurrency and rate limits per user account; ",
        "multi-year runs may experience extended queue wait times or rate throttling."
      ),
      n_years, min(years), max(years), approx_days
    ))
  }

  if (n_requests > warn_max_requests) {
    msgs <- c(msgs, sprintf(
      paste0(
        "High API request count (%d requests in '%s' mode). ",
        "Consider downloading in smaller year batches (e.g., 5-10 years per run) ",
        "to avoid CDS queue saturation."
      ),
      n_requests, batch_mode
    ))
  }

  if (area_sq_deg > warn_large_area && batch_mode == "yearly") {
    msgs <- c(msgs, sprintf(
      paste0(
        "Large spatial extent (%.1f sq. deg) combined with 'yearly' batching (%d fields/request). ",
        "CDS may reject large payloads; set BATCH_MODE <- 'monthly' for continental or global extents."
      ),
      area_sq_deg, fields_per_req
    ))
  }

  if (length(msgs) > 0) {
    for (m in msgs) {
      warning("[CDS SERVER LIMITATION] ", m, call. = FALSE, immediate. = TRUE)
    }
    if (interactive() && confirm) {
      ans <- tolower(trimws(readline(
        prompt = sprintf("Proceed with submitting up to %d CDS API requests? [y/N]: ", n_requests)
      )))
      if (!ans %in% c("y", "yes")) {
        stop("Download cancelled by user.")
      }
    }
  }

  invisible(n_requests)
}

# ==============================================================================
# 3. STANDARDIZE NETCDF ARCHIVES TO GEOTIFF
# ==============================================================================

standardize_agera5_zips <- function(zip_files, out_dir, aoi_extent = NULL, spatial_res = 0.1) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  temp_extract_dir <- tempfile(pattern = "agera5_nc_")
  dir.create(temp_extract_dir, showWarnings = FALSE)
  on.exit(unlink(temp_extract_dir, recursive = TRUE), add = TRUE)

  # Base grid aligned to origin (0, 0) to match Rice_HeatStress_Daily.R / Generic_RiceStress_Daily.R
  rst_base <- rast()
  res(rst_base) <- spatial_res
  if (!is.null(aoi_extent)) {
    rst_base <- crop(rst_base, aoi_extent, snap = "out")
  }

  for (zf in zip_files) {
    if (!file.exists(zf)) next
    files_in_zip <- unzip(zf, list = TRUE)$Name
    nc_in_zip <- files_in_zip[grepl("\\.nc$", files_in_zip)]
    if (length(nc_in_zip) == 0) next

    # Check if all daily GeoTIFFs from this archive already exist
    needs_extract <- FALSE
    for (f_nc in nc_in_zip) {
      var_short <- if (grepl("Max-24h", f_nc, ignore.case = TRUE)) {
        "tmax"
      } else if (grepl("Min-24h", f_nc, ignore.case = TRUE)) {
        "tmin"
      } else {
        "var"
      }
      dt_match <- regmatches(f_nc, regexpr("[0-9]{8}", f_nc))
      if (length(dt_match) == 1) {
        formatted_date <- format(as.Date(dt_match, "%Y%m%d"), "%Y.%m.%d")
        out_tif <- file.path(out_dir, sprintf("AGERA5.6arcmin.%s.%s.tif", var_short, formatted_date))
        if (!file.exists(out_tif)) {
          needs_extract <- TRUE
          break
        }
      } else {
        needs_extract <- TRUE
        break
      }
    }

    if (!needs_extract) next

    unzip(zf, exdir = temp_extract_dir)
    extracted_nc <- dir(temp_extract_dir, pattern = "\\.nc$", full.names = TRUE)

    for (nc in extracted_nc) {
      bn <- basename(nc)
      var_short <- if (grepl("Max-24h", bn, ignore.case = TRUE)) {
        "tmax"
      } else if (grepl("Min-24h", bn, ignore.case = TRUE)) {
        "tmin"
      } else {
        "var"
      }

      dt_match <- regmatches(bn, regexpr("[0-9]{8}", bn))
      if (length(dt_match) == 1) {
        formatted_date <- format(as.Date(dt_match, "%Y%m%d"), "%Y.%m.%d")
      } else {
        r_tmp <- rast(nc)
        formatted_date <- format(time(r_tmp)[1], "%Y.%m.%d")
      }

      target_tif <- file.path(out_dir, sprintf("AGERA5.6arcmin.%s.%s.tif", var_short, formatted_date))
      if (!file.exists(target_tif)) {
        r <- try(rast(nc), silent = TRUE)
        if (!inherits(r, "try-error")) {
          target_grid <- if (!is.null(aoi_extent)) rst_base else crop(rst_base, ext(r[[1]]), snap = "out")
          r <- resample(r[[1]], target_grid, method = "near")

          vals <- values(r, mat = FALSE)
          first_val <- na.omit(vals)[1]
          if (!is.na(first_val) && first_val > 150) {
            r <- (r - 273.15) * 100
          } else {
            r <- r * 100
          }
          r <- round(r, 0)
          writeRaster(
            r,
            filename = target_tif,
            wopt = list(datatype = "INT2S", gdal = "COMPRESS=LZW"),
            overwrite = TRUE
          )
        }
      }
      unlink(nc)
    }
  }
}

# ==============================================================================
# 4. MAIN EXECUTION PIPELINE
# ==============================================================================

run_agera5_pipeline <- function(aoi_extent = AOI_EXTENT,
                                aoi_vector = AOI_VECTOR,
                                cropcal = CROPCAL_DATA,
                                cropcal_filter = CROPCAL_FILTER,
                                active_season = ACTIVE_SEASON,
                                buffer_deg = AOI_BUFFER_DEG,
                                years = DOWNLOAD_YEARS,
                                months = DOWNLOAD_MONTHS,
                                batch_mode = BATCH_MODE,
                                variable = AGERA5_VARIABLE,
                                statistics = AGERA5_STATISTICS,
                                version = AGERA5_VERSION,
                                raw_zip_dir = RAW_ZIP_DIR,
                                standard_dir = STANDARD_DIR,
                                spatial_res = SPATIAL_RES,
                                cds_token = CDS_API_TOKEN,
                                confirm_on_warning = CONFIRM_ON_WARNING) {
  # 1. Resolve Area of Interest
  aoi <- resolve_aoi_extent(
    aoi_extent     = aoi_extent,
    aoi_vector     = aoi_vector,
    cropcal        = cropcal,
    cropcal_filter = cropcal_filter,
    active_season  = active_season,
    buffer_deg     = buffer_deg
  )

  message("========================================================")
  message("  Generic AgERA5 Daily Downloader & Standardizer")
  message("  AOI Source  : ", aoi$label)
  message(sprintf(
    "  SpatExtent  : xmin=%.2f, xmax=%.2f, ymin=%.2f, ymax=%.2f",
    xmin(aoi$buffered_extent), xmax(aoi$buffered_extent),
    ymin(aoi$buffered_extent), ymax(aoi$buffered_extent)
  ))
  message("  Years       : ", min(years), " to ", max(years), " (", length(years), " yrs)")
  message("  Months      : ", paste(sprintf("%02d", months), collapse = ", "))
  message("  Batch Mode  : ", batch_mode, " (AgERA5 v", version, ")")
  message("  Output Dir  : ", standard_dir)
  message("========================================================")

  # 2. Check CDS Server Limitations
  check_cds_request_limits(
    years             = years,
    months            = months,
    statistics        = statistics,
    aoi_extent        = aoi$buffered_extent,
    batch_mode        = batch_mode,
    warn_max_years    = WARN_MAX_YEARS,
    warn_max_requests = WARN_MAX_REQUESTS,
    warn_large_area   = WARN_LARGE_AREA_SQDEG,
    confirm           = confirm_on_warning
  )

  # 3. Prompt for User's CDS Personal Access Token
  user_token <- prompt_cds_token(cds_token)

  # 4. Invoke Generic Python Downloader
  py_args <- c(
    DOWNLOADER_PY,
    if (nzchar(user_token)) c("--token", user_token) else character(0),
    "--cds-url", CDS_API_URL,
    "--extent",
    as.character(xmin(aoi$buffered_extent)),
    as.character(xmax(aoi$buffered_extent)),
    as.character(ymin(aoi$buffered_extent)),
    as.character(ymax(aoi$buffered_extent)),
    "--buffer", "0",
    "--years", paste(years, collapse = ","),
    "--months", paste(months, collapse = ","),
    "--batch-mode", batch_mode,
    "--variable", variable,
    "--statistics", statistics,
    "--version", version,
    "--output-dir", raw_zip_dir,
    "--tag", aoi$tag,
    "--yes"
  )

  status <- system2(PYTHON_BIN, args = py_args)
  if (!is.null(status) && status != 0) {
    warning("Python downloader exited with status code: ", status)
  }

  # 5. Standardize Downloaded Zips to Daily GeoTIFFs
  year_regex <- paste0("_(", paste(years, collapse = "|"), ")_.*", aoi$tag, "\\.zip$")
  zip_files <- sort(dir(raw_zip_dir, pattern = year_regex, full.names = TRUE))
  if (length(zip_files) > 0) {
    message("Standardizing ", length(zip_files), " zip archive(s) to daily GeoTIFFs...")
    standardize_agera5_zips(zip_files, out_dir = standard_dir, aoi_extent = aoi$buffered_extent, spatial_res = spatial_res)
    message("Standardization complete. Daily GeoTIFFs saved in: ", standard_dir)
  }

  invisible(standard_dir)
}

if (sys.nframe() == 0L) {
  run_agera5_pipeline()
}
