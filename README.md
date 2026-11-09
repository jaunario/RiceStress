# RiceStress

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24071671.svg)](https://doi.org/10.5281/zenodo.24071671)
[![License: CC0-1.0](https://img.shields.io/badge/License-CC0_1.0-lightgrey.svg)](https://creativecommons.org/publicdomain/zero/1.0/)

Spatial modeling scripts for assessing seasonal rice heat stress incidence (daytime and nighttime heat stress across reproductive and flowering stages) using Copernicus AgERA5 daily temperature rasters and rice crop calendar polygons.

---

## Repository Contents

| File | Description |
| :--- | :--- |
| `Download_Standardize_AgERA5.R` | Generic R pipeline to download AgERA5 v2.0 daily temperature NetCDF archives from the Copernicus Climate Data Store (CDS) and standardize them into `0.1°` (`6 arcmin`) GeoTIFF rasters (`AGERA5.6arcmin.<tmax\|tmin>.YYYY.MM.DD.tif`, °C $\times 100$, `INT2S`). |
| `download_agera5.py` | Python downloader (`cdsapi`) invoked by `Download_Standardize_AgERA5.R` or run standalone via CLI. Prompts for the user's CDS Personal Access Token, resolves Areas of Interest (AOI) from bounding boxes, vector files, or filtered Rice Crop Calendar polygons, warns of CDS server limitations, and reuses completed jobs in the user's CDS job history. |
| `Rice_HeatStress_Daily.R` | Seasonal rice heat stress mapping engine computing annual threshold-exceedance day counts and multi-year stress incidence frequency from daily weather rasters. |

---

## Software Environment & Package Versions

The scripts were developed and verified under the following software versions:

* **R**: `4.6.1 (2026-06-24)`
  * `terra`: `1.9.34`
* **Python**: `3.14.8`
  * `cdsapi`: `0.7.7`

Install required packages prior to running:

```r
# In R
install.packages("terra")
```

```bash
# In terminal
pip install "cdsapi>=0.7.7"
```

---

## Input Data Sources

### 1. Rice Crop Calendar (RiceAtlas)
If you do not have a custom rice crop calendar shapefile, download **RiceAtlas** from Harvard Dataverse:

* **Dataset (Harvard Dataverse)**: https://doi.org/10.7910/DVN/JE6R2R
* **Reference**: Laborte, A. G., Gutierrez, M. A., Balanza, J. G., Saito, K., Zwart, S. J., Boschetti, M., Murty, M. V. R., Villano, L., Aunario, J. K., Reinke, R., Koo, J., Hijmans, R. J., & Nelson, A. (2017). *RiceAtlas, a spatial database of global rice calendars and production*. **Scientific Data**, 4, 170074. https://doi.org/10.1038/sdata.2017.74

### 2. Daily Weather Data (Copernicus AgERA5 v2.0)
Daily minimum (`24_hour_minimum`) and maximum (`24_hour_maximum`) 2 m air temperature rasters are retrieved from the Copernicus Climate Data Store dataset `sis-agrometeorological-indicators` (version `2_0`).

* Register for a free Copernicus CDS account and obtain your **Personal Access Token** at: https://cds.climate.copernicus.eu/profile

---

## Usage Workflow

### Step 1: Download & Standardize Daily AgERA5 Rasters

Open `Download_Standardize_AgERA5.R` (or `source()` it in R) and configure your Area of Interest (AOI) using one of three supported methods:

1. **`SpatExtent` object** (or numeric `c(xmin, xmax, ymin, ymax)`):
   ```r
   AOI_EXTENT <- ext(120.8, 121.8, 13.8, 14.8)
   ```
2. **`SpatVector` object** (or path to `.shp` / `.gpkg` / `.geojson`):
   ```r
   AOI_VECTOR <- vect("path/to/study_area.shp")
   ```
3. **Rice Crop Calendar** (loaded as `SpatVector` or shapefile path, optionally filtered by attribute and active season):
   ```r
   source("Download_Standardize_AgERA5.R")

   run_agera5_pipeline(
     cropcal        = "path/to/RiceCalendar.shp",
     cropcal_filter = list(COUNTRY = "Philippines", SUB_REGION = "Laguna"),
     active_season  = 1,
     years          = 2020:2024,
     months         = 1:12,
     batch_mode     = "monthly",
     raw_zip_dir    = "./agera5_raw",
     standard_dir   = "./agera5_daily"
   )
   ```

When `CDS_API_TOKEN` is left as `NULL`, the script prompts interactively for your Copernicus CDS Personal Access Token. If the requested time period or number of API requests exceeds recommended CDS queue thresholds, the script warns of server limitations before queuing requests.

#### Standalone Python CLI Alternative
`download_agera5.py` can also be executed directly from the command line:

```bash
python3 download_agera5.py \
  --cropcal "path/to/RiceCalendar.shp" \
  --filter COUNTRY=Philippines --filter SUB_REGION=Laguna \
  --years 2020:2024 \
  --months 1:12 \
  --batch-mode monthly \
  --output-dir ./agera5_raw
```

---

### Step 2: Run Seasonal Rice Heat Stress Assessment

Configure the user settings at the top of `Rice_HeatStress_Daily.R` to point `SRC_DIRECTORY` to the standardized daily GeoTIFF folder (`./agera5_daily`) produced in Step 1:

```r
HEAT_STRESS_TYPE <- "daytime_heat"  # "daytime_heat" (DTH) or "nighttime_heat" (NTH)
CROP_PHASE       <- "reproductive"  # "reproductive", "flowering", or "custom"
SRC_DIRECTORY    <- "./agera5_daily"
ANALYSIS_YEARS   <- 2020:2024
ANALYSIS_SEASONS <- 1:3
CROPCAL_FILE     <- "path/to/RiceCalendar.shp"
RICEMASK_FILE    <- NULL            # Optional path to rice mask GeoTIFF
AOI_FILTER       <- list(field = "SUB_REGION", value = "Laguna") # NULL for full extent
```

#### Heat Stress Presets & Stage Windows

| Parameter | Daytime Heat (`DTH`) | Nighttime Heat (`NTH`) |
| :--- | :--- | :--- |
| **Daily Variable** | `tmax` (`24_hour_maximum`) | `tmin` (`24_hour_minimum`) |
| **Temperature Threshold** | $> 35\text{ }^\circ\text{C}$ | $> 25\text{ }^\circ\text{C}$ |
| **Season Stress Criterion** | $\ge 10\text{ days}$ | $\ge 15\text{ days}$ |

* **Reproductive Phase (`reproductive`)**: $\text{Harvest Date} - 56$ to $\text{Harvest Date} - 26$ days.
* **Flowering Phase (`flowering`)**: $\text{Harvest Date} - 34$ to $\text{Harvest Date} - 26$ days.

#### Output Products
Written to `OUTPUT_DIR` (`./out/<DTH|NTH>_AgERA5_Daily/`):
* `<DTH|NTH>_SEASON<1..3>_XTHRES_COUNT_Y<YYYY>.tif`: Annual count of days exceeding the temperature threshold during the target crop stage (`INT2S`).
* `<DTH|NTH>_SEASON<1..3>_AgERA5_STRESSYRS_COUNT_<START>-<END>.tif`: Total number of years in the analysis period where threshold exceedance days met or exceeded the critical day-count threshold (`INT2S`).
* `Masked/<DTH|NTH>_SEASON<1..3>_AgERA5_STRESSYRS_COUNT_<START>-<END>_masked.tif`: Rice-mask-filtered stressed years count (generated when `RICEMASK_FILE` is supplied).

---

## Citation

If you use this software in your research, please cite it as:

```bibtex
@software{aunario_ricestress_2026,
  author       = {Aunario, Jorrel Khalil S.},
  title        = {RiceStress: Seasonal Rice Heat Stress Assessment Using Daily Temperature Data},
  month        = nov,
  year         = 2026,
  publisher    = {Zenodo},
  version      = {v1.0.0},
  doi          = {10.5281/zenodo.24071671},
  url          = {https://doi.org/10.5281/zenodo.24071671}
}
```

---

## License

This work is dedicated to the public domain under the [CC0 1.0 Universal (CC0 1.0) Public Domain Dedication](LICENSE).

