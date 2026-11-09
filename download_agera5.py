#!/usr/bin/env python3
"""Generic AgERA5 Daily Downloader (Copernicus CDS API).

Features:
1. Prompts the user for their Copernicus CDS Personal Access Token (rather than
   relying on pre-existing developer credentials).
2. Supports Area of Interest (AOI) selection via:
   - Extent object / bounding box (xmin, xmax, ymin, ymax)
   - Vector dataset / SpatVector file (.shp, .gpkg, .geojson)
   - Rice Crop Calendar shapefile with attribute and active-season filtering
     (e.g., RiceAtlas: Laborte et al., 2017, https://doi.org/10.1038/sdata.2017.74;
     downloadable from Harvard Dataverse at https://doi.org/10.7910/DVN/JE6R2R)
3. Supports flexible time periods (year ranges, specific months, or date ranges)
   with 'yearly' or 'monthly' request batching.
4. Warns of Copernicus CDS server limitations (queue throttling, per-request
   field/volume caps, and multi-year rate limits) when the requested period or
   spatial-temporal volume is large.
5. Reuses local zip cache and existing completed jobs on the user's CDS account
   (via get_jobs / get_receipt / get_results) before submitting new requests.

Environment & Package Versions:
- Python : 3.14.8
- cdsapi : 0.7.7
- R      : 4.6.1 (2026-06-24) with terra 1.9.34 (optional fallback for non-.shp vector formats)
"""

import argparse
from datetime import date, datetime
import getpass
import os
from pathlib import Path
import struct
import subprocess
import sys
from typing import Dict, List, Optional, Sequence, Tuple, Union

import cdsapi

DEFAULT_CDS_URL = "https://cds.climate.copernicus.eu/api"

# Thresholds for Copernicus CDS server limitation warnings
WARN_MAX_YEARS = 2
WARN_MAX_REQUESTS = 12
WARN_MAX_TOTAL_DAYS = 730  # ~2 years of daily data
WARN_LARGE_AREA_SQ_DEG = 2500.0  # e.g., > 50 deg x 50 deg


class SpatExtent:
    """Lightweight Extent object matching terra::SpatExtent (xmin, xmax, ymin, ymax)."""

    def __init__(self, xmin: float, xmax: float, ymin: float, ymax: float):
        if xmin > xmax:
            xmin, xmax = xmax, xmin
        if ymin > ymax:
            ymin, ymax = ymax, ymin
        self.xmin = float(xmin)
        self.xmax = float(xmax)
        self.ymin = float(ymin)
        self.ymax = float(ymax)

    def buffer(self, deg: float) -> "SpatExtent":
        """Pad extent by `deg` degrees and clamp to valid global lon/lat bounds."""
        return SpatExtent(
            max(-180.0, round(self.xmin - deg, 4)),
            min(180.0, round(self.xmax + deg, 4)),
            max(-90.0, round(self.ymin - deg, 4)),
            min(90.0, round(self.ymax + deg, 4)),
        )

    def to_cds_area(self) -> List[float]:
        """Convert to Copernicus CDS bounding box format: [North, West, South, East]."""
        return [
            round(self.ymax, 2),
            round(self.xmin, 2),
            round(self.ymin, 2),
            round(self.xmax, 2),
        ]

    @property
    def area_sq_deg(self) -> float:
        return max(0.0, (self.xmax - self.xmin) * (self.ymax - self.ymin))

    @property
    def approx_cells_01deg(self) -> int:
        nx = max(1, int(round((self.xmax - self.xmin) / 0.1)))
        ny = max(1, int(round((self.ymax - self.ymin) / 0.1)))
        return nx * ny

    def __repr__(self) -> str:
        return (
            f"SpatExtent(xmin={self.xmin:.4f}, xmax={self.xmax:.4f}, "
            f"ymin={self.ymin:.4f}, ymax={self.ymax:.4f})"
        )


# ==============================================================================
# 1. SHAPEFILE / SPATVECTOR & CROP CALENDAR AOI RESOLUTION
# ==============================================================================


def _read_dbf_records(dbf_path: Path) -> Tuple[List[str], List[Dict[str, object]]]:
    """Read field names and attribute records from an ESRI Shapefile .dbf file."""
    with open(dbf_path, "rb") as f:
        header = f.read(32)
        num_records, header_len, record_len = struct.unpack("<4xIHH20x", header)
        num_fields = (header_len - 33) // 32

        fields = []
        for _ in range(num_fields):
            field_desc = f.read(32)
            name = field_desc[:11].split(b"\x00", 1)[0].decode("latin1").strip()
            ftype = chr(field_desc[11])
            flen = field_desc[16]
            fdec = field_desc[17]
            fields.append((name, ftype, flen, fdec))

        f.seek(header_len)
        records = []
        for _ in range(num_records):
            raw = f.read(record_len)
            if not raw or raw[0:1] == b"*":
                continue
            pos = 1
            rec: Dict[str, object] = {}
            for name, ftype, flen, fdec in fields:
                val_bytes = raw[pos : pos + flen]
                pos += flen
                txt = val_bytes.decode("latin1", errors="replace").strip()
                if ftype in ("N", "F"):
                    if txt == "":
                        rec[name] = None
                    elif fdec > 0 or "." in txt:
                        try:
                            rec[name] = float(txt)
                        except ValueError:
                            rec[name] = txt
                    else:
                        try:
                            rec[name] = int(txt)
                        except ValueError:
                            rec[name] = txt
                else:
                    rec[name] = txt
            records.append(rec)

    return [f[0] for f in fields], records


def _read_shp_bboxes(shp_path: Path) -> List[Optional[SpatExtent]]:
    """Read per-feature bounding boxes from an ESRI .shp/.shx pair."""
    shx_path = shp_path.with_suffix(".shx")
    bboxes: List[Optional[SpatExtent]] = []

    with open(shx_path, "rb") as f_shx, open(shp_path, "rb") as f_shp:
        shx_header = f_shx.read(100)
        if len(shx_header) < 100:
            raise ValueError(f"Invalid .shx file: {shx_path}")

        while True:
            idx_rec = f_shx.read(8)
            if len(idx_rec) < 8:
                break
            offset_words, content_len_words = struct.unpack(">II", idx_rec)
            byte_offset = offset_words * 2
            f_shp.seek(byte_offset + 8)
            shape_type_bytes = f_shp.read(4)
            if len(shape_type_bytes) < 4:
                bboxes.append(None)
                continue
            (shape_type,) = struct.unpack("<i", shape_type_bytes)
            if shape_type == 0:
                bboxes.append(None)
            elif shape_type in (1, 11, 21):
                x, y = struct.unpack("<2d", f_shp.read(16))
                bboxes.append(SpatExtent(x, x, y, y))
            elif shape_type in (3, 5, 8, 13, 15, 18, 23, 25, 28, 31):
                xmin, ymin, xmax, ymax = struct.unpack("<4d", f_shp.read(32))
                bboxes.append(SpatExtent(xmin, xmax, ymin, ymax))
            else:
                bboxes.append(None)

    return bboxes


def resolve_vector_or_cropcal_extent(
    vector_path: Union[str, Path],
    filters: Optional[Dict[str, Union[str, Sequence[str]]]] = None,
    active_season: Optional[int] = None,
) -> Tuple[SpatExtent, int]:
    """Compute bounding box from a vector file or filtered Rice Crop Calendar."""
    path = Path(vector_path)
    if not path.exists():
        raise FileNotFoundError(
            f"Vector/Crop Calendar file not found: {path}\n"
            "If using the global Rice Crop Calendar, download RiceAtlas (Laborte et al., 2017) "
            "from Harvard Dataverse: https://doi.org/10.7910/DVN/JE6R2R"
        )

    # Direct pure-Python reader for .shp files
    if path.suffix.lower() == ".shp":
        dbf_path = path.with_suffix(".dbf")
        bboxes = _read_shp_bboxes(path)
        fields, records = _read_dbf_records(dbf_path)

        if filters:
            for k in filters:
                if k not in fields:
                    raise ValueError(
                        f"Filter field '{k}' not found in {path.name}. "
                        f"Available fields: {', '.join(fields)}"
                    )

        matched_extents: List[SpatExtent] = []
        for rec, bbox in zip(records, bboxes):
            if bbox is None:
                continue
            keep = True
            if filters:
                for field, target in filters.items():
                    rec_val = str(rec.get(field, "")).strip().lower()
                    if isinstance(target, (list, tuple, set)):
                        allowed = {str(t).strip().lower() for t in target}
                        if rec_val not in allowed:
                            keep = False
                            break
                    else:
                        if rec_val != str(target).strip().lower():
                            keep = False
                            break
            if keep and active_season is not None:
                pl_col = f"PLANT_PK{active_season}"
                hv_col = f"HARV_PK{active_season}"
                pl_val = rec.get(pl_col, 0) or 0
                hv_val = rec.get(hv_col, 0) or 0
                if not (float(pl_val) > 0 and float(hv_val) > 0):
                    keep = False

            if keep:
                matched_extents.append(bbox)

        if not matched_extents:
            raise ValueError(
                f"No features matched in {path.name} for filters={filters}, "
                f"active_season={active_season}."
            )

        xmin = min(e.xmin for e in matched_extents)
        xmax = max(e.xmax for e in matched_extents)
        ymin = min(e.ymin for e in matched_extents)
        ymax = max(e.ymax for e in matched_extents)
        return SpatExtent(xmin, xmax, ymin, ymax), len(matched_extents)

    # Fallback via R terra::vect() for non-shapefile formats (.gpkg, .geojson, etc.)
    r_expr = (
        f"suppressPackageStartupMessages(library(terra)); "
        f"v <- vect('{path.as_posix()}'); "
        f"e <- ext(v); "
        f"cat(sprintf('%f,%f,%f,%f,%d\\n', e$xmin, e$xmax, e$ymin, e$ymax, nrow(v)))"
    )
    res = subprocess.run(
        ["Rscript", "-e", r_expr], capture_output=True, text=True, check=True
    )
    parts = res.stdout.strip().splitlines()[-1].split(",")
    return (
        SpatExtent(float(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])),
        int(parts[4]),
    )


def resolve_aoi(
    extent: Optional[Union[SpatExtent, Sequence[float]]] = None,
    vector_path: Optional[Union[str, Path]] = None,
    cropcal_path: Optional[Union[str, Path]] = None,
    cropcal_filters: Optional[Dict[str, Union[str, Sequence[str]]]] = None,
    active_season: Optional[int] = None,
    buffer_deg: float = 0.2,
) -> Tuple[SpatExtent, str]:
    """Resolve user AOI from an extent object, a SpatVector/vector file, or Crop Calendar."""
    if extent is not None:
        if isinstance(extent, SpatExtent):
            raw_ext = extent
        elif len(extent) == 4:
            raw_ext = SpatExtent(extent[0], extent[1], extent[2], extent[3])
        else:
            raise ValueError("extent must be a SpatExtent or [xmin, xmax, ymin, ymax].")
        buffered = raw_ext.buffer(buffer_deg)
        return buffered, f"Extent {raw_ext} (buffered by {buffer_deg} deg)"

    if vector_path is not None:
        raw_ext, n_feat = resolve_vector_or_cropcal_extent(
            vector_path, filters=cropcal_filters, active_season=active_season
        )
        buffered = raw_ext.buffer(buffer_deg)
        return (
            buffered,
            f"SpatVector '{Path(vector_path).name}' ({n_feat} features, buffered by {buffer_deg} deg)",
        )

    if cropcal_path is not None:
        raw_ext, n_feat = resolve_vector_or_cropcal_extent(
            cropcal_path, filters=cropcal_filters, active_season=active_season
        )
        buffered = raw_ext.buffer(buffer_deg)
        desc = (
            f"Rice Crop Calendar '{Path(cropcal_path).name}' "
            f"(filters={cropcal_filters or 'ALL'}, season={active_season or 'ALL'}, "
            f"{n_feat} units, buffered by {buffer_deg} deg)"
        )
        return buffered, desc

    raise ValueError(
        "No Area of Interest specified. Provide --extent, --vector, or --cropcal."
    )


# ==============================================================================
# 2. CDS PERSONAL ACCESS TOKEN AUTHENTICATION
# ==============================================================================


def get_cds_client(
    token: Optional[str] = None,
    cds_url: str = DEFAULT_CDS_URL,
    save_cdsapirc: Optional[Path] = None,
) -> cdsapi.Client:
    """Authenticate with Copernicus CDS using the user's Personal Access Token."""
    api_token = (token or os.environ.get("CDS_API_TOKEN") or "").strip()
    default_rc = Path.home() / ".cdsapirc"

    if not api_token:
        if sys.stdin.isatty():
            print("\n=== Copernicus CDS Authentication ===")
            print(f"CDS API Endpoint: {cds_url}")
            print(
                "Find your Personal Access Token at: https://cds.climate.copernicus.eu/profile"
            )
            prompt_msg = (
                "Enter your CDS Personal Access Token (or press Enter to use ~/.cdsapirc): "
                if default_rc.exists()
                else "Enter your CDS Personal Access Token: "
            )
            api_token = getpass.getpass(prompt_msg).strip()
        elif not default_rc.exists():
            raise RuntimeError(
                "CDS Personal Access Token is required in non-interactive mode. "
                "Pass --token <TOKEN> or set the CDS_API_TOKEN environment variable."
            )

    if not api_token:
        if default_rc.exists():
            return cdsapi.Client()
        raise ValueError(
            "No CDS Personal Access Token provided. Aborting download."
        )

    if save_cdsapirc is not None:
        rc_path = Path(save_cdsapirc).expanduser()
        rc_path.write_text(f"url: {cds_url}\nkey: {api_token}\n", encoding="utf-8")
        os.chmod(rc_path, 0o600)
        print(f"Saved CDS credentials to: {rc_path}")

    return cdsapi.Client(url=cds_url, key=api_token)


# ==============================================================================
# 3. TIME PERIOD PARSING & SERVER LIMITATION WARNINGS
# ==============================================================================


def parse_int_sequence(spec: str, min_val: int, max_val: int, name: str) -> List[int]:
    """Parse integer range/list expressions like '1993:2023', '1993-2023', or '2020,2021,2022'."""
    values: List[int] = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if ":" in part or "-" in part:
            sep = ":" if ":" in part else "-"
            start_s, end_s = part.split(sep, 1)
            start_i, end_i = int(start_s), int(end_s)
            if start_i > end_i:
                start_i, end_i = end_i, start_i
            values.extend(range(start_i, end_i + 1))
        else:
            values.append(int(part))

    unique_vals = sorted(set(values))
    for v in unique_vals:
        if v < min_val or v > max_val:
            raise ValueError(f"{name} value {v} is outside valid range [{min_val}, {max_val}].")
    return unique_vals


def check_cds_server_limitations(
    years: Sequence[int],
    months: Sequence[int],
    statistics: Sequence[str],
    aoi_extent: SpatExtent,
    batch_mode: str,
    auto_confirm: bool = False,
) -> None:
    """Warn the user when the requested time period or volume risks CDS server throttling."""
    n_years = len(years)
    n_months_per_year = len(months)
    total_months = n_years * n_months_per_year
    approx_total_days = int(round(total_months * 30.4375))
    n_requests = n_years if batch_mode == "yearly" else total_months
    fields_per_request = (
        int(round(n_months_per_year * 30.4375)) * len(statistics)
        if batch_mode == "yearly"
        else 31 * len(statistics)
    )

    warnings_triggered: List[str] = []

    if n_years > WARN_MAX_YEARS or approx_total_days > WARN_MAX_TOTAL_DAYS:
        warnings_triggered.append(
            f"- Long time period requested ({n_years} years, ~{approx_total_days} days across "
            f"{years[0]}-{years[-1]}). Copernicus CDS enforces strict per-user queue concurrency "
            f"and rate limits; multi-year batches can experience long queue wait times (hours) "
            f"or temporary rate throttling."
        )

    if n_requests > WARN_MAX_REQUESTS:
        warnings_triggered.append(
            f"- High API request count ({n_requests} separate CDS retrieve jobs). "
            f"Submitting >{WARN_MAX_REQUESTS} jobs sequentially may trigger CDS queue delays. "
            f"Consider downloading in smaller year chunks (e.g., 5-10 years at a time) or using "
            f"'yearly' batching for regional AOIs."
        )

    if aoi_extent.area_sq_deg > WARN_LARGE_AREA_SQ_DEG and batch_mode == "yearly":
        warnings_triggered.append(
            f"- Large spatial extent ({aoi_extent.area_sq_deg:.1f} sq. deg, "
            f"~{aoi_extent.approx_cells_01deg:,} grid cells at 0.1 deg) combined with 'yearly' "
            f"batching ({fields_per_request} fields/request). CDS may reject large payloads or "
            f"time out; '--batch-mode monthly' is strongly recommended for continental/global extents."
        )

    if warnings_triggered:
        print(
            "\n======================================================================",
            file=sys.stderr,
        )
        print(
            "WARNING: COPERNICUS CDS SERVER & API REQUEST LIMITATIONS",
            file=sys.stderr,
        )
        print(
            "======================================================================",
            file=sys.stderr,
        )
        for w in warnings_triggered:
            print(w, file=sys.stderr)
        print(
            "----------------------------------------------------------------------",
            file=sys.stderr,
        )
        print(
            f"Summary: {n_requests} CDS API requests | Years: {years[0]}-{years[-1]} ({n_years} yrs) | "
            f"Months/yr: {n_months_per_year} | Batch mode: {batch_mode}",
            file=sys.stderr,
        )
        print(
            "======================================================================\n",
            file=sys.stderr,
        )

        if sys.stdin.isatty() and not auto_confirm:
            ans = (
                input(f"Proceed with submitting up to {n_requests} requests to CDS? [y/N]: ")
                .strip()
                .lower()
            )
            if ans not in ("y", "yes"):
                print("Download cancelled by user.")
                sys.exit(0)


# ==============================================================================
# 4. CDS SERVER-SIDE JOB CACHE REUSE & DOWNLOAD EXECUTION
# ==============================================================================


def normalize_param(val: object) -> Union[List[str], str]:
    """Normalize request parameters for exact comparison against CDS receipts."""
    if isinstance(val, (list, tuple)):
        return sorted([str(v) for v in val])
    return str(val)


def find_existing_server_job(
    client: cdsapi.Client, target_req: Dict[str, object], cached_jobs: List[Dict[str, object]]
) -> Optional[str]:
    """Scan existing successful jobs on the user's CDS account for matching parameters."""
    target_year = normalize_param(target_req.get("year"))
    target_month = normalize_param(target_req.get("month"))
    target_day = normalize_param(target_req.get("day"))
    target_var = normalize_param(target_req.get("variable"))
    target_stat = normalize_param(target_req.get("statistic"))
    target_area = [round(float(x), 2) for x in target_req.get("area", [])]  # type: ignore[arg-type]

    for job in cached_jobs:
        if job.get("status") != "successful":
            continue
        job_id = job.get("jobID")
        if not job_id:
            continue
        try:
            receipt = client.client.get_receipt(job_id)
            req = receipt.get("request", {})
            if (
                normalize_param(req.get("year")) == target_year
                and normalize_param(req.get("month")) == target_month
                and normalize_param(req.get("day")) == target_day
                and normalize_param(req.get("variable")) == target_var
                and normalize_param(req.get("statistic")) == target_stat
            ):
                req_area = [round(float(x), 2) for x in req.get("area", [])]
                if req_area == target_area:
                    return str(job_id)
        except Exception:
            continue
    return None


def download_agera5(
    client: cdsapi.Client,
    aoi_extent: SpatExtent,
    years: Sequence[int],
    months: Sequence[int],
    output_dir: Union[str, Path],
    region_tag: str = "aoi",
    variable: str = "2m_temperature",
    statistics: Sequence[str] = ("24_hour_minimum", "24_hour_maximum"),
    batch_mode: str = "monthly",
    version: str = "2_0",
    check_server_cache: bool = True,
) -> List[Path]:
    """Execute AgERA5 downloads for the resolved AOI and time period."""
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    cds_area = aoi_extent.to_cds_area()
    days_list = [f"{d:02d}" for d in range(1, 32)]
    months_str = [f"{m:02d}" for m in months]
    stats_list = list(statistics)

    server_jobs: List[Dict[str, object]] = []
    if check_server_cache:
        print("Checking recent completed jobs on your Copernicus CDS account...")
        try:
            jobs_resp = client.client.get_jobs()
            server_jobs = jobs_resp.json.get("jobs", [])
            print(f"Found {len(server_jobs)} past jobs in CDS account history.")
        except Exception as exc:
            print(f"Note: Could not query CDS job history ({exc}). Proceeding with standard requests.")

    downloaded_files: List[Path] = []

    # Build batch schedule (yearly vs monthly)
    batches: List[Tuple[str, List[str], Path]] = []
    clean_tag = "".join(c if c.isalnum() or c in ("-", "_") else "_" for c in region_tag).strip("_")
    for yr in years:
        if batch_mode == "yearly":
            zip_file = out_path / f"agera5_{variable}_{yr}_{clean_tag}.zip"
            batches.append((str(yr), months_str, zip_file))
        else:
            for m_str in months_str:
                zip_file = out_path / f"agera5_{variable}_{yr}_{m_str}_{clean_tag}.zip"
                batches.append((str(yr), [m_str], zip_file))

    for idx, (yr_str, req_months, zip_path) in enumerate(batches, start=1):
        label = f"{yr_str}" if batch_mode == "yearly" else f"{yr_str}-{req_months[0]}"

        # 1. Local disk cache check
        if zip_path.exists() and zip_path.stat().st_size > 1024:
            print(f"[{idx}/{len(batches)}] [{label}] Local file exists: {zip_path.name}. Skipping.")
            downloaded_files.append(zip_path)
            continue

        req: Dict[str, object] = {
            "variable": variable,
            "statistic": stats_list if len(stats_list) > 1 else stats_list[0],
            "year": yr_str,
            "month": req_months if len(req_months) > 1 else req_months[0],
            "day": days_list,
            "version": version,
            "area": cds_area,
            "format": "zip",
        }

        # 2. Server-side completed job cache check
        if server_jobs:
            existing_job_id = find_existing_server_job(client, req, server_jobs)
            if existing_job_id:
                print(
                    f"[{idx}/{len(batches)}] [{label}] Found matching completed CDS job "
                    f"({existing_job_id}). Downloading directly..."
                )
                try:
                    results = client.client.get_results(existing_job_id)
                    results.download(str(zip_path))
                    print(f"[{idx}/{len(batches)}] [{label}] Downloaded: {zip_path.name}")
                    downloaded_files.append(zip_path)
                    continue
                except Exception as exc:
                    print(
                        f"[{idx}/{len(batches)}] [{label}] Cached job download failed ({exc}), "
                        f"submitting new request..."
                    )

        # 3. Submit new CDS request
        print(f"[{idx}/{len(batches)}] [{label}] Submitting AgERA5 request to CDS (area={cds_area})...")
        try:
            client.retrieve("sis-agrometeorological-indicators", req, str(zip_path))
            print(f"[{idx}/{len(batches)}] [{label}] Successfully downloaded: {zip_path.name}")
            downloaded_files.append(zip_path)
        except Exception as exc:
            print(f"[{idx}/{len(batches)}] [{label}] Error downloading: {exc}", file=sys.stderr)

    return downloaded_files


# ==============================================================================
# 5. CLI ENTRY POINT
# ==============================================================================


def parse_filter_args(filter_items: Optional[Sequence[str]]) -> Optional[Dict[str, Union[str, List[str]]]]:
    """Parse CLI '--filter FIELD=VALUE' arguments into a dictionary."""
    if not filter_items:
        return None
    out: Dict[str, Union[str, List[str]]] = {}
    for item in filter_items:
        if "=" not in item:
            raise ValueError(f"Invalid --filter '{item}'. Expected format: FIELD=VALUE")
        key, val = item.split("=", 1)
        key = key.strip()
        vals = [v.strip() for v in val.split(",") if v.strip()]
        out[key] = vals if len(vals) > 1 else vals[0]
    return out


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generic AgERA5 Daily Downloader (Extent / SpatVector / Rice Crop Calendar)."
    )

    # Authentication
    auth_grp = parser.add_argument_group("Authentication (Copernicus CDS)")
    auth_grp.add_argument(
        "--token",
        type=str,
        default=None,
        help="Copernicus CDS Personal Access Token. If omitted, prompts interactively.",
    )
    auth_grp.add_argument(
        "--cds-url",
        type=str,
        default=DEFAULT_CDS_URL,
        help=f"CDS API URL (default: {DEFAULT_CDS_URL}).",
    )
    auth_grp.add_argument(
        "--save-cdsapirc",
        type=str,
        default=None,
        help="Optional path to save the entered CDS token in .cdsapirc format.",
    )

    # Area of Interest (AOI)
    aoi_grp = parser.add_argument_group("Area of Interest (Extent / SpatVector / Crop Calendar)")
    aoi_grp.add_argument(
        "--extent",
        nargs=4,
        type=float,
        metavar=("XMIN", "XMAX", "YMIN", "YMAX"),
        help="Bounding box extent in decimal degrees: xmin xmax ymin ymax.",
    )
    aoi_grp.add_argument(
        "--vector",
        type=str,
        default=None,
        help="Path to a vector file (.shp, .gpkg, .geojson) to extract AOI extent.",
    )
    aoi_grp.add_argument(
        "--cropcal",
        type=str,
        default=None,
        help=(
            "Path to a Rice Crop Calendar shapefile (.shp) to extract AOI extent. "
            "RiceAtlas (Laborte et al., 2017; https://doi.org/10.1038/sdata.2017.74) "
            "can be downloaded from Harvard Dataverse: https://doi.org/10.7910/DVN/JE6R2R"
        ),
    )
    aoi_grp.add_argument(
        "--filter",
        action="append",
        metavar="FIELD=VALUE",
        help="Attribute filter(s) on the Crop Calendar or vector file (e.g., --filter SUB_REGION=Laguna or --filter COUNTRY=Philippines).",
    )
    aoi_grp.add_argument(
        "--season",
        type=int,
        choices=[1, 2, 3],
        default=None,
        help="Optional Rice Crop Calendar season (1, 2, or 3) requiring active planting/harvest dates (PLANT_PKn > 0).",
    )
    aoi_grp.add_argument(
        "--buffer",
        type=float,
        default=0.2,
        help="Spatial buffer in degrees added around the resolved extent (default: 0.2).",
    )

    # Time Period
    time_grp = parser.add_argument_group("Time Period")
    time_grp.add_argument(
        "--years",
        type=str,
        default=None,
        help="Year range or list, e.g., '2020:2023', '1993-2023', or '2021,2022,2023'.",
    )
    time_grp.add_argument("--start-year", type=int, default=None, help="Start year (e.g., 1993).")
    time_grp.add_argument("--end-year", type=int, default=None, help="End year (e.g., 2023).")
    time_grp.add_argument(
        "--months",
        type=str,
        default="1:12",
        help="Month range or list, e.g., '1:12' or '5,6,7,8,9,10' (default: '1:12').",
    )
    time_grp.add_argument(
        "--batch-mode",
        choices=["monthly", "yearly"],
        default="monthly",
        help="Submit 1 CDS request per month ('monthly', default to stay within CDS cost limits) or per year ('yearly').",
    )

    # Dataset & Output
    out_grp = parser.add_argument_group("Dataset & Output")
    out_grp.add_argument(
        "--variable",
        type=str,
        default="2m_temperature",
        help="AgERA5 variable name (default: 2m_temperature).",
    )
    out_grp.add_argument(
        "--statistics",
        nargs="+",
        default=["24_hour_minimum", "24_hour_maximum"],
        help="AgERA5 statistic(s) (default: 24_hour_minimum 24_hour_maximum).",
    )
    out_grp.add_argument(
        "--version",
        type=str,
        default="2_0",
        help="AgERA5 dataset version on CDS (default: 2_0).",
    )
    out_grp.add_argument(
        "--output-dir",
        type=str,
        default="./agera5_raw",
        help="Directory to store downloaded AgERA5 .zip archives (default: ./agera5_raw).",
    )
    out_grp.add_argument(
        "--tag",
        type=str,
        default=None,
        help="Short identifier tag for output filenames (auto-derived from AOI if omitted).",
    )
    out_grp.add_argument(
        "--yes",
        action="store_true",
        help="Skip interactive confirmation prompt when CDS server limitation warnings trigger.",
    )
    out_grp.add_argument(
        "--no-server-cache",
        action="store_true",
        help="Disable checking CDS account job history for previously completed identical requests.",
    )

    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    # 1. Resolve Time Period
    if args.years:
        years = parse_int_sequence(args.years, 1979, date.today().year, "Year")
    elif args.start_year is not None and args.end_year is not None:
        years = parse_int_sequence(
            f"{args.start_year}:{args.end_year}", 1979, date.today().year, "Year"
        )
    elif args.start_year is not None:
        years = [args.start_year]
    else:
        parser.error("Specify time period via --years (e.g. 2020:2023) or --start-year / --end-year.")

    months = parse_int_sequence(args.months, 1, 12, "Month")

    # 2. Resolve AOI
    filters = parse_filter_args(args.filter)
    aoi_extent, aoi_desc = resolve_aoi(
        extent=args.extent,
        vector_path=args.vector,
        cropcal_path=args.cropcal,
        cropcal_filters=filters,
        active_season=args.season,
        buffer_deg=args.buffer,
    )

    if args.tag:
        region_tag = args.tag
    elif filters:
        first_val = next(iter(filters.values()))
        region_tag = first_val[0] if isinstance(first_val, list) else str(first_val)
    elif args.vector:
        region_tag = Path(args.vector).stem
    else:
        region_tag = (
            f"bbox_{aoi_extent.xmin:.1f}_{aoi_extent.xmax:.1f}_"
            f"{aoi_extent.ymin:.1f}_{aoi_extent.ymax:.1f}"
        )

    print("======================================================================")
    print("  Generic AgERA5 Daily Downloader")
    print(f"  AOI Source   : {aoi_desc}")
    print(f"  SpatExtent   : {aoi_extent}")
    print(f"  CDS Box      : {aoi_extent.to_cds_area()} [North, West, South, East]")
    print(f"  Years        : {years[0]} to {years[-1]} ({len(years)} year(s))")
    print(f"  Months       : {', '.join(f'{m:02d}' for m in months)}")
    print(f"  Variable     : {args.variable} ({', '.join(args.statistics)})")
    print(f"  Output Dir   : {args.output_dir}")
    print("======================================================================")

    # 3. Check CDS Server Limitations
    check_cds_server_limitations(
        years=years,
        months=months,
        statistics=args.statistics,
        aoi_extent=aoi_extent,
        batch_mode=args.batch_mode,
        auto_confirm=args.yes,
    )

    # 4. Prompt for User's CDS Personal Access Token & Initialize Client
    client = get_cds_client(
        token=args.token,
        cds_url=args.cds_url,
        save_cdsapirc=Path(args.save_cdsapirc) if args.save_cdsapirc else None,
    )

    # 5. Run Downloads
    download_agera5(
        client=client,
        aoi_extent=aoi_extent,
        years=years,
        months=months,
        output_dir=args.output_dir,
        region_tag=region_tag,
        variable=args.variable,
        statistics=args.statistics,
        batch_mode=args.batch_mode,
        version=args.version,
        check_server_cache=not args.no_server_cache,
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())

