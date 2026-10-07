# ==========================================================================
# PUBLIC RELEASE. Generated from the canonical analysis.R by
# scripts/make_public_code.py - do not edit this copy by hand; edit
# analysis.R and regenerate. Analysis logic is identical to the canonical
# script; only machine-specific paths and internal notes differ.
# Input data files are not distributed (see README): the release holds
# the pipeline and the manuscript's published derived tables only.
# Unzip derived_tables.zip to output/tables/ to compare a rerun.
# License: MIT (full text in README.md).
# Run from the repository root: Rscript analysis.R
# ==========================================================================
# ==========================================================================
# Prehospital ambulance access to trauma, STEMI, and stroke care in
# Saudi Arabia: a national two-leg spatial model of the SRCA network
#
# Single-file analysis pipeline. Sections:
#   0. Setup and parameters (incl. SRCA 2025 report observed moments)
#   1. Supply: 509 SRCA ambulance centres; destination sets for trauma
#      (MoH general/tertiary hospitals), STEMI (73 cath centres), stroke
#      (71 SRCA-receiving hospitals; optional reperfusion-capable Tier 2)
#   2. Governorates, regions, census population (demand)
#   3. Gridded population (WorldPop scaled to GASTAT 2022 totals)
#   4. Travel times (dodgr street-network routing, no external server)
#      Leg A: nearest SRCA centre -> population cell   (station-to-scene)
#      Leg B: population cell -> nearest facility      (scene-to-door)
#   5. Calibration (demand-weighted, one moment) and held-out validation
#   6. Full-chain call-to-door times and standards compliance
#   7. 8-minute response KPI: structural vs operational decomposition,
#      incl. a KPI-matched operational-delay scenario (lower bound)
#   8. Bypass analysis (mothership vs nearest-ED penalty)
#   9. Inequality (weighted Gini; Theil between/within regions + governorates)
#  9b. Spatial clustering test (global Moran's I + LISA, permutation)
#  10. Sensitivity analyses (on-scene time, speed, denominator, >=100-bed
#      trauma set, Leg-B speed factor, operational delay, urban-rural)
#  11. Tables and figures (5 main exhibits + supplement CSVs)
#  12. Decision analysis: reach vs delivery leg decomposition; greedy
#      capability-upgrade site rankings (STEMI cath, stroke receiving)
#      across the Health Holding Company (MoH) general-hospital network
#
# Model: call-to-door time for each 1-km populated cell =
#   fixed pre-departure intervals (queue 0:48 + mobilization 0:27, SRCA
#   2025 report, Table 17/2) + alpha * LegA + on-scene time (10/15/20 min
#   scenarios) + alpha * LegB, where alpha is a single speed factor
#   calibrated so the demand-weighted mean of alpha*LegA equals SRCA's
#   observed mean travel-to-scene (7:46). Remaining report moments (45.43%
#   of responses <= 8 min; stroke-pathway response 11 min) are held out
#   for validation, never fitted.
#
# Open decisions are parameters, not code changes:
#   PAR$trauma_set   -- "general" (default) | "all" | "beds100"
#   PAR$scene_min    -- 15 (default) | 10 | 20
#   Denominator      -- total residents primary; Saudi-only is a sensitivity
#   Stroke Tier 2    -- auto-runs iff a reperfusion list is supplied (see S1)
#
# Expensive steps cache to data/processed/ as .rds; delete a cache file to
# force recomputation. Section 4 needs the OSM road extract (see its header);
# everything else is plain CRAN R.
# ==========================================================================

# ---- 0. Setup and parameters ---------------------------------------------

pkgs <- c("sf", "terra", "dplyr", "tidyr", "readr", "stringr", "dodgr",
          "osmextract", "units", "RANN", "curl", "ggplot2", "patchwork",
          "viridis", "jsonlite", "spdep")   # after 13: `::` only
new <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs[1:13], library, character.only = TRUE))

root <- "."   # project root (srca/); auto-detected below when the script is
              # source()'d or Rscript'ed from elsewhere
if (!dir.exists(file.path(root, "raw"))) {
  f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)  # source() path
  if (is.null(f) || !nzchar(f)) {
    a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(a)) f <- sub("^--file=", "", a[1])               # Rscript path
  }
  cands <- character()
  if (!is.null(f) && nzchar(f)) {
    d <- dirname(normalizePath(f))
    cands <- c(d, dirname(d))
  }
  for (cand in cands)
    if (dir.exists(file.path(cand, "raw"))) { root <- cand; break }
  if (!dir.exists(file.path(root, "raw")))
    stop("raw/ not found. setwd() to the srca project root or run via ",
         "source(\"<full path to analysis.R>\").")
}
raw  <- file.path(root, "raw")
prc  <- file.path(root, "data/processed")
figs <- file.path(root, "output/figures")
tabs <- file.path(root, "output/tables")
for (d in c(prc, figs, tabs)) dir.create(d, recursive = TRUE,
                                         showWarnings = FALSE)

PAR <- list(
  crs_m       = 32638,          # UTM 38N (metric, covers KSA core)
  knn_candidates = 10,          # Euclidean candidates per origin
  pbf_path    = "raw/gcc-states-latest.osm.pbf",          # Geofabrik Gulf
  hw_keep     = c("motorway", "motorway_link", "trunk", "trunk_link",
                  "primary", "primary_link", "secondary", "secondary_link",
                  "tertiary", "tertiary_link", "unclassified",
                  "residential", "track"),
  track_speed = 30,             # km/h on unpaved tracks (AccessMod convention)

  # --- SRCA Annual Report 2025 observed moments (report page refs) ---
  obs = list(
    queue_min        = 48 / 60,   # call queuing            (p72, Table 17/2)
    mobilize_min     = 27 / 60,   # crew mobilization       (p72, Table 17/2)
    travel_scene_min = 7 + 46/60, # travel to scene, mean   (p72, Table 17/2)
    resp_mean_min    = 10 + 31/60,# overall response, mean  (p72, Table 18/2)
    resp8_share      = 45.43,     # % responses <= 8 min    (p32/46/67, KPI)
    stroke_resp_min  = 11,        # stroke-pathway response (p104)
    door_cath_min    = 37         # arrival-to-cath claim   (p99)
  ),

  # --- model scenarios and standards ---
  scene_min   = 15,             # on-scene minutes (primary; 10/20 sensitivity)
  std = list(
    response  = 8,              # SRCA KPI / NFPA 1710
    trauma    = 60,             # golden-hour convention (Branas 2005;
                                #   causal caveat: Newgard 2010)
    stemi_door  = 90,           # AHA/ACC FMC-to-device
    stemi_wire  = 120,          # ESC 2023 diagnosis-to-wire
    stroke    = 60              # Adeoye 2014 convention
  ),
  thresholds  = c(30, 60, 90, 120),   # minutes, for coverage tables

  # --- destination-set toggles (open decisions; edit and rerun) ---
  trauma_set  = "general",      # "general": MoH secondary with General/
                                #   Emergency scope + all tertiary (default)
                                # "all":     every MoH secondary + tertiary
                                # "beds100": general subset with >= 100 beds
  alpha_sens  = c(0.9, 1.1),    # +/-10% speed-factor sensitivity (lights &
                                #   sirens literature: ~2 min urban / ~9 min
                                #   rural savings; Brown 2000, Petzall 2011)
  scene_sens  = c(10, 20)
)
FIXED_PRE <- PAR$obs$queue_min + PAR$obs$mobilize_min   # 1.25 min

cache <- function(name, expr) {
  path <- file.path(prc, paste0(name, ".rds"))
  if (file.exists(path)) return(readRDS(path))
  val <- force(expr)
  saveRDS(val, path)
  val
}

# ---- 1. Supply: SRCA centres and condition-specific destinations ---------

srca_c <- read_csv(file.path(raw, "SRCA Centers/srca_centers.csv"),
                   show_col_types = FALSE) |>
  filter(!is.na(lat), !is.na(lon))
srca_sf <- st_as_sf(srca_c, coords = c("lon", "lat"), crs = 4326) |>
  st_transform(PAR$crs_m)

cath <- read_csv(file.path(raw, "Cath Centers/cath_centers.csv"),
                 show_col_types = FALSE) |>
  filter(!is.na(lat), !is.na(lon))
cath_sf <- st_as_sf(cath, coords = c("lon", "lat"), crs = 4326) |>
  st_transform(PAR$crs_m)

stroke <- read_csv(file.path(raw, "Stroke Centers/stroke_centers.csv"),
                   show_col_types = FALSE) |>
  filter(!is.na(lat), !is.na(lon))
stroke_sf <- st_as_sf(stroke, coords = c("lon", "lat"), crs = 4326) |>
  st_transform(PAR$crs_m)

# Stroke Tier 2 (clinician-verified reperfusion-capable subset). Two ways to
# supply it, either of which activates the Tier-2 arm on the next run:
#   (a) add a reperfusion_capable column (TRUE/FALSE) to stroke_centers.csv
#   (b) place a stroke_centers_tier2.csv (same schema) beside it
tier2_path <- file.path(raw, "Stroke Centers/stroke_centers_tier2.csv")
stroke_t2_sf <- if ("reperfusion_capable" %in% names(stroke)) {
  filter(stroke_sf, reperfusion_capable %in% c(TRUE, "TRUE", "yes", "1"))
} else if (file.exists(tier2_path)) {
  read_csv(tier2_path, show_col_types = FALSE) |>
    filter(!is.na(lat), !is.na(lon)) |>
    st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
    st_transform(PAR$crs_m)
} else NULL
has_t2 <- !is.null(stroke_t2_sf) && nrow(stroke_t2_sf) > 0

# Trauma destinations from the MoH facility list (2,480 rows; PHC excluded).
# NOTE (scope statement for the paper): cath and stroke lists span sectors
# because SRCA publishes its receiving hospitals; no equivalent hospital-level
# trauma destination table exists, so the trauma arm is scoped to the MoH
# hospital network and stated as such in Methods.
prov <- read_csv(file.path(raw, "list_of_healthcare_providers.csv"),
                 show_col_types = FALSE) |>
  filter(!is.na(lat), !is.na(lon), type %in% c("2ry", "3ry")) |>
  mutate(is_general = type == "3ry" |
           (type == "2ry" & str_detect(scope, "General|Emergency")))
trauma <- switch(PAR$trauma_set,
  all     = prov,
  general = filter(prov, is_general),
  beds100 = filter(prov, is_general,
                   type == "3ry" | replace_na(bed_capacity >= 100, FALSE)),
  stop("unknown PAR$trauma_set"))
trauma_sf <- st_as_sf(trauma, coords = c("lon", "lat"), crs = 4326) |>
  st_transform(PAR$crs_m)

message(sprintf(paste0(
  "Supply: %d SRCA centres | %d trauma hospitals (set '%s') | ",
  "%d cath centres | %d stroke Tier-1%s"),
  nrow(srca_sf), nrow(trauma_sf), PAR$trauma_set, nrow(cath_sf),
  nrow(stroke_sf),
  if (has_t2) sprintf(" | %d stroke Tier-2", nrow(stroke_t2_sf)) else
    " | Tier-2 not supplied (skipped)"))

# ---- 2. Governorates, regions, census population (demand) ----------------

gov <- st_read(file.path(raw, "governorate/Governorate.gpkg"), quiet = TRUE) |>
  st_transform(PAR$crs_m)
reg <- st_read(file.path(raw, "Regions/Regions.shp"), quiet = TRUE) |>
  st_transform(PAR$crs_m)

norm_ar <- function(x) {   # GASTAT <-> shapefile Arabic name matching
  x |>
    str_remove_all("\\(.*?\\)") |>
    str_replace_all("[أإآ]", "ا") |>
    str_replace_all("ة", "ه") |>
    str_replace_all("ى", "ي") |>
    str_remove_all("[ً-ْ]") |>
    str_remove_all("\\s")
}
gov$gov_key <- norm_ar(gov$Gov_AR)

# Region keys need one extra rule: the census writes Eastern Province as
# "المنطقة الشرقية" while the MoH provider file writes "الشرقية"; strip the
# generic leading "المنطقه" (post-norm_ar form) so the two sources match.
norm_reg <- function(x) norm_ar(x) |> str_remove("^المنطقه")

pop_city <- read_csv(
  file.path(raw, "PopulationbyNationalitybyRegionGovernorateCityandNationalityARCSV.csv"),
  show_col_types = FALSE)
names(pop_city) <- c("nationality", "governorate", "city", "region",
                     "gender", "pop")

pop_gov <- pop_city |>
  mutate(gov_key = norm_ar(governorate),
         saudi   = nationality == "سعودي") |>
  group_by(gov_key) |>
  summarise(pop_total = sum(pop),
            pop_saudi = sum(pop[saudi]), .groups = "drop")

# Governorate -> administrative region lookup (census), English names via the
# MoH provider file's clean region_ar/region_en pairs.
reg_names <- prov |> distinct(region_ar, region_en) |>
  mutate(reg_key = norm_reg(region_ar)) |> select(reg_key, region_en)
gov_region <- pop_city |>
  mutate(gov_key = norm_ar(governorate), reg_key = norm_reg(region)) |>
  count(gov_key, reg_key) |>
  group_by(gov_key) |> slice_max(n, n = 1, with_ties = FALSE) |> ungroup() |>
  left_join(reg_names, by = "reg_key") |> select(gov_key, region_en)

gov_pop <- gov |>
  left_join(pop_gov,    by = "gov_key") |>
  left_join(gov_region, by = "gov_key")

unmatched <- gov_pop |> st_drop_geometry() |>
  filter(is.na(pop_total)) |> select(Gov_AR, Gov_EN, gov_key)
if (nrow(unmatched) > 0) {
  write_csv(unmatched, file.path(prc, "unmatched_governorates.csv"))
  warning(sprintf(
    "%d governorates unmatched to census names -> data/processed/unmatched_governorates.csv",
    nrow(unmatched)))
}
reg_missing <- gov_pop |> st_drop_geometry() |>
  filter(is.na(region_en), !is.na(pop_total)) |> select(Gov_EN, gov_key)
if (nrow(reg_missing) > 0) {
  warning(sprintf(
    "%d governorates have census population but no region_en (region join miss): %s",
    nrow(reg_missing), paste(reg_missing$Gov_EN, collapse = ", ")))
}

# ---- 3. Gridded population -----------------------------------------------
# WorldPop 2020 constrained (100 m, settled areas only) aggregated to 1 km
# supplies the within-governorate distribution; cells rescaled so each
# governorate sums to its GASTAT 2022 census totals (total and Saudi).
# Denominator: TOTAL residents primary -- SRCA responds universally
# (opposite of the cluster-enrolment paper); Saudi-only is a sensitivity.

wp_path <- file.path(raw, "sau_ppp_2020_constrained.tif")
if (!file.exists(wp_path))
  stop("WorldPop raster not found at raw/sau_ppp_2020_constrained.tif")

grid <- cache("grid_population", {
  r  <- terra::aggregate(terra::rast(wp_path), fact = 10, fun = "sum",
                         na.rm = TRUE)   # 100 m -> 1 km
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "wp"
  df <- df[df$wp > 0, ]
  st_as_sf(df, coords = c("x", "y"), crs = 4326) |>
    st_transform(PAR$crs_m) |>
    st_join(gov_pop[, c("Gov_ID", "Gov_EN", "gov_key", "region_en",
                        "pop_total", "pop_saudi")],
            join = st_within) |>
    filter(!is.na(Gov_ID)) |>
    group_by(Gov_ID) |>
    mutate(w = wp / sum(wp)) |>
    ungroup() |>
    mutate(cell_total = w * pop_total,
           cell_saudi = w * pop_saudi,
           cell_id    = row_number()) |>
    select(cell_id, Gov_ID, Gov_EN, region_en, cell_total, cell_saudi)
})
message(sprintf("Population grid: %s cells | residents %.0f | Saudi %.0f",
                format(nrow(grid), big.mark = ","),
                sum(grid$cell_total), sum(grid$cell_saudi)))

# ---- 4. Travel times (dodgr street-network routing) ----------------------
# One-time input: the Geofabrik Gulf extract at PAR$pbf_path
#   https://download.geofabrik.de/asia/gcc-states-latest.osm.pbf  (~250 MB)
# Routing runs entirely in R (dodgr): car travel times from OSM road
# geometry, the 'motorcar' weighting profile, and maxspeed/oneway tags.
# No server, no Docker, no Java. The ~1 GB derived graph is cached beside
# the pbf after the first run and loads in seconds on reruns.

if (!file.exists(path.expand(PAR$pbf_path)))
  stop("Road extract not found at ", PAR$pbf_path,
       " - download it (see Section 4 header), then rerun.")

graph_path <- file.path(dirname(path.expand(PAR$pbf_path)), "dodgr_graph.rds")
graph <- if (file.exists(graph_path)) readRDS(graph_path) else {
  roads <- cache("roads_saudi", {
    q <- sprintf("SELECT * FROM lines WHERE highway IN (%s)",
                 paste(sprintf("'%s'", PAR$hw_keep), collapse = ", "))
    rd <- osmextract::oe_read(path.expand(PAR$pbf_path), layer = "lines",
                              extra_tags = c("oneway", "maxspeed", "lanes"),
                              query = q, quiet = TRUE)
    ksa <- st_union(gov) |> st_buffer(20000) |> st_transform(4326)
    st_filter(rd, ksa)
  })
  # Motorcar profile with unpaved tracks passable at PAR$track_speed. Must go
  # through a profile FILE: a profile data.frame makes dodgr skip travel-time
  # computation ("graph has no time column").
  wpf <- file.path(prc, "wt_profiles.json")
  if (!file.exists(wpf)) {
    dodgr::write_dodgr_wt_profile(file.path(prc, "wt_profiles"))
    j <- jsonlite::fromJSON(wpf)
    i <- j$weighting_profiles$name == "motorcar" &
         j$weighting_profiles$way  == "track"
    j$weighting_profiles$value[i]     <- 0.5
    j$weighting_profiles$max_speed[i] <- PAR$track_speed
    jsonlite::write_json(j, wpf, pretty = TRUE)
  }
  g <- dodgr::weight_streetnet(roads, wt_profile = "motorcar",
                               wt_profile_file = wpf)
  # Keep every component with >= 100 edges (inhabited islands, e.g. Farasan),
  # not just the largest; cross-component routes return NA (unreachable).
  keep <- names(which(table(g$component) >= 100))
  g <- g[g$component %in% as.integer(keep), ]
  saveRDS(g, graph_path)
  g
}
verts   <- dodgr::dodgr_vertices(graph)
vert_xy <- as.matrix(verts[, c("x", "y")])

snap_verts <- function(x, label) {
  ll <- st_coordinates(st_transform(x, 4326))
  nn <- RANN::nn2(vert_xy, ll, k = 1)
  d_km <- nn$nn.dists[, 1] * 111
  message(sprintf("  snap %-8s: median %.2f km, max %.1f km, >5 km: %d of %d",
                  label, median(d_km), max(d_km), sum(d_km > 5), nrow(x)))
  verts$id[nn$nn.idx[, 1]]
}

# Travel time between each origin and its nearest destination (minutes).
# Euclidean k-NN pruning, then chunked dodgr many-to-many time matrices.
# direction: "to_dest"  routes origin -> destination (scene -> hospital);
#            "from_dest" routes destination -> origin (station -> scene) --
#            the two differ on one-way streets, so ambulance legs use the
#            true driving direction.
tt_nearest <- function(origins, dests, direction = c("to_dest", "from_dest"),
                       k = PAR$knn_candidates, chunk = 1000, label = "dests") {
  direction <- match.arg(direction)
  k   <- min(k, nrow(dests))
  nn  <- RANN::nn2(st_coordinates(dests), st_coordinates(origins), k = k)
  o_v <- snap_verts(origins, "origins")
  d_v <- snap_verts(dests, label)
  out_t <- rep(NA_real_, nrow(origins))
  out_i <- rep(NA_integer_, nrow(origins))
  starts <- seq(1, nrow(origins), by = chunk)
  for (s in starts) {
    idx  <- s:min(s + chunk - 1, nrow(origins))
    cand <- sort(unique(as.vector(nn$nn.idx[idx, , drop = FALSE])))
    tmat <- if (direction == "to_dest") {
      dodgr::dodgr_times(graph, from = unique(o_v[idx]),
                         to = unique(d_v[cand])) / 60
    } else {
      dodgr::dodgr_times(graph, from = unique(d_v[cand]),
                         to = unique(o_v[idx])) / 60
    }
    for (j in seq_along(idx)) {
      tj <- if (direction == "to_dest") {
        tmat[match(o_v[idx[j]], rownames(tmat)),
             match(d_v[nn$nn.idx[idx[j], ]], colnames(tmat))]
      } else {
        tmat[match(d_v[nn$nn.idx[idx[j], ]], rownames(tmat)),
             match(o_v[idx[j]], colnames(tmat))]
      }
      tj[!is.finite(tj)] <- NA
      if (all(is.na(tj))) next
      b <- which.min(tj)
      out_t[idx[j]] <- tj[b]
      out_i[idx[j]] <- nn$nn.idx[idx[j], b]
    }
    message(sprintf("  ...%d / %d origins", max(idx), nrow(origins)))
  }
  list(minutes = out_t, dest_row = out_i)
}

tt <- cache("travel_times_grid", {
  out <- list(
    legA   = tt_nearest(grid, srca_sf, direction = "from_dest",
                        label = "srca"),
    trauma = tt_nearest(grid, trauma_sf, label = "trauma"),
    cath   = tt_nearest(grid, cath_sf,   label = "cath"),
    strk1  = tt_nearest(grid, stroke_sf, label = "stroke1")
  )
  if (has_t2) out$strk2 <- tt_nearest(grid, stroke_t2_sf, label = "stroke2")
  out
})
if (has_t2 && is.null(tt$strk2)) {
  tt$strk2 <- tt_nearest(grid, stroke_t2_sf, label = "stroke2")
  saveRDS(tt, file.path(prc, "travel_times_grid.rds"))
}

# Destination-set sensitivity leg: >=100-bed general/tertiary hospitals only
# (half the default trauma set is <100 beds; reviewers will question small
# general hospitals as trauma destinations, so this bounds the effect).
trauma_b100_sf <- prov |>
  filter(is_general, type == "3ry" | replace_na(bed_capacity >= 100, FALSE)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_transform(PAR$crs_m)
if (is.null(tt$trauma_b100)) {
  message(sprintf("Routing trauma >=100-bed sensitivity leg (n=%d)...",
                  nrow(trauma_b100_sf)))
  tt$trauma_b100 <- tt_nearest(grid, trauma_b100_sf, label = "trauma_b100")
  saveRDS(tt, file.path(prc, "travel_times_grid.rds"))
}

acc <- grid |>
  mutate(t_legA   = tt$legA$minutes,   srca_row   = tt$legA$dest_row,
         t_trauma = tt$trauma$minutes, trauma_row = tt$trauma$dest_row,
         t_cath   = tt$cath$minutes,   cath_row   = tt$cath$dest_row,
         t_strk1  = tt$strk1$minutes,  strk1_row  = tt$strk1$dest_row,
         t_strk2  = if (has_t2) tt$strk2$minutes else NA_real_)

# Population-weighted summary helpers
wq <- function(x, w, p) {
  ok <- !is.na(x) & !is.na(w); x <- x[ok]; w <- w[ok]
  if (!length(x)) return(NA_real_)
  o <- order(x); x <- x[o]; w <- w[o]
  x[which.max(cumsum(w) / sum(w) >= p)]
}
cov_pct <- function(x, w, thr) sum(w[!is.na(x) & x <= thr]) / sum(w) * 100

# ---- 5. Calibration and held-out validation ------------------------------
# One moment fitted: alpha scales unadjusted car travel times so the
# DEMAND-weighted mean of alpha*LegA equals the observed mean travel-to-
# scene (7:46). Demand weight = cell population (call volume proxy: SRCA
# means are call-weighted and therefore urban-skewed; an unweighted mean
# over cells would compare a population-weighted model to a call-weighted
# statistic). All other observed moments are held out as validation.

alpha <- PAR$obs$travel_scene_min /
  weighted.mean(acc$t_legA, acc$cell_total, na.rm = TRUE)
message(sprintf(
  "Calibration: alpha = %.3f (raw demand-weighted LegA mean %.2f min -> %.2f)",
  alpha, weighted.mean(acc$t_legA, acc$cell_total, na.rm = TRUE),
  PAR$obs$travel_scene_min))

acc <- acc |> mutate(t_resp = FIXED_PRE + alpha * t_legA)

# Held-out checks. Note the report's own components (0:48 + 0:27 + 7:46 =
# 9:01) undershoot its overall mean response (10:31); the residual ~1.5 min
# is unmodelled dispatch variance -- reported, not fitted.
validation <- tibble(
  moment   = c("mean travel-to-scene (min) [FITTED]",
               "share of responses <= 8 min (%) [held out]",
               "stroke-pathway mean response (min) [held out]",
               "overall mean response (min) [held out]"),
  observed = c(PAR$obs$travel_scene_min, PAR$obs$resp8_share,
               PAR$obs$stroke_resp_min, PAR$obs$resp_mean_min),
  modelled = c(weighted.mean(alpha * acc$t_legA, acc$cell_total, na.rm = TRUE),
               cov_pct(acc$t_resp, acc$cell_total, PAR$std$response),
               weighted.mean(acc$t_resp, acc$cell_total, na.rm = TRUE),
               weighted.mean(acc$t_resp, acc$cell_total, na.rm = TRUE))
)
write_csv(validation, file.path(tabs, "validation_moments.csv"))
print(validation)

# ---- 6. Full-chain call-to-door times and standards compliance -----------

acc <- acc |>
  mutate(
    t_door_trauma = t_resp + PAR$scene_min + alpha * t_trauma,
    t_door_cath   = t_resp + PAR$scene_min + alpha * t_cath,
    t_door_strk1  = t_resp + PAR$scene_min + alpha * t_strk1,
    t_door_strk2  = t_resp + PAR$scene_min + alpha * t_strk2,
    # STEMI call-to-wire adds the report's observed 37-min arrival-to-cath
    t_wire_cath   = t_door_cath + PAR$obs$door_cath_min
  )

natl_row <- function(x, w, std, label) {
  tibble(condition = label,
         pop_covered_pct = sum(w[!is.na(x)]) / sum(w) * 100,
         mean   = weighted.mean(x, w, na.rm = TRUE),
         median = wq(x, w, .5),
         p90    = wq(x, w, .9),
         std_min = std,
         within_std_pct = cov_pct(x, w, std)) |>
    bind_cols(as_tibble(setNames(
      lapply(PAR$thresholds, function(th) cov_pct(x, w, th)),
      paste0("le_", PAR$thresholds, "min"))))
}

w <- acc$cell_total
natl <- bind_rows(
  natl_row(acc$t_resp,        w, PAR$std$response,  "Response (call to scene)"),
  natl_row(acc$t_door_trauma, w, PAR$std$trauma,    "Trauma: call to hospital door"),
  natl_row(acc$t_door_cath,   w, PAR$std$stemi_door,"STEMI: call to cath-centre door"),
  natl_row(acc$t_wire_cath,   w, PAR$std$stemi_wire,"STEMI: call to wire (door + 37 min)"),
  natl_row(acc$t_door_strk1,  w, PAR$std$stroke,    "Stroke: call to receiving-hospital door (Tier 1)"),
  if (has_t2) natl_row(acc$t_door_strk2, w, PAR$std$stroke,
                       "Stroke: call to reperfusion-capable door (Tier 2)")
)
write_csv(natl, file.path(tabs, "table1_national.csv"))
print(natl, width = Inf)

# Governorate table (Table 2 / supplement)
gov_tab <- acc |> st_drop_geometry() |>
  group_by(region_en, Gov_EN) |>
  summarise(
    pop            = sum(cell_total),
    resp_med       = wq(t_resp, cell_total, .5),
    resp_le8_pct   = cov_pct(t_resp, cell_total, PAR$std$response),
    trauma_med     = wq(t_door_trauma, cell_total, .5),
    trauma_std_pct = cov_pct(t_door_trauma, cell_total, PAR$std$trauma),
    stemi_med      = wq(t_door_cath, cell_total, .5),
    stemi_std_pct  = cov_pct(t_door_cath, cell_total, PAR$std$stemi_door),
    stroke_med     = wq(t_door_strk1, cell_total, .5),
    stroke_std_pct = cov_pct(t_door_strk1, cell_total, PAR$std$stroke),
    .groups = "drop") |>
  arrange(region_en, desc(pop))
write_csv(gov_tab, file.path(tabs, "table2_governorate.csv"))

reg_tab <- acc |> st_drop_geometry() |>
  group_by(region_en) |>
  summarise(
    pop            = sum(cell_total),
    resp_le8_pct   = cov_pct(t_resp, cell_total, PAR$std$response),
    trauma_std_pct = cov_pct(t_door_trauma, cell_total, PAR$std$trauma),
    stemi_std_pct  = cov_pct(t_door_cath, cell_total, PAR$std$stemi_door),
    stroke_std_pct = cov_pct(t_door_strk1, cell_total, PAR$std$stroke),
    .groups = "drop") |> arrange(desc(pop))

# SRCA CAD regional aggregates (Annual Report 2025 Tables 22/2, 23/2, 28/2;
# provenance in raw/srca_report_regional_2025.csv header). The report has no
# per-region response-time table, so these support demand-vs-access summaries,
# not regional validation. Checksums guard against transcription drift.
srca_reg <- read_csv(file.path(raw, "srca_report_regional_2025.csv"),
                     comment = "#", show_col_types = FALSE)
stopifnot(sum(srca_reg$missions_2025)      == 1416037,
          sum(srca_reg$transports_2025)    == 568827,
          sum(srca_reg$launch_points_2025) == 520,
          all(reg_tab$region_en %in% srca_reg$region_en))
reg_tab <- reg_tab |>
  left_join(srca_reg, by = "region_en") |>
  mutate(pop_share_pct     = pop / sum(pop) * 100,
         mission_share_pct = missions_2025 / sum(missions_2025) * 100,
         missions_per_1000 = missions_2025 / pop * 1000)

# Demand-weighted national coverage: reweight regional shares by observed 2025
# mission volume instead of residential population. Population weighting over
# a partition reproduces the cell-level national shares exactly, so the
# population-weighted row doubles as an internal-consistency check.
wt_row <- function(wcol, label) {
  w <- reg_tab[[wcol]]   # capture before summarise() masks the columns
  reg_tab |>
    summarise(region_en = label,
              across(c(resp_le8_pct, trauma_std_pct, stemi_std_pct,
                       stroke_std_pct), ~ weighted.mean(.x, w)),
              pop = sum(pop),
              missions_2025 = sum(missions_2025),
              transports_2025 = sum(transports_2025),
              launch_points_2025 = sum(launch_points_2025),
              pop_share_pct = 100, mission_share_pct = 100,
              missions_per_1000 = sum(missions_2025) / sum(pop) * 1000)
}
reg_tab <- bind_rows(reg_tab,
                     wt_row("pop",           "National (population-weighted)"),
                     wt_row("missions_2025", "National (mission-weighted)"))
write_csv(reg_tab, file.path(tabs, "tableS_regional.csv"))
print(reg_tab |> select(region_en, resp_le8_pct, missions_per_1000,
                        pop_share_pct, mission_share_pct), n = 15)

# Urban-rural gradient. Degree-of-Urbanisation density conventions on 1-km
# cells (EC/GHSL DEGURBA): >= 1,500 residents/km2 urban centre; 300-1,500
# urban cluster; < 300 rural. Ties the access deficit to GASTAT's finding
# that 39.7% of serious road-traffic accidents occur outside cities.
acc <- acc |>
  mutate(dens_class = cut(cell_total, c(-Inf, 300, 1500, Inf),
    labels = c("Rural (<300 per km2)", "Peri-urban (300-1,500 per km2)",
               "Urban (>=1,500 per km2)")))
ur_tab <- acc |> st_drop_geometry() |>
  group_by(dens_class) |>
  summarise(
    pop            = sum(cell_total),
    resp_med       = wq(t_resp, cell_total, .5),
    resp_le8_pct   = cov_pct(t_resp, cell_total, PAR$std$response),
    trauma_med     = wq(t_door_trauma, cell_total, .5),
    trauma_std_pct = cov_pct(t_door_trauma, cell_total, PAR$std$trauma),
    stemi_std_pct  = cov_pct(t_door_cath, cell_total, PAR$std$stemi_door),
    stroke_std_pct = cov_pct(t_door_strk1, cell_total, PAR$std$stroke),
    .groups = "drop") |>
  mutate(pop_pct = pop / sum(pop) * 100, .after = pop)
write_csv(ur_tab, file.path(tabs, "tableS_urban_rural.csv"))
print(ur_tab)

# ---- 7. 8-minute KPI: structural vs operational decomposition ------------
# How much of the KPI shortfall is geography (station placement vs where
# people live) and how much is pre-departure process (queue + mobilization)?
#   structural ceiling: share reachable with alpha*LegA <= 8 (instant dispatch)
#   operational model : share with FIXED_PRE + alpha*LegA <= 8
#   observed KPI      : 45.43%
#
# Operational-delay scenario: the modelled KPI exceeds the observed one
# because unit unavailability, out-of-position starts, simultaneous calls,
# and dispatch variance are outside a geographic model. delta_op is the
# uniform extra pre-travel delay that reconciles the modelled response
# distribution with the observed KPI; door-time compliance recomputed under
# it (Section 10) is the operational LOWER bound, the base case the
# geographic UPPER bound.

delta_op <- uniroot(
  function(d) cov_pct(FIXED_PRE + d + alpha * acc$t_legA, w,
                      PAR$std$response) - PAR$obs$resp8_share,
  c(0, 30))$root
acc <- acc |> mutate(t_resp_op = t_resp + delta_op)
message(sprintf(
  "Operational delay reconciling modelled with observed KPI: %.2f min",
  delta_op))

# The last two rows are the operational lever: the 8-minute share if the
# implied per-call delay were cut by 0.5 or 1 minute (Results, Comparing Levers).
kpi <- tibble(
  quantity = c("Structural ceiling (travel only <= 8 min)",
               "Modelled KPI (with observed queue + mobilization)",
               "Observed SRCA KPI 2025",
               "Implied additional operational delay (min)",
               "Modelled KPI with operational delay reduced by 0.5 min",
               "Modelled KPI with operational delay reduced by 1 min"),
  value = c(cov_pct(alpha * acc$t_legA, w, PAR$std$response),
            cov_pct(acc$t_resp, w, PAR$std$response),
            PAR$obs$resp8_share,
            delta_op,
            cov_pct(acc$t_resp_op - 0.5, w, PAR$std$response),
            cov_pct(acc$t_resp_op - 1, w, PAR$std$response)),
  unit = c("percent", "percent", "percent", "minutes", "percent", "percent"))
write_csv(kpi, file.path(tabs, "tableS_kpi_decomposition.csv"))
print(kpi)

# ---- 8. Bypass analysis (mothership vs nearest-ED penalty) ---------------
# Extra transport minutes to reach definitive care (cath / stroke centre)
# beyond the nearest trauma-set hospital: the geographic price of a
# mothership (direct transport) policy vs drip-and-ship.

acc <- acc |>
  mutate(bypass_cath  = alpha * (t_cath  - t_trauma),
         bypass_strk1 = alpha * (t_strk1 - t_trauma),
         bypass_strk2 = alpha * (t_strk2 - t_trauma))

bypass <- tibble(
  condition = c("STEMI (cath centre)", "Stroke Tier 1",
                if (has_t2) "Stroke Tier 2 (reperfusion-capable)"),
  median_extra_min = c(wq(acc$bypass_cath, w, .5),
                       wq(acc$bypass_strk1, w, .5),
                       if (has_t2) wq(acc$bypass_strk2, w, .5)),
  p90_extra_min    = c(wq(acc$bypass_cath, w, .9),
                       wq(acc$bypass_strk1, w, .9),
                       if (has_t2) wq(acc$bypass_strk2, w, .9)),
  extra_le15_pct   = c(cov_pct(acc$bypass_cath, w, 15),
                       cov_pct(acc$bypass_strk1, w, 15),
                       if (has_t2) cov_pct(acc$bypass_strk2, w, 15)),
  extra_le30_pct   = c(cov_pct(acc$bypass_cath, w, 30),
                       cov_pct(acc$bypass_strk1, w, 30),
                       if (has_t2) cov_pct(acc$bypass_strk2, w, 30)))
write_csv(bypass, file.path(tabs, "tableS_bypass.csv"))

# ---- 9. Inequality: weighted Gini and Theil decomposition ----------------

gini_w <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0; x <- x[ok]; w <- w[ok]
  o <- order(x); x <- x[o]; w <- w[o]
  p <- cumsum(w) / sum(w)
  L <- cumsum(x * w) / sum(x * w)
  sum(diff(c(0, p)) * (L + c(0, head(L, -1)))) |> (\(a) 1 - a)()
}
theil_w <- function(x, w) {   # Theil-T; x > 0
  ok <- !is.na(x) & !is.na(w) & w > 0 & x > 0; x <- x[ok]; w <- w[ok]
  s <- w / sum(w); mu <- sum(s * x)
  sum(s * (x / mu) * log(x / mu))
}
theil_decomp <- function(x, w, g) {
  ok <- !is.na(x) & !is.na(w) & w > 0 & x > 0 & !is.na(g)
  x <- x[ok]; w <- w[ok]; g <- g[ok]
  s  <- w / sum(w); mu <- sum(s * x)
  by_g <- tibble(x, w, g) |> group_by(g) |>
    summarise(sw = sum(w), mug = weighted.mean(x, w),
              Tg = theil_w(x, w), .groups = "drop") |>
    mutate(sh = sw / sum(sw))
  between <- with(by_g, sum(sh * (mug / mu) * log(mug / mu)))
  within  <- with(by_g, sum(sh * (mug / mu) * Tg))
  tibble(total = between + within, between = between, within = within,
         between_pct = between / (between + within) * 100)
}

eq_tab <- bind_rows(lapply(
  list("Response"       = acc$t_resp,
       "Trauma door"    = acc$t_door_trauma,
       "STEMI door"     = acc$t_door_cath,
       "Stroke door T1" = acc$t_door_strk1) |>
    (\(l) if (has_t2) c(l, list("Stroke door T2" = acc$t_door_strk2)) else l)(),
  function(x) {
    td_reg <- theil_decomp(x, w, acc$region_en)
    td_gov <- theil_decomp(x, w, acc$Gov_EN)
    tibble(gini               = gini_w(x, w),
           theil_total        = td_reg$total,
           between_region_pct = td_reg$between_pct,
           between_gov_pct    = td_gov$between_pct)
  }),
  .id = "metric")
write_csv(eq_tab, file.path(tabs, "tableS_inequality.csv"))
print(eq_tab)

# ---- 9b. Spatial clustering of access deficits (the paper's one test) ----
# H1: governorate-level compliance is spatially clustered (a contiguous
# deficit belt), not randomly arranged. Permutation inference is the only
# valid frame here: modelled compliance values are population quantities
# (a census of cells, nothing sampled), so classical sampling-based
# p-values are undefined; the null that IS well-defined is random spatial
# rearrangement of governorate values. Global Moran's I (queen contiguity,
# row-standardised weights, 9,999 permutations, one-sided for positive
# autocorrelation); LISA Low-Low clusters (BH-adjusted across 150
# governorates) name the belt. Governorates enter unweighted (areal units).

gov_sp <- gov |> select(Gov_EN) |> inner_join(gov_tab, by = "Gov_EN")
message(sprintf("Spatial test frame: %d of %d governorates matched",
                nrow(gov_sp), nrow(gov_tab)))
set.seed(2026)
nb <- spdep::poly2nb(gov_sp, queen = TRUE)
lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
sp_metrics <- c(resp_le8_pct   = "Response <= 8 min",
                trauma_std_pct = "Trauma door <= 60 min",
                stemi_std_pct  = "STEMI door <= 90 min",
                stroke_std_pct = "Stroke door <= 60 min")

sp_tab <- bind_rows(lapply(names(sp_metrics), function(v) {
  mc <- spdep::moran.mc(gov_sp[[v]], lw, nsim = 9999, zero.policy = TRUE)
  tibble(metric = sp_metrics[[v]], morans_I = unname(mc$statistic),
         p_perm = mc$p.value)
}))
write_csv(sp_tab, file.path(tabs, "tableS_spatial_moran.csv"))
print(sp_tab)

lisa_lowlow <- bind_rows(lapply(names(sp_metrics), function(v) {
  x    <- gov_sp[[v]]
  lm_  <- spdep::localmoran_perm(x, lw, nsim = 9999, zero.policy = TRUE)
  z    <- x - mean(x)
  lagz <- spdep::lag.listw(lw, z, zero.policy = TRUE)
  padj <- p.adjust(lm_[, "Pr(folded) Sim"], "BH")
  keep <- which(z < 0 & lagz < 0 & padj < 0.05)
  if (!length(keep)) return(NULL)
  tibble(metric         = sp_metrics[[v]],
         region_en      = gov_sp$region_en[keep],
         governorate    = gov_sp$Gov_EN[keep],
         pop            = gov_sp$pop[keep],
         compliance_pct = x[keep])
}))
write_csv(lisa_lowlow, file.path(tabs, "tableS_lisa_lowlow.csv"))
message(sprintf("LISA Low-Low (deficit belt) governorate-metric pairs: %d",
                nrow(lisa_lowlow)))

# ---- 10. Sensitivity analyses --------------------------------------------
# Grid of scenarios x conditions, reporting % within the primary standard.

sens_row <- function(scenario, scene, a, wgt, t_tr, t_ca, t_s1,
                     extra = 0, b = a) {
  resp <- FIXED_PRE + extra + a * acc$t_legA
  tibble(scenario = scenario,
         trauma_std_pct = cov_pct(resp + scene + b * t_tr, wgt, PAR$std$trauma),
         stemi_std_pct  = cov_pct(resp + scene + b * t_ca, wgt, PAR$std$stemi_door),
         stroke_std_pct = cov_pct(resp + scene + b * t_s1, wgt, PAR$std$stroke),
         resp_le8_pct   = cov_pct(resp, wgt, PAR$std$response))
}

sens <- bind_rows(
  sens_row("Base (scene 15, alpha, total pop)", PAR$scene_min, alpha, w,
           acc$t_trauma, acc$t_cath, acc$t_strk1),
  sens_row("On-scene 10 min", PAR$scene_sens[1], alpha, w,
           acc$t_trauma, acc$t_cath, acc$t_strk1),
  sens_row("On-scene 20 min", PAR$scene_sens[2], alpha, w,
           acc$t_trauma, acc$t_cath, acc$t_strk1),
  sens_row("Speed +10% (alpha x 0.9)", PAR$scene_min, alpha * PAR$alpha_sens[1],
           w, acc$t_trauma, acc$t_cath, acc$t_strk1),
  sens_row("Speed -10% (alpha x 1.1)", PAR$scene_min, alpha * PAR$alpha_sens[2],
           w, acc$t_trauma, acc$t_cath, acc$t_strk1),
  sens_row("Saudi-only denominator", PAR$scene_min, alpha, acc$cell_saudi,
           acc$t_trauma, acc$t_cath, acc$t_strk1),
  sens_row(sprintf("Operational-delay scenario (+%.1f min, KPI-matched)",
                   delta_op), PAR$scene_min, alpha, w,
           acc$t_trauma, acc$t_cath, acc$t_strk1, extra = delta_op),
  sens_row("Leg B uncalibrated (scene-to-door speed factor 1.0)",
           PAR$scene_min, alpha, w,
           acc$t_trauma, acc$t_cath, acc$t_strk1, b = 1),
  sens_row(sprintf("Trauma: >=100-bed hospitals only (n=%d)",
                   nrow(trauma_b100_sf)), PAR$scene_min, alpha, w,
           tt$trauma_b100$minutes, acc$t_cath, acc$t_strk1),
  if (has_t2)
    sens_row("Stroke Tier 2 destinations", PAR$scene_min, alpha, w,
             acc$t_trauma, acc$t_cath, acc$t_strk2)
)
write_csv(sens, file.path(tabs, "tableS_sensitivity.csv"))
print(sens)

# ---- 11. Figures ---------------------------------------------------------

theme_map <- theme_void(base_size = 9) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 10, hjust = 0))
time_bins <- function(x) cut(x, c(0, 30, 60, 90, 120, Inf),
  labels = c("≤30", "31–60", "61–90", "91–120", ">120"),
  include.lowest = TRUE)
# RdYlBu (5-class, reversed): fast = blue, slow = red; CVD-safe
bin_cols <- setNames(c("#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D7191C"),
                     c("≤30", "31–60", "61–90",
                       "91–120", ">120"))

gov_bg  <- st_geometry(gov)

# Display layer: 10-km blocks (pop-weighted mean), projected for legibility
crs_disp <- 32638                                     # UTM 38N
gov_disp <- st_transform(gov_bg, crs_disp)
xy_disp  <- st_coordinates(st_transform(acc, crs_disp))
blk <- acc |>
  st_drop_geometry() |>
  mutate(bx = floor(xy_disp[, 1] / 1e4), by = floor(xy_disp[, 2] / 1e4)) |>
  group_by(bx, by) |>
  summarise(
    t_door_trauma = weighted.mean(t_door_trauma, cell_total, na.rm = TRUE),
    t_door_cath   = weighted.mean(t_door_cath,   cell_total, na.rm = TRUE),
    t_door_strk1  = weighted.mean(t_door_strk1,  cell_total, na.rm = TRUE),
    .groups = "drop") |>
  mutate(x = bx * 1e4 + 5e3, y = by * 1e4 + 5e3)

cities <- tibble::tribble(
  ~name,     ~lon,  ~lat,
  "Riyadh",  46.71, 24.63,
  "Jeddah",  39.19, 21.49,
  "Makkah",  39.83, 21.42,
  "Madinah", 39.61, 24.47,
  "Dammam",  50.10, 26.43,
  "Abha",    42.51, 18.22,
  "Tabuk",   36.58, 28.38) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_transform(crs_disp)
cities <- cities |>
  mutate(x = st_coordinates(cities)[, 1], y = st_coordinates(cities)[, 2])

panel_map <- function(tcol, fac_layer, title) {
  d <- blk |> mutate(bin = time_bins(.data[[tcol]])) |> filter(!is.na(bin))
  ggplot() +
    geom_sf(data = gov_disp, fill = "grey92", color = "grey78",
            linewidth = 0.12) +
    geom_tile(data = d, aes(x, y, fill = bin), width = 1e4, height = 1e4) +
    geom_sf(data = st_transform(fac_layer, crs_disp), shape = 3, size = 0.7,
            color = "black", stroke = 0.5) +
    geom_text(data = st_drop_geometry(cities), aes(x, y, label = name),
              size = 2.3, color = "grey10", fontface = "bold", vjust = -0.7) +
    scale_fill_manual(values = bin_cols, name = "Call-to-door (min)",
                      drop = FALSE) +
    labs(title = title) + theme_map
}

fig1 <- panel_map("t_door_trauma", trauma_sf,
                  sprintf("A  Trauma (n=%d hospitals)", nrow(trauma_sf))) +
        panel_map("t_door_cath", cath_sf,
                  sprintf("B  STEMI (n=%d cath centres)", nrow(cath_sf))) +
        panel_map("t_door_strk1", stroke_sf,
                  sprintf("C  Stroke Tier 1 (n=%d)", nrow(stroke_sf))) +
        patchwork::plot_layout(ncol = 3, guides = "collect") &
        theme(legend.position = "bottom")
ggsave(file.path(figs, "fig1_maps_call_to_door.png"), fig1,
       width = 13, height = 5.6, dpi = 300, bg = "white")

# Fig 2 -- calibration/validation + KPI decomposition. Panel A: one facet
# per moment (own y scale, so minutes and percent never share an axis).
f2a <- validation |>
  mutate(moment = str_wrap(moment, 22)) |>
  pivot_longer(c(observed, modelled)) |>
  mutate(name = factor(name, c("observed", "modelled"))) |>
  ggplot(aes(name, value, fill = name)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f", value)), vjust = -0.35, size = 2.3) +
  facet_wrap(~moment, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c(observed = "grey35", modelled = "steelblue")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "A  Calibration (1 fitted moment) and held-out validation",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(strip.text = element_text(size = 6.6),
        axis.text.x = element_text(size = 7.5))

resp_grid <- tibble(t = seq(0, 30, 0.25)) |>
  rowwise() |>
  mutate(model  = cov_pct(acc$t_resp, w, t),
         travel = cov_pct(alpha * acc$t_legA, w, t),
         oper   = cov_pct(acc$t_resp_op, w, t)) |>
  ungroup()
lt <- c("Travel only (structural ceiling)"   = 2,
        "With queue + mobilization"          = 1,
        "With operational delay (KPI-matched)" = 3)
f2b <- ggplot(resp_grid, aes(t)) +
  geom_line(aes(y = travel, linetype = "Travel only (structural ceiling)")) +
  geom_line(aes(y = model,  linetype = "With queue + mobilization")) +
  geom_line(aes(y = oper,
                linetype = "With operational delay (KPI-matched)")) +
  geom_vline(xintercept = 8, color = "red3", linewidth = 0.3) +
  geom_hline(yintercept = PAR$obs$resp8_share, color = "grey40",
             linetype = "dotted") +
  annotate("text", x = 29, y = PAR$obs$resp8_share - 4, hjust = 1, size = 2.6,
           label = sprintf("Observed KPI %.1f%%", PAR$obs$resp8_share)) +
  scale_linetype_manual(values = lt, breaks = names(lt), name = NULL) +
  labs(title = "B  Population within t minutes of ambulance response",
       x = "minutes from call", y = "% of population") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom") +
  guides(linetype = guide_legend(nrow = 3))
ggsave(file.path(figs, "fig2_calibration_kpi.png"), f2a + f2b,
       width = 10, height = 4.4, dpi = 300, bg = "white")

# Fig 3 -- Lorenz curves (population share vs cumulative time burden)
lorenz_df <- function(x, wgt, label) {
  ok <- !is.na(x) & wgt > 0
  x <- x[ok]; wgt <- wgt[ok]; o <- order(x)
  tibble(p = cumsum(wgt[o]) / sum(wgt),
         L = cumsum(x[o] * wgt[o]) / sum(x * wgt),
         metric = sprintf("%s (Gini %.2f)", label, gini_w(x, wgt)))
}
lor <- bind_rows(
  lorenz_df(acc$t_resp,        w, "Response"),
  lorenz_df(acc$t_door_trauma, w, "Trauma"),
  lorenz_df(acc$t_door_cath,   w, "STEMI"),
  lorenz_df(acc$t_door_strk1,  w, "Stroke T1"),
  if (has_t2) lorenz_df(acc$t_door_strk2, w, "Stroke T2"))
fig3 <- ggplot(lor, aes(p, L, color = metric)) +
  geom_abline(slope = 1, linetype = 3, color = "grey60") +
  geom_line(linewidth = 0.6) +
  scale_color_viridis_d(end = 0.85, name = NULL) +
  coord_equal() +
  labs(x = "cumulative population share (fastest first)",
       y = "cumulative share of prehospital time",
       title = "Concentration of prehospital time burden") +
  theme_minimal(base_size = 9) + theme(legend.position = "bottom")
ggsave(file.path(figs, "fig3_lorenz.png"), fig3,
       width = 5.5, height = 5.8, dpi = 300, bg = "white")

# Fig 4 (supplement) -- governorate choropleth, median call-to-door by cond.
gov_map <- gov |>
  left_join(gov_tab |> select(Gov_EN, trauma_med, stemi_med, stroke_med),
            by = "Gov_EN")
f4p <- function(col, title) {
  ggplot(gov_map) +
    geom_sf(aes(fill = .data[[col]]), color = "grey70", linewidth = 0.1) +
    scale_fill_viridis_c(direction = -1, name = "median min", na.value = "white") +
    labs(title = title) + theme_map
}
fig4 <- f4p("trauma_med", "A  Trauma") + f4p("stemi_med", "B  STEMI") +
        f4p("stroke_med", "C  Stroke T1") +
        patchwork::plot_layout(ncol = 3)
ggsave(file.path(figs, "figS_governorate_choropleth.png"), fig4,
       width = 12, height = 4.8, dpi = 300, bg = "white")

# Fig S2 (supplement) -- facility locations over region/governorate boundaries;
# receiving facilities colored by sector (MoH network = Health Holding Company
# clusters, private, university, military, National Guard, other government).
reg_disp <- st_transform(st_geometry(reg), crs_disp)
sector_cols <- c("MoH" = "#0072B2", "Private" = "#E69F00",
                 "University" = "#009E73", "National Guard" = "#CC79A7",
                 "Other government" = "#56B4E9", "Military" = "#D55E00")
fac_panel <- function(layer, title, by_sector = FALSE, pt = 0.7, col = "#0072B2",
                      legend = TRUE) {
  p <- ggplot() +
    geom_sf(data = gov_disp, fill = "grey96", color = "grey82", linewidth = 0.1) +
    geom_sf(data = reg_disp, fill = NA, color = "grey35", linewidth = 0.35)
  d <- st_transform(layer, crs_disp)
  if (by_sector) {
    d$sector <- factor(d$sector, levels = names(sector_cols))
    p <- p + geom_sf(data = d, aes(color = sector), size = 1.1) +
      scale_color_manual(values = sector_cols, name = "Sector", drop = FALSE,
                         guide = if (legend) "legend" else "none")
  } else {
    p <- p + geom_sf(data = d, color = col, size = pt)
  }
  p + labs(title = title) + theme_map
}
figS2 <- (fac_panel(srca_sf, sprintf("A  SRCA ambulance stations (n=%d)", nrow(srca_sf)),
                    pt = 0.4, col = "grey25") +
          fac_panel(trauma_sf, sprintf("B  MoH general hospitals (n=%d)", nrow(trauma_sf)))) /
         (fac_panel(cath_sf, sprintf("C  Catheterization centers (n=%d)", nrow(cath_sf)),
                    by_sector = TRUE) +
          fac_panel(stroke_sf, sprintf("D  Stroke receiving hospitals (n=%d)", nrow(stroke_sf)),
                    by_sector = TRUE, legend = FALSE)) +
         patchwork::plot_layout(guides = "collect") &
         theme(legend.position = "bottom")
ggsave(file.path(figs, "figS_facility_maps.png"), figS2,
       width = 11, height = 9.5, dpi = 300, bg = "white")

# Cell-level export for the supplement / reproducibility deposit
acc |> st_drop_geometry() |>
  select(cell_id, Gov_EN, region_en, cell_total, cell_saudi, dens_class,
         t_legA, t_resp, t_resp_op, starts_with("t_door"),
         starts_with("bypass")) |>
  write_csv(file.path(prc, "cell_level_results.csv"))

# ---- 12. Decision analysis: leg decomposition + capability-upgrade sites -
# Added 2026-10-07. Turns the descriptive access surface into two decision
# products:
#   (a) tableS_transport_legs.csv -- population-weighted REACH leg (call to
#       scene) vs DELIVERY leg (scene to door) per condition. Reach is
#       operations-bound (dispatch/deployment, Section 7); delivery is
#       capability-location-bound, which motivates (b).
#   (b) tableS_upgrade_sites_stemi.csv / _stroke.csv -- greedy maximal-
#       coverage ranking of capability upgrades across the MoH general-
#       hospital network (operated by the Health Holding Company clusters;
#       cluster_en carried through from the provider list). Each round picks
#       the hospital whose upgrade newly covers the most residents within
#       the 90-min STEMI / 60-min stroke call-to-door window (base case).
#       STEMI candidates: the >=100-bed general/tertiary subset (cath-lab
#       feasibility), minus hospitals already on the SRCA cath list
#       (matched exactly by moh_facility_id). Stroke candidates: all
#       general hospitals minus those already on the SRCA stroke list.
# Routing for (b) is direction-true (cell -> hospital): the street graph is
# reversed by renaming its from_/to_ columns, so one Dijkstra per CANDIDATE
# (n=250) replaces one per cell (n=67,075). Cached sparse at <= 80 raw
# minutes; the widest per-cell budget is (90 - 15 - 1.25)/alpha = 78.7.

# (a) reach vs delivery legs
leg_row <- function(x, label) tibble(
  leg          = label,
  median_min   = wq(x, w, .5),
  mean_min     = weighted.mean(x, w, na.rm = TRUE),
  p90_min      = wq(x, w, .9),
  le_15min_pct = cov_pct(x, w, 15),
  le_30min_pct = cov_pct(x, w, 30),
  le_60min_pct = cov_pct(x, w, 60))
legs_tab <- bind_rows(
  leg_row(alpha * acc$t_legA,   "Drive to scene (calibrated, excl. pre-departure)"),
  leg_row(acc$t_resp,           "Reach: call to scene (1.25 min pre-departure + drive)"),
  leg_row(alpha * acc$t_trauma, "Deliver: scene to trauma door (transport only)"),
  leg_row(alpha * acc$t_cath,   "Deliver: scene to cath-centre door (transport only)"),
  leg_row(alpha * acc$t_strk1,  "Deliver: scene to stroke Tier-1 door (transport only)"))
write_csv(legs_tab, file.path(tabs, "tableS_transport_legs.csv"))
print(legs_tab)

# (b) candidate travel times (cell -> each general hospital), sparse cache.
# Candidate set is fixed to the full general-hospital network regardless of
# PAR$trauma_set, so the cache stays valid when the trauma toggle changes;
# an identity attribute guards against a stale cache anyway.
cand_sf <- prov |> filter(is_general) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_transform(PAR$crs_m)
cand_id <- as.character(cand_sf$facility_id)

cand_tt <- cache("legB_candidate_times", {
  grev <- graph
  flip <- c(from_id = "to_id",   from_lon = "to_lon",   from_lat = "to_lat",
            to_id   = "from_id", to_lon   = "from_lon", to_lat   = "from_lat")
  names(grev)[match(names(flip), names(grev))] <- unname(flip)
  cand_v <- snap_verts(cand_sf, "cand")
  cell_v <- snap_verts(grid, "cells")
  message(sprintf("  reverse-graph routing: %d candidates x %s cells...",
                  nrow(cand_sf), format(nrow(grid), big.mark = ",")))
  tm <- dodgr::dodgr_times(grev, from = unique(cand_v),
                           to = unique(cell_v)) / 60
  ri <- match(cand_v, rownames(tm))
  ci <- match(cell_v, colnames(tm))
  sp <- lapply(seq_len(nrow(cand_sf)), function(k) {
    tk <- tm[ri[k], ci]
    keep <- which(is.finite(tk) & tk <= 80)
    if (!length(keep)) return(NULL)
    data.frame(cand = k, cell = keep, t_raw = unname(tk[keep]))
  })
  res <- do.call(rbind, sp)
  attr(res, "cand_facility_id") <- as.character(cand_sf$facility_id)
  res
})
stopifnot(identical(attr(cand_tt, "cand_facility_id"), cand_id))

cath_moh   <- as.character(stats::na.omit(cath$moh_facility_id))
stroke_moh <- as.character(stats::na.omit(stroke$moh_facility_id))
b100_ids   <- as.character(trauma_b100_sf$facility_id)
km_to_nearest <- function(ex_sf)
  round(apply(st_distance(cand_sf, ex_sf), 1, min) / 1000, 1)
km_cath   <- km_to_nearest(cath_sf)
km_stroke <- km_to_nearest(stroke_sf)

greedy_upgrade <- function(idx, door_base, window, km_existing, top_n = 15,
                           min_gain = 10000) {
  budget  <- (window - PAR$scene_min - acc$t_resp) / alpha  # raw-minute budget
  covered <- !is.na(door_base) & door_base <= window
  base_pct <- sum(w[covered]) / sum(w) * 100
  el <- cand_tt[cand_tt$cand %in% idx, ]
  el <- el[!is.na(budget[el$cell]) & el$t_raw <= budget[el$cell], , drop = FALSE]
  rows <- list(tibble(
    rank = 0L, facility_id = NA_character_,
    hospital = "(baseline: current receiving network, no upgrade)",
    type = NA_character_, beds = NA_real_, governorate = NA_character_,
    region = NA_character_, cluster = NA_character_,
    km_to_nearest_existing = NA_real_, newly_covered_pop = 0,
    cum_within_pct = base_pct))
  cum <- base_pct
  for (r in seq_len(top_n)) {
    el <- el[!covered[el$cell], , drop = FALSE]
    if (!nrow(el)) break
    gain <- tapply(w[el$cell], el$cand, sum)
    g    <- max(gain)
    if (g < min_gain) break
    best  <- as.integer(names(gain)[which.max(gain)])
    newly <- unique(el$cell[el$cand == best])
    covered[newly] <- TRUE
    cum <- cum + g / sum(w) * 100
    rows[[r + 1]] <- tibble(
      rank = r, facility_id = cand_id[best],
      hospital = cand_sf$name_en[best], type = cand_sf$type[best],
      beds = cand_sf$bed_capacity[best],
      governorate = cand_sf$governorate_en[best],
      region = cand_sf$region_en[best], cluster = cand_sf$cluster_en[best],
      km_to_nearest_existing = km_existing[best],
      newly_covered_pop = g, cum_within_pct = cum)
  }
  bind_rows(rows)
}

idx_stemi  <- which(cand_id %in% b100_ids & !(cand_id %in% cath_moh))
idx_stroke <- which(!(cand_id %in% stroke_moh))
message(sprintf(
  "Upgrade candidates: STEMI %d (>=100-bed, non-cath) | stroke %d (general, non-receiving)",
  length(idx_stemi), length(idx_stroke)))

up_stemi  <- greedy_upgrade(idx_stemi,  acc$t_door_cath,  PAR$std$stemi_door, km_cath)
up_stroke <- greedy_upgrade(idx_stroke, acc$t_door_strk1, PAR$std$stroke,     km_stroke)
write_csv(up_stemi,  file.path(tabs, "tableS_upgrade_sites_stemi.csv"))
write_csv(up_stroke, file.path(tabs, "tableS_upgrade_sites_stroke.csv"))
print(up_stemi,  n = 20, width = Inf)
print(up_stroke, n = 20, width = Inf)

# (c) rank stability of the greedy shortlists under the +/-10% speed
# scenarios of Section 10 (alpha x 0.9 / x 1.1): response and transport legs
# both rescale; the candidate time table is in raw minutes and unchanged.
# Reports how much of each shortlist, and of its order, survives.
greedy_picks <- function(idx, door, window, a, t_resp, top_n = 15,
                         min_gain = 10000) {
  budget  <- (window - PAR$scene_min - t_resp) / a
  covered <- !is.na(door) & door <= window
  start   <- sum(w[covered]) / sum(w) * 100
  el <- cand_tt[cand_tt$cand %in% idx, ]
  el <- el[!is.na(budget[el$cell]) & el$t_raw <= budget[el$cell], , drop = FALSE]
  picks <- character(); cum <- start
  for (r in seq_len(top_n)) {
    el <- el[!covered[el$cell], , drop = FALSE]
    if (!nrow(el)) break
    gain <- tapply(w[el$cell], el$cand, sum)
    if (max(gain) < min_gain) break
    best <- as.integer(names(gain)[which.max(gain)])
    covered[unique(el$cell[el$cand == best])] <- TRUE
    cum <- cum + max(gain) / sum(w) * 100
    picks <- c(picks, cand_id[best])
  }
  list(picks = picks, start = start, cum = cum)
}
base_lists <- list(STEMI  = up_stemi$facility_id[-1],
                   Stroke = up_stroke$facility_id[-1])
stab_row <- function(cn, sc, run) {
  b <- base_lists[[cn]]; p <- run$picks
  tibble(condition = cn, scenario = sc,
         start_pct = run$start, cum_within_pct = run$cum,
         gain_pts = run$cum - run$start,
         top5_retained  = length(intersect(b[1:5], p[1:5])),
         top10_retained = length(intersect(b[1:10], p[1:10])),
         top15_retained = length(intersect(b, p)),
         same_rank1 = identical(p[1], b[1]),
         rank1_hospital = cand_sf$name_en[match(p[1], cand_id)])
}
scen <- list("Base" = 1,
             "Speed +10% (alpha x 0.9)" = PAR$alpha_sens[1],
             "Speed -10% (alpha x 1.1)" = PAR$alpha_sens[2])
stab <- bind_rows(lapply(names(scen), function(sc) {
  a <- alpha * scen[[sc]]
  t_resp_a <- FIXED_PRE + a * acc$t_legA
  bind_rows(
    stab_row("STEMI", sc, greedy_picks(idx_stemi,
      t_resp_a + PAR$scene_min + a * acc$t_cath, PAR$std$stemi_door, a, t_resp_a)),
    stab_row("Stroke", sc, greedy_picks(idx_stroke,
      t_resp_a + PAR$scene_min + a * acc$t_strk1, PAR$std$stroke, a, t_resp_a)))
}))
stopifnot(all(stab$top15_retained[stab$scenario == "Base"] == 15))
write_csv(stab, file.path(tabs, "tableS_upgrade_stability.csv"))
print(stab, width = Inf)

message("Done. Tables -> output/tables | Figures -> output/figures | ",
        sprintf("alpha = %.3f | trauma set '%s' (n=%d)%s",
                alpha, PAR$trauma_set, nrow(trauma_sf),
                if (has_t2) "" else " | stroke Tier 2 pending clinician list"))
