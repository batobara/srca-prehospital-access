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
#   5. Calibration (demand-weighted, one moment)
#   6. Full-chain call-to-door times and standards compliance
#   7. 8-minute response KPI: structural vs operational decomposition,
#      incl. a KPI-matched residual-delay scenario (lower bound)
#   8. Bypass analysis (mothership vs nearest-ED penalty)
#   9. Inequality in driving time (weighted Gini; nested Theil: regions,
#      governorates within regions, within governorates)
#  9b. Spatial clustering, descriptive (global Moran's I + LISA)
#  10. Sensitivity analyses (on-scene time, speed, denominator, >=100-bed
#      trauma set, Leg-B speed factor, residual delay, trauma tiers)
#  11. Tables and figures (5 main exhibits + supplement CSVs)
#  12. Decision analysis: reach vs delivery leg decomposition; greedy
#      capability-upgrade site rankings (STEMI cath, stroke receiving)
#      across the Health Holding Company (MoH) general-hospital network;
#      greedy new-station siting for the 8-minute response (structural
#      lever, compared with cutting the per-call residual delta)
#
# Model: call-to-door time for each 1-km populated cell =
#   fixed pre-departure intervals (queue 0:48 + mobilization 0:27, SRCA
#   2025 report, Table 17/2) + alpha * LegA + on-scene time (10/15/20 min
#   scenarios) + alpha * LegB, where alpha is a single speed factor
#   calibrated so the demand-weighted mean of alpha*LegA equals SRCA's
#   observed mean travel-to-scene (7:46). delta, a uniform per-call residual
#   (operational and model-related), is then solved so the modelled 8-min
#   share equals the observed 45.43%; the report's mean response (10:31) and
#   stroke-pathway response (11 min) serve as consistency checks.
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
    stemi_resp_min   = 12,        # STEMI-pathway response  (p99, Table 35/2)
    stemi_hosp_min   = 37         # STEMI-pathway mean time to a specialised
                                  #   cardiac hospital, 1,680 patients
                                  #   (p99, Table 35/2): consistency check
  ),

  # --- model scenarios and standards ---
  scene_min   = 15,             # on-scene minutes (primary; 10/20 sensitivity)
  # Door-to-device for STEMI is ASSUMED: SRCA's pathway takes patients
  # straight to the cath lab (ED bypass), and no Saudi door-to-device time
  # for ED-bypass ambulance patients is published. 45 min and the STARS-2
  # national median door-to-balloon (63 min; mostly ED arrivals, 8.5% by
  # EMS) are sensitivity rows. SRCA's 37 min (Table 35/2) is the time to
  # REACH a cardiac hospital, not a hospital interval.
  door_device_min  = 30,
  door_device_sens = c(45, 63),
  std = list(
    response  = 8,              # SRCA KPI / NFPA 1710
    trauma    = 60,             # golden-hour convention (Branas 2005;
                                #   causal caveat: Newgard 2010)
    stemi_door  = 90,           # call-to-door (secondary; not a guideline metric)
    stemi_fmc   = 90,           # ACC/AHA: first medical contact (ambulance
                                #   arrival) to device <= 90 min (primary)
    stemi_fmc_strategy = 120,   # FMC-to-device > 120 min favors fibrinolysis
    stroke    = 60              # call-to-door access convention (Adeoye 2014);
                                #   not an AHA/ASA treatment target
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

# Snap only to each road component's strongly connected core: vertices that
# can both reach and be reached from a central hub of their component. A
# nearest vertex on a one-way stub is unreachable, which silently dropped
# Tathlith General Hospital and three SRCA stations before this guard
# (2026-10-08). Hub search: one forward and one reversed Dijkstra per
# component; a hub whose core covers < 90% of its component is replaced.
core_ids <- cache("graph_core_vertices", {
  flip <- c(from_id = "to_id",   from_lon = "to_lon",   from_lat = "to_lat",
            to_id   = "from_id", to_lon   = "from_lon", to_lat   = "from_lat")
  rows <- split(seq_len(nrow(graph)), graph$component)
  unlist(lapply(rows, function(r) {
    gc <- graph[r, ]
    grc <- gc; names(grc)[match(names(flip), names(grc))] <- unname(flip)
    vc <- dodgr::dodgr_vertices(gc)
    ord <- order((vc$x - median(vc$x))^2 + (vc$y - median(vc$y))^2)
    best <- character()
    for (h in vc$id[head(ord, 5)]) {
      fwd <- dodgr::dodgr_times(gc,  from = h, to = vc$id)
      bwd <- dodgr::dodgr_times(grc, from = h, to = vc$id)
      core <- vc$id[is.finite(fwd[1, ]) & is.finite(bwd[1, ])]
      if (length(core) > length(best)) best <- core
      if (length(best) >= 0.9 * nrow(vc)) break
    }
    best
  }), use.names = FALSE)
})
message(sprintf("Routable core: %s of %s vertices",
                format(length(core_ids), big.mark = ","),
                format(nrow(verts), big.mark = ",")))
verts   <- verts[verts$id %in% core_ids, ]
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

# Trauma-capability tiers (sensitivity). No national trauma-centre
# designation list exists; two Riyadh centres align with ACS Level I and
# tertiary hospitals manage major trauma elsewhere (Khan 2025). Tiers:
#   tertiary-capable  : public hospitals SRCA designates for BOTH STEMI and
#                       stroke (same facility: MoH facility ID, otherwise
#                       identical coordinates), plus the Level I-equivalent
#                       centres;
#   Level I-equivalent: King Saud Medical City and King Abdulaziz Medical
#                       City Riyadh.
fac_cols <- c("name_en", "sector", "region_en", "governorate_en", "lat", "lon",
              "moh_facility_id")
pub_cath   <- cath   |> filter(sector != "Private") |> select(all_of(fac_cols))
pub_stroke <- stroke |> filter(sector != "Private") |> select(all_of(fac_cols))
on_both <- with(pub_cath,
  (!is.na(moh_facility_id) &
     moh_facility_id %in% stats::na.omit(pub_stroke$moh_facility_id)) |
  paste(lat, lon) %in% paste(pub_stroke$lat, pub_stroke$lon))
level1 <- bind_rows(cath, stroke) |> select(all_of(fac_cols)) |>
  filter(name_en %in% c("King Saud Medical City",
                        "King Abdulaziz Medical City Riyadh")) |>
  distinct(name_en, .keep_all = TRUE)
stopifnot(nrow(level1) == 2)
tier3 <- bind_rows(pub_cath[on_both, ], level1) |>
  distinct(lat, lon, .keep_all = TRUE)
tier3_sf  <- st_as_sf(tier3,  coords = c("lon", "lat"), crs = 4326) |>
  st_transform(PAR$crs_m)
level1_sf <- st_as_sf(level1, coords = c("lon", "lat"), crs = 4326) |>
  st_transform(PAR$crs_m)
# Broad tertiary tier: the GENEROUS end of the tier range. The tertiary-
# capable tier (both SRCA designations at one facility) is the conservative
# end; the broad tier also counts towns that hold the two designations
# between two hospitals. Members: the tertiary-capable tier; the registry's
# general tertiary hospitals (medical cities and specialist hospitals;
# single-specialty centres excluded); and every public stroke receiving
# hospital whose governorate (boundary polygon) has a public catheterization
# centre in the SAME TOWN, i.e. both hospitals fall in one connected cluster
# of urban or peri-urban cells (>= 300 residents per km2; cells touching at
# edges or corners). A pair is counted at the stroke receiving hospital: the
# hospital with emergency care and CT that an injured patient would reach; a
# cardiac centre alone is not a trauma destination. Points within 0.5 km of
# an earlier member are one campus and dropped. Rule adopted 2026-10-08
# (reviewer round 3): a 5-km cutoff missed Taif's pair at 5.48 km, and any
# cutoff can be nudged; the governorate is the paper's own unit.
reg3 <- prov |>
  filter(type == "3ry", scope %in% c("Medical City", "Specialist")) |>
  transmute(name_en, sector = "MoH", region_en, governorate_en, lat, lon,
            moh_facility_id = facility_id)
to_m <- function(d) st_as_sf(d, coords = c("lon", "lat"), crs = 4326) |>
  st_transform(PAR$crs_m)
gov_of <- function(d) gov$Gov_EN[st_nearest_feature(to_m(d), gov)]
# 8-connected clusters of urban/peri-urban cells on the 30-arc-second lattice
ll_cell <- st_coordinates(st_transform(grid, 4326))
lat_ij  <- cbind(round((ll_cell[, 1] - min(ll_cell[, 1])) * 120),
                 round((ll_cell[, 2] - min(ll_cell[, 2])) * 120))
dense   <- which(grid$cell_total >= 300)
nb_d    <- spdep::dnearneigh(lat_ij[dense, ], 0, 1.5)        # 8 neighbours
cluster <- rep(NA_integer_, nrow(grid))
cluster[dense] <- spdep::n.comp.nb(nb_d)$comp.id
town_of <- function(d)                                       # NA = not in a dense cell
  cluster[RANN::nn2(st_coordinates(grid), st_coordinates(to_m(d)), k = 1)$nn.idx[, 1]]
pc_gov <- gov_of(pub_cath); ps_gov <- gov_of(pub_stroke)
pc_town <- town_of(pub_cath); ps_town <- town_of(pub_stroke)
d_sc <- units::drop_units(st_distance(to_m(pub_stroke), to_m(pub_cath))) / 1000
pair_k <- sapply(seq_len(nrow(pub_stroke)), function(i) {
  j <- which(pc_gov == ps_gov[i])
  if (length(j)) j[which.min(d_sc[i, j])] else NA_integer_
})
has_pair <- !is.na(pair_k)
km_pair  <- ifelse(has_pair, d_sc[cbind(seq_along(pair_k), pair_k)], NA_real_)
same_town <- has_pair & ((km_pair <= 0.5) %in% TRUE |     # NA cluster = not same town
  (ps_town == pc_town[pair_k]) %in% TRUE)
stopifnot(!anyNA(same_town))
broad3 <- bind_rows(tier3, reg3, pub_stroke[same_town, ])
dmat <- units::drop_units(st_distance(to_m(broad3)))
keep <- rep(TRUE, nrow(broad3))
for (i in seq_len(nrow(broad3))[-1])
  keep[i] <- !any(dmat[i, seq_len(i - 1)][keep[seq_len(i - 1)]] <= 500)
added <- which(same_town)[keep[nrow(tier3) + nrow(reg3) + seq_len(sum(same_town))]]
broad3 <- broad3[keep, ]
broad3_sf <- to_m(broad3)
# Every separate-facility pair (stroke hospital vs nearest same-governorate
# cath centre, > 0.5 km apart), with road distance (shortest path,
# stroke -> cath). The coordinate version stays private (data/processed);
# the published table carries names, distances and the town rule only.
pr <- which(has_pair & km_pair > 0.5)
road_km <- vapply(pr, function(i) {
  dd <- dodgr::dodgr_dists(graph, from = snap_verts(to_m(pub_stroke[i, ]), "pair_s"),
                           to = snap_verts(to_m(pub_cath[pair_k[i], ]), "pair_c"),
                           shortest = TRUE)
  dd[1, 1] / 1000
}, numeric(1))
pairs_tab <- tibble(
  governorate      = ps_gov[pr],
  cath_hospital    = pub_cath$name_en[pair_k[pr]],
  cath_lat         = pub_cath$lat[pair_k[pr]], cath_lon = pub_cath$lon[pair_k[pr]],
  stroke_hospital  = pub_stroke$name_en[pr],
  stroke_lat       = pub_stroke$lat[pr], stroke_lon = pub_stroke$lon[pr],
  straight_km      = km_pair[pr], road_km = road_km,
  town_rule        = ifelse(same_town[pr], "same town", "different towns"),
  added_to_broad_tier = pr %in% added) |>
  arrange(straight_km)
write_csv(pairs_tab, file.path(prc, "broad_tier_pairs.csv"))
write_csv(select(pairs_tab, -ends_with("_lat"), -ends_with("_lon")),
          file.path(tabs, "tableS_broad_tier_pairs.csv"))
message(sprintf("Trauma tiers: broad tertiary n=%d | tertiary-capable n=%d | Level I-equivalent n=%d",
                nrow(broad3_sf), nrow(tier3_sf), nrow(level1_sf)))
# routing cache keyed to the member set, so a rule change cannot reuse it
broad_key <- paste(sprintf("%.6f,%.6f", broad3$lat, broad3$lon), collapse = ";")
if (is.null(tt$trauma_broad3) || !identical(attr(tt$trauma_broad3, "key"), broad_key)) {
  tt$trauma_broad3 <- tt_nearest(grid, broad3_sf, label = "broad3")
  attr(tt$trauma_broad3, "key") <- broad_key
  saveRDS(tt, file.path(prc, "travel_times_grid.rds"))
}
if (is.null(tt$trauma_tier3)) {
  tt$trauma_tier3 <- tt_nearest(grid, tier3_sf, label = "tier3")
  saveRDS(tt, file.path(prc, "travel_times_grid.rds"))
}
if (is.null(tt$trauma_level1)) {
  tt$trauma_level1 <- tt_nearest(grid, level1_sf, label = "level1")
  saveRDS(tt, file.path(prc, "travel_times_grid.rds"))
}

acc <- grid |>
  mutate(t_legA   = tt$legA$minutes,   srca_row   = tt$legA$dest_row,
         t_trauma = tt$trauma$minutes, trauma_row = tt$trauma$dest_row,
         t_cath   = tt$cath$minutes,   cath_row   = tt$cath$dest_row,
         t_strk1  = tt$strk1$minutes,  strk1_row  = tt$strk1$dest_row,
         t_strk2  = if (has_t2) tt$strk2$minutes else NA_real_)

# Reachability guard (2026-10-08 defect: Tathlith General Hospital and three
# SRCA stations had snapped to unreachable one-way vertices). A facility that
# is the straight-line nearest for residents yet receives none by road must
# share a campus (<= 0.5 km) with another facility of its set, or be routable
# from (station: to) its most populous straight-line cell: finite road time.
reach_rows <- function(set, fac_sf, rows) {
  eu <- RANN::nn2(st_coordinates(fac_sf), st_coordinates(grid), k = 1)$nn.idx[, 1]
  lv <- seq_len(nrow(fac_sf))
  pe <- tapply(grid$cell_total, factor(eu, levels = lv), sum)
  pr <- tapply(grid$cell_total, factor(rows, levels = lv), sum)
  pe[is.na(pe)] <- 0; pr[is.na(pr)] <- 0
  d <- units::drop_units(st_distance(fac_sf)); diag(d) <- Inf
  near <- apply(d, 1, min) <= 500
  routable <- rep(NA, length(lv))
  for (i in which(pe > 0 & pr == 0 & !near)) {
    cell <- which(eu == i)[which.max(grid$cell_total[eu == i])]
    fv <- snap_verts(fac_sf[i, ], set); cv <- snap_verts(grid[cell, ], "cell")
    tm <- if (set == "srca") dodgr::dodgr_times(graph, from = fv, to = cv)
          else dodgr::dodgr_times(graph, from = cv, to = fv)
    routable[i] <- is.finite(tm[1, 1])
  }
  nmcol <- intersect(c("name_en", "center_en", "name"), names(fac_sf))[1]
  tibble(set = set, name = fac_sf[[nmcol]], euclid_pop = as.numeric(pe),
         routed_pop = as.numeric(pr), co_located = near, routable = routable,
         ok = !(as.numeric(pe) > 0 & as.numeric(pr) == 0) | near | routable %in% TRUE)
}
reach <- bind_rows(
  reach_rows("srca",   srca_sf,   tt$legA$dest_row),
  reach_rows("trauma", trauma_sf, tt$trauma$dest_row),
  reach_rows("cath",   cath_sf,   tt$cath$dest_row),
  reach_rows("stroke", stroke_sf, tt$strk1$dest_row),
  reach_rows("tier3",  tier3_sf,  tt$trauma_tier3$dest_row),
  reach_rows("broad3", broad3_sf, tt$trauma_broad3$dest_row))
write_csv(reach, file.path(prc, "reachability_check.csv"))
if (any(!reach$ok))
  warning(sprintf("%d facilities unreachable: %s", sum(!reach$ok),
                  paste(reach$name[!reach$ok], collapse = "; ")))

# Population-weighted summary helpers
wq <- function(x, w, p) {
  ok <- !is.na(x) & !is.na(w); x <- x[ok]; w <- w[ok]
  if (!length(x)) return(NA_real_)
  o <- order(x); x <- x[o]; w <- w[o]
  x[which.max(cumsum(w) / sum(w) >= p)]
}
cov_pct <- function(x, w, thr) sum(w[!is.na(x) & x <= thr]) / sum(w) * 100

# ---- 5. Calibration ------------------------------------------------------
# One moment fitted: alpha scales unadjusted car travel times so the
# DEMAND-weighted mean of alpha*LegA equals the observed mean travel-to-
# scene (7:46). Demand weight = cell population (call volume proxy: SRCA
# means are call-weighted and therefore urban-skewed; an unweighted mean
# over cells would compare a population-weighted model to a call-weighted
# statistic). The residual delta and the consistency checks follow in Section 7.

alpha <- PAR$obs$travel_scene_min /
  weighted.mean(acc$t_legA, acc$cell_total, na.rm = TRUE)
message(sprintf(
  "Calibration: alpha = %.3f (raw demand-weighted LegA mean %.2f min -> %.2f)",
  alpha, weighted.mean(acc$t_legA, acc$cell_total, na.rm = TRUE),
  PAR$obs$travel_scene_min))

acc <- acc |> mutate(t_resp = FIXED_PRE + alpha * t_legA)

# The model-vs-report comparison (validation_moments.csv) is written in
# Section 7, once the residual delta is solved: one fitted parameter (alpha),
# one solved residual (delta), and consistency checks. The base-model mean
# response (FIXED_PRE + fitted 7:46 = 9:01) equals the sum of the report's
# own intervals by construction, so it is not an independent test.

# ---- 6. Full-chain call-to-door times and standards compliance -----------

acc <- acc |>
  mutate(
    t_door_trauma = t_resp + PAR$scene_min + alpha * t_trauma,
    t_door_cath   = t_resp + PAR$scene_min + alpha * t_cath,
    t_door_strk1  = t_resp + PAR$scene_min + alpha * t_strk1,
    t_door_strk2  = t_resp + PAR$scene_min + alpha * t_strk2,
    # STEMI guideline metric: first medical contact (ambulance arrival at the
    # patient) to device = on scene + transport + assumed door-to-device.
    # Excludes the response leg, so the residual delay does not touch it.
    t_fmc_device  = PAR$scene_min + alpha * t_cath + PAR$door_device_min
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
  natl_row(acc$t_fmc_device,  w, PAR$std$stemi_fmc, "STEMI: first medical contact to device"),
  natl_row(acc$t_door_cath,   w, PAR$std$stemi_door,"STEMI: call to cath-centre door"),
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
    stemi_med      = wq(t_fmc_device, cell_total, .5),
    stemi_std_pct  = cov_pct(t_fmc_device, cell_total, PAR$std$stemi_fmc),
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
    stemi_std_pct  = cov_pct(t_fmc_device, cell_total, PAR$std$stemi_fmc),
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

# Urban-rural gradient. Density thresholds ADAPTED from the Degree of
# Urbanisation (EC/GHSL DEGURBA) on 1-km cells: >= 1,500 residents/km2 urban;
# 300-1,500 peri-urban; < 300 rural. Cell density only: DEGURBA's contiguity
# and minimum-cluster-size rules are not applied. Ties the access deficit to GASTAT's finding
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
    stemi_std_pct  = cov_pct(t_fmc_device, cell_total, PAR$std$stemi_fmc),
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
# Residual-delay scenario: the modelled KPI exceeds the observed one
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
               "Residual per-call delay delta (min)",
               "Modelled KPI with residual delay reduced by 0.5 min",
               "Modelled KPI with residual delay reduced by 1 min"),
  value = c(cov_pct(alpha * acc$t_legA, w, PAR$std$response),
            cov_pct(acc$t_resp, w, PAR$std$response),
            PAR$obs$resp8_share,
            delta_op,
            cov_pct(acc$t_resp_op - 0.5, w, PAR$std$response),
            cov_pct(acc$t_resp_op - 1, w, PAR$std$response)),
  unit = c("percent", "percent", "percent", "minutes", "percent", "percent"))
# The report's own intervals (0:48 + 0:27 + 7:46 = 9:01) fall short of its
# reported mean response (10:31): a model-free residual of about 1.5 min,
# close to delta. Recorded beside the KPI decomposition.
kpi <- bind_rows(kpi, tibble(
  quantity = "SRCA report: mean response minus sum of reported intervals (min)",
  value = PAR$obs$resp_mean_min - FIXED_PRE - PAR$obs$travel_scene_min,
  unit = "minutes"))
write_csv(kpi, file.path(tabs, "tableS_kpi_decomposition.csv"))
print(kpi)

resp_mean_base <- weighted.mean(acc$t_resp, w, na.rm = TRUE)
validation <- tibble(
  moment   = c("Mean travel to scene (min)",
               "Responses within 8 min (%)",
               "Mean response (min)",
               "Stroke-pathway mean response (min)",
               "STEMI-pathway mean response (min)",
               "STEMI-pathway time to cardiac hospital (min)",
               "STEMI-pathway time to cardiac hospital (min)"),
  role     = c("fitted (alpha)",
               "base model; the gap defines delta",
               "consistency check (base model + delta)",
               "consistency check (base model + delta)",
               "consistency check (base model + delta)",
               "consistency check (model median from call, base + delta)",
               "consistency check (model median from first medical contact)"),
  model_label = c("Base model", "Base model", "Model + \u03b4", "Model + \u03b4",
                  "Model + \u03b4", "Median from call + \u03b4",
                  "Median from contact"),
  observed = c(PAR$obs$travel_scene_min, PAR$obs$resp8_share,
               PAR$obs$resp_mean_min, PAR$obs$stroke_resp_min,
               PAR$obs$stemi_resp_min, PAR$obs$stemi_hosp_min,
               PAR$obs$stemi_hosp_min),
  modelled = c(weighted.mean(alpha * acc$t_legA, w, na.rm = TRUE),
               cov_pct(acc$t_resp, w, PAR$std$response),
               resp_mean_base + delta_op,
               resp_mean_base + delta_op,
               resp_mean_base + delta_op,
               wq(acc$t_door_cath, w, .5) + delta_op,
               wq(PAR$scene_min + alpha * acc$t_cath, w, .5)))
write_csv(validation, file.path(tabs, "validation_moments.csv"))
print(validation)

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
# Theil-T with a nested three-level decomposition: between regions, between
# governorates within regions, within governorates (0 * log 0 = 0, so cells
# with zero driving time keep their weight in the mean).
theil_nested <- function(x, w, reg, gov) {
  ok <- !is.na(x) & !is.na(w) & w > 0 & x >= 0 & !is.na(reg) & !is.na(gov)
  x <- x[ok]; w <- w[ok]; reg <- reg[ok]; gov <- paste(reg[ok], gov[ok])
  tl  <- function(r) ifelse(r > 0, r * log(r), 0)
  s   <- w / sum(w); mu <- sum(s * x)
  mean_by  <- function(g) tapply(x * w, g, sum) / tapply(w, g, sum)
  share_by <- function(g) tapply(w, g, sum) / sum(w)
  total   <- sum(s * tl(x / mu))
  b_reg   <- sum(share_by(reg) * tl(mean_by(reg) / mu))
  b_gov   <- sum(share_by(gov) * tl(mean_by(gov) / mu))   # all governorates
  tibble(total = total, between_region = b_reg,
         between_gov_within_region = b_gov - b_reg, within_gov = total - b_gov)
}

# Inequality is measured on the DRIVING components only. Relative measures
# fall when a constant is added, and the outcomes carry different constants
# (1.25 min for response; 16.25 min for door times; 45 min for STEMI
# first-medical-contact-to-device), so totals would bias the comparison.
drive <- list(
  "Response"       = alpha * acc$t_legA,
  "Trauma door"    = alpha * (acc$t_legA + acc$t_trauma),
  "STEMI FMC-to-device" = alpha * acc$t_cath,
  "Stroke door T1" = alpha * (acc$t_legA + acc$t_strk1))
eq_tab <- bind_rows(lapply(drive, function(x) {
  td <- theil_nested(x, w, acc$region_en, acc$Gov_EN)
  tibble(gini = gini_w(x, w), theil_total = td$total,
         between_region_pct = td$between_region / td$total * 100,
         between_gov_within_region_pct = td$between_gov_within_region / td$total * 100,
         within_gov_pct = td$within_gov / td$total * 100,
         between_gov_pct = (td$between_region + td$between_gov_within_region) /
           td$total * 100)
}), .id = "metric")
write_csv(eq_tab, file.path(tabs, "tableS_inequality.csv"))
print(eq_tab, width = Inf)

# ---- 9b. Spatial clustering of access deficits (descriptive, supplement) --
# Is governorate-level compliance spatially clustered? Reported descriptively
# in the supplement: neighboring governorates share nearest facilities, so
# positive autocorrelation is expected by construction. Permutation inference is the only
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
                stemi_std_pct  = "STEMI FMC-to-device <= 90 min",
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
                     extra = 0, b = a, dd = PAR$door_device_min) {
  resp <- FIXED_PRE + extra + a * acc$t_legA
  fmc  <- scene + b * t_ca + dd
  tibble(scenario = scenario,
         trauma_std_pct = cov_pct(resp + scene + b * t_tr, wgt, PAR$std$trauma),
         stemi_std_pct  = cov_pct(fmc, wgt, PAR$std$stemi_fmc),
         stemi_120_pct  = cov_pct(fmc, wgt, PAR$std$stemi_fmc_strategy),
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
  sens_row(sprintf("Residual-delay scenario (+%.1f min, KPI-matched)",
                   delta_op), PAR$scene_min, alpha, w,
           acc$t_trauma, acc$t_cath, acc$t_strk1, extra = delta_op),
  sens_row("Leg B uncalibrated (scene-to-door speed factor 1.0)",
           PAR$scene_min, alpha, w,
           acc$t_trauma, acc$t_cath, acc$t_strk1, b = 1),
  sens_row("STEMI: door-to-device 45 min",
           PAR$scene_min, alpha, w, acc$t_trauma, acc$t_cath, acc$t_strk1,
           dd = PAR$door_device_sens[1]),
  sens_row("STEMI: door-to-device 63 min (STARS-2 median)",
           PAR$scene_min, alpha, w, acc$t_trauma, acc$t_cath, acc$t_strk1,
           dd = PAR$door_device_sens[2]),
  sens_row(sprintf("Trauma: >=100-bed hospitals only (n=%d)",
                   nrow(trauma_b100_sf)), PAR$scene_min, alpha, w,
           tt$trauma_b100$minutes, acc$t_cath, acc$t_strk1),
  sens_row(sprintf("Trauma: broad tertiary tier (n=%d)",
                   nrow(broad3_sf)), PAR$scene_min, alpha, w,
           tt$trauma_broad3$minutes, acc$t_cath, acc$t_strk1),
  sens_row(sprintf("Trauma: tertiary-capable hospitals only (n=%d)",
                   nrow(tier3_sf)), PAR$scene_min, alpha, w,
           tt$trauma_tier3$minutes, acc$t_cath, acc$t_strk1),
  sens_row(sprintf("Trauma: Level I-equivalent centers only (n=%d)",
                   nrow(level1_sf)), PAR$scene_min, alpha, w,
           tt$trauma_level1$minutes, acc$t_cath, acc$t_strk1),
  if (has_t2)
    sens_row("Stroke Tier 2 destinations", PAR$scene_min, alpha, w,
             acc$t_trauma, acc$t_cath, acc$t_strk2)
)
write_csv(sens, file.path(tabs, "tableS_sensitivity.csv"))
print(sens)

# Trauma access under the three destination definitions, by region
acc <- acc |> mutate(
  t_door_broad3 = t_resp + PAR$scene_min + alpha * tt$trauma_broad3$minutes,
  t_door_tier3  = t_resp + PAR$scene_min + alpha * tt$trauma_tier3$minutes,
  t_door_level1 = t_resp + PAR$scene_min + alpha * tt$trauma_level1$minutes)
tier_sum <- function(d) d |> summarise(
  pop = sum(cell_total),
  general_hospital_pct  = cov_pct(t_door_trauma, cell_total, PAR$std$trauma),
  broad_tertiary_pct    = cov_pct(t_door_broad3, cell_total, PAR$std$trauma),
  tertiary_capable_pct  = cov_pct(t_door_tier3,  cell_total, PAR$std$trauma),
  level1_equivalent_pct = cov_pct(t_door_level1, cell_total, PAR$std$trauma),
  .groups = "drop")
acc_df <- st_drop_geometry(acc)
tier_tab <- bind_rows(
  acc_df |> group_by(region_en) |> tier_sum() |> arrange(desc(pop)),
  tier_sum(acc_df) |> mutate(region_en = "National", .before = 1))
write_csv(tier_tab, file.path(tabs, "tableS_trauma_tiers.csv"))
print(tier_tab, width = Inf)

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
    t_fmc_device  = weighted.mean(t_fmc_device,  cell_total, na.rm = TRUE),
    t_door_strk1  = weighted.mean(t_door_strk1,  cell_total, na.rm = TRUE),
    pop           = sum(cell_total),
    .groups = "drop") |>
  mutate(x = bx * 1e4 + 5e3, y = by * 1e4 + 5e3,
         # opacity by block population (log scale): sparsely populated desert
         # blocks recede so the map reads like the population-weighted shares
         pop_alpha = c(0.18, 0.55, 1)[findInterval(pop, c(1e3, 1e4)) + 1])

cities <- tibble::tribble(
  ~name,     ~lon,  ~lat,  ~hj,  ~vj,
  "Riyadh",  46.71, 24.63, 0.5, -0.7,
  "Jeddah",  39.19, 21.49, 1.1,  0.5,    # left of the point (Makkah is 60 km east)
  "Makkah",  39.83, 21.42, -0.1, 0.5,    # right of the point
  "Madinah", 39.61, 24.47, 0.5, -0.7,
  "Dammam",  50.10, 26.43, 0.5, -0.7,
  "Abha",    42.51, 18.22, 0.5, -0.7,
  "Tabuk",   36.58, 28.38, 0.5, -0.7) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_transform(crs_disp)
cities <- cities |>
  mutate(x = st_coordinates(cities)[, 1], y = st_coordinates(cities)[, 2])

panel_map <- function(tcol, fac_layer, title) {
  d <- blk |> mutate(bin = time_bins(.data[[tcol]])) |> filter(!is.na(bin))
  ggplot() +
    geom_sf(data = gov_disp, fill = "grey92", color = "grey78",
            linewidth = 0.12) +
    # show.legend = TRUE keeps keys for unused bins (panel B has no cell
    # <= 30 min), so the three panel legends stay identical and collect
    geom_tile(data = d, aes(x, y, fill = bin, alpha = pop_alpha),
              width = 1e4, height = 1e4, show.legend = TRUE) +
    scale_alpha_identity() +
    geom_sf(data = st_transform(fac_layer, crs_disp), shape = 3, size = 0.7,
            color = "black", stroke = 0.5) +
    geom_text(data = st_drop_geometry(cities),
              aes(x, y, label = name, hjust = hj, vjust = vj),
              size = 2.3, color = "grey10", fontface = "bold") +
    scale_fill_manual(values = bin_cols, name = "Minutes",
                      drop = FALSE) +
    labs(title = title) + theme_map
}

fig1 <- panel_map("t_door_trauma", trauma_sf,
                  sprintf("A  Trauma, call to door\n(n=%d general hospitals)", nrow(trauma_sf))) +
        panel_map("t_fmc_device", cath_sf,
                  sprintf("B  STEMI, first medical contact to device\n(n=%d catheterization centers)", nrow(cath_sf))) +
        panel_map("t_door_strk1", stroke_sf,
                  sprintf("C  Stroke, call to door\n(n=%d receiving hospitals)", nrow(stroke_sf))) +
        patchwork::plot_layout(ncol = 3, guides = "collect") &
        theme(legend.position = "bottom")
ggsave(file.path(figs, "fig1_maps_call_to_door.png"), fig1,
       width = 13, height = 5.6, dpi = 300, bg = "white")

# Fig 2 -- calibration/validation + KPI decomposition. Panel A: one facet
# per moment (own y scale, so minutes and percent never share an axis).
# Panel A omits the fitted moment (travel to scene, equal by construction).
lab_lv <- c("SRCA", "Base model", "Model + \u03b4", "Median from call + \u03b4",
            "Median from contact")
f2a <- validation |>
  filter(role != "fitted (alpha)") |>
  mutate(moment = factor(str_wrap(moment, 18), unique(str_wrap(moment, 18)))) |>
  (\(v) bind_rows(
    distinct(v, moment, observed) |> transmute(moment, lab = "SRCA", value = observed),
    transmute(v, moment, lab = model_label, value = modelled)))() |>
  mutate(lab = factor(lab, lab_lv)) |>
  ggplot(aes(lab, value, fill = lab)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f", value)), vjust = -0.35, size = 2.3) +
  facet_wrap(~moment, scales = "free", nrow = 1) +
  scale_fill_manual(values = c("SRCA" = "grey35", "Base model" = "steelblue",
                               "Model + \u03b4" = "#7FA7CF",
                               "Median from call + \u03b4" = "#7FA7CF",
                               "Median from contact" = "#B8CFE6")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "A  Model against SRCA 2025 moments", x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(strip.text = element_text(size = 6.5),
        axis.text.x = element_text(size = 6.5, angle = 30, hjust = 1))

resp_grid <- tibble(t = seq(0, 30, 0.25)) |>
  rowwise() |>
  mutate(model  = cov_pct(acc$t_resp, w, t),
         travel = cov_pct(alpha * acc$t_legA, w, t),
         oper   = cov_pct(acc$t_resp_op, w, t)) |>
  ungroup()
lt <- c("Travel only (structural ceiling)"   = 2,
        "With queue + mobilization"          = 1,
        "With residual delay (KPI-matched)" = 3)
f2b <- ggplot(resp_grid, aes(t)) +
  geom_line(aes(y = travel, linetype = "Travel only (structural ceiling)")) +
  geom_line(aes(y = model,  linetype = "With queue + mobilization")) +
  geom_line(aes(y = oper,
                linetype = "With residual delay (KPI-matched)")) +
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
ggsave(file.path(figs, "fig2_calibration_kpi.png"),
       f2a + f2b + patchwork::plot_layout(widths = c(1.7, 1)),
       width = 12.5, height = 4.6, dpi = 300, bg = "white")

# Fig 3 -- where inequality in driving time lies: nested Theil shares
# (between regions / between governorates within regions / within
# governorates), driving components only; Gini in the axis label. No
# in-figure title (the legend carries it).
th_long <- eq_tab |>
  transmute(metric = sprintf("%s\n(Gini %.2f)", c(
              "Response" = "Response: drive to scene",
              "Trauma door" = "Trauma: drive to scene + to door",
              "STEMI FMC-to-device" = "STEMI: drive scene to door",
              "Stroke door T1" = "Stroke: drive to scene + to door")[metric], gini),
            `Between regions` = between_region_pct,
            `Between governorates within regions` = between_gov_within_region_pct,
            `Within governorates` = within_gov_pct) |>
  pivot_longer(-metric, names_to = "level", values_to = "pct") |>
  mutate(level = factor(level, c("Within governorates",
                                 "Between governorates within regions",
                                 "Between regions")),
         metric = factor(metric, rev(unique(metric))))
fig3 <- ggplot(th_long, aes(pct, metric, fill = level)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = ifelse(pct >= 4, sprintf("%.0f%%", pct), "")),
            position = position_stack(vjust = 0.5), size = 2.6,
            color = "white", fontface = "bold") +
  scale_fill_manual(values = c("Within governorates" = "grey55",
                               "Between governorates within regions" = "#2C7BB6",
                               "Between regions" = "#D7191C"),
                    breaks = c("Between regions",
                               "Between governorates within regions",
                               "Within governorates"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.03)),
                     labels = function(v) paste0(v, "%")) +
  labs(x = "share of Theil T inequality in driving time (ambulance drive legs only)",
       y = NULL) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank()) +
  guides(fill = guide_legend(nrow = 2))
ggsave(file.path(figs, "fig3_theil.png"), fig3,
       width = 6.2, height = 3.6, dpi = 300, bg = "white")

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
fig4 <- f4p("trauma_med", "A  Trauma") + f4p("stemi_med", "B  STEMI (contact to device)") +
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
         t_legA, t_resp, t_resp_op, starts_with("t_door"), t_fmc_device,
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
#       the STEMI first-medical-contact-to-device window (90 min; ACC/AHA;
#       30-min door-to-device assumed) or the 60-min stroke call-to-door
#       window (base case).
#       STEMI candidates: the >=100-bed general/tertiary subset (cath-lab
#       feasibility), minus hospitals already on the SRCA cath list
#       (matched exactly by moh_facility_id). Stroke candidates: all
#       general hospitals minus those already on the SRCA stroke list.
# Routing for (b) is direction-true (cell -> hospital): the street graph is
# reversed by renaming its from_/to_ columns, so one Dijkstra per CANDIDATE
# (n=250) replaces one per cell (n=67,075). Cached sparse at <= 80 raw
# minutes; the widest per-cell budget is STEMI's (90 - 15 - 30)/(0.9 alpha) = 53.

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
# Governorate from the boundary polygons, not the registry label (the
# registry files Al Majardah General Hospital under Bariq; it lies in Al
# Majaridah). Nearest polygon for points on a coastline or border.
cand_gov <- gov$Gov_EN[st_nearest_feature(cand_sf, st_transform(gov, PAR$crs_m))]
message(sprintf("Upgrade candidates whose registry governorate differs from the polygon: %d",
                sum(cand_gov != cand_sf$governorate_en, na.rm = TRUE)))

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

# fixed = the outcome's minutes outside the transport leg: on scene + assumed
# door-to-device for STEMI first-medical-contact-to-device; response + on
# scene for stroke and trauma call-to-door.
greedy_upgrade <- function(idx, door_base, window, km_existing, fixed,
                           top_n = 15, min_gain = 10000) {
  budget  <- (window - fixed) / alpha                       # raw-minute budget
  budget  <- rep_len(budget, nrow(acc))
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
      governorate = cand_gov[best],
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

fixed_stemi <- PAR$scene_min + PAR$door_device_min
up_stemi  <- greedy_upgrade(idx_stemi,  acc$t_fmc_device, PAR$std$stemi_fmc, km_cath,
                            fixed = fixed_stemi)
up_stroke <- greedy_upgrade(idx_stroke, acc$t_door_strk1, PAR$std$stroke,    km_stroke,
                            fixed = acc$t_resp + PAR$scene_min)
# Trauma-centre upgrades, ranked from the BROAD baseline (reviewer round 3):
# only MoH secondary general hospitals (>= 100 beds) that are tertiary under
# no definition used here are candidates (not a broad-tier member by MoH ID
# or within 0.5 km of one); each round adds the hospital whose upgrade newly
# covers the most residents within 60 minutes of a broad-tier hospital.
broad_moh <- as.character(stats::na.omit(broad3$moh_facility_id))
km_broad3 <- km_to_nearest(broad3_sf)
idx_trauma <- which(cand_sf$type == "2ry" &
                    replace_na(cand_sf$bed_capacity >= 100, FALSE) &
                    !(cand_id %in% broad_moh) & km_broad3 > 0.5)
up_trauma <- greedy_upgrade(idx_trauma, acc$t_door_broad3, PAR$std$trauma, km_broad3,
                            fixed = acc$t_resp + PAR$scene_min)
message(sprintf("Trauma-centre upgrade candidates: %d (2ry, >=100 beds, outside the broad tier)",
                length(idx_trauma)))
write_csv(up_stemi,  file.path(tabs, "tableS_upgrade_sites_stemi.csv"))
write_csv(up_stroke, file.path(tabs, "tableS_upgrade_sites_stroke.csv"))
write_csv(up_trauma, file.path(tabs, "tableS_upgrade_sites_trauma.csv"))
print(up_trauma, n = 20, width = Inf)
print(up_stemi,  n = 20, width = Inf)
print(up_stroke, n = 20, width = Inf)

# (c) rank stability of the greedy shortlists under the +/-10% speed
# scenarios of Section 10 (alpha x 0.9 / x 1.1): response and transport legs
# both rescale; the candidate time table is in raw minutes and unchanged.
# Reports how much of each shortlist, and of its order, survives.
greedy_picks <- function(idx, door, window, a, fixed, top_n = 15,
                         min_gain = 10000) {
  budget  <- rep_len((window - fixed) / a, nrow(acc))
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
                   Stroke = up_stroke$facility_id[-1],
                   Trauma = up_trauma$facility_id[-1])
stab_row <- function(cn, sc, run) {
  b <- base_lists[[cn]]; p <- run$picks
  tibble(condition = cn, scenario = sc,
         start_pct = run$start, cum_within_pct = run$cum,
         gain_pts = run$cum - run$start,
         top5_retained  = length(intersect(head(b, 5), head(p, 5))),
         top10_retained = length(intersect(head(b, 10), head(p, 10))),
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
      fixed_stemi + a * acc$t_cath, PAR$std$stemi_fmc, a, fixed_stemi)),
    stab_row("Stroke", sc, greedy_picks(idx_stroke,
      t_resp_a + PAR$scene_min + a * acc$t_strk1, PAR$std$stroke, a,
      t_resp_a + PAR$scene_min)),
    stab_row("Trauma", sc, greedy_picks(idx_trauma,
      t_resp_a + PAR$scene_min + a * tt$trauma_broad3$minutes, PAR$std$trauma, a,
      t_resp_a + PAR$scene_min)))
}))
stopifnot(all(stab$top15_retained[stab$scenario == "Base"] ==
              lengths(base_lists)[stab$condition[stab$scenario == "Base"]]))
write_csv(stab, file.path(tabs, "tableS_upgrade_stability.csv"))
print(stab, width = Inf)

# (d) New-station siting for the 8-minute response: the structural lever,
# compared like with like against cutting the per-call residual delta. Same
# greedy maximal-covering heuristic; candidate sites are the 67,075 populated
# 1-km cells (no facility type implied). A cell is covered when
# FIXED_PRE + delta + alpha * drive <= 8, i.e. the operational scenario, so
# gains are measured from the observed 45.4% just as the delta lever is.
# Exact local routing: a drive of thr_raw minutes cannot exceed
# v_max * thr_raw km, so each 0.2-degree tile is routed on the subgraph of
# its 3 x 3 tile block (>= 19 km margin, above that bound).
thr_raw <- (PAR$std$response - FIXED_PRE - delta_op) / alpha   # raw minutes
v_max   <- max(graph$d / pmax(graph$time, 1e-9)) * 3.6         # km/h
stopifnot(v_max * thr_raw / 60 < 19)
covered_op <- !is.na(acc$t_resp_op) & acc$t_resp_op <= PAR$std$response

site_cov <- cache("station_site_coverage", {
  dodgr::dodgr_cache_off()      # thousands of one-off subgraphs; no reuse
  cell_v  <- snap_verts(grid, "cells")
  vxy     <- verts[match(cell_v, verts$id), c("x", "y")]
  tkey    <- function(x, y) paste(floor(x / 0.2), floor(y / 0.2))
  cell_tile <- tkey(vxy$x, vxy$y)
  e_tile  <- tkey(graph$from_lon, graph$from_lat)
  e_rows  <- split(seq_len(nrow(graph)), e_tile)
  c_rows  <- split(seq_len(nrow(grid)), cell_tile)
  out <- list()
  for (tk in names(c_rows)) {
    xy  <- as.integer(strsplit(tk, " ")[[1]])
    blk <- as.vector(outer(xy[1] + (-1:1), xy[2] + (-1:1), paste))
    tgt <- unlist(c_rows[intersect(blk, names(c_rows))], use.names = FALSE)
    tgt <- tgt[!covered_op[tgt]]                    # only uncovered cells gain
    if (!length(tgt)) next
    sub <- graph[unlist(e_rows[intersect(blk, names(e_rows))], use.names = FALSE), ]
    sv  <- dodgr::dodgr_vertices(sub)$id
    src <- c_rows[[tk]]
    src <- src[cell_v[src] %in% sv]; tgt <- tgt[cell_v[tgt] %in% sv]
    if (!length(src) || !length(tgt)) next
    tm <- dodgr::dodgr_times(sub, from = unique(cell_v[src]),
                             to = unique(cell_v[tgt])) / 60
    hit <- which(is.finite(tm) & tm <= thr_raw, arr.ind = TRUE)
    if (!nrow(hit)) next
    fv <- rownames(tm)[hit[, 1]]; tv <- colnames(tm)[hit[, 2]]
    pairs <- merge(data.frame(fv, tv),
                   data.frame(fv = cell_v[src], site = src))
    pairs <- merge(pairs, data.frame(tv = cell_v[tgt], cell = tgt))
    out[[tk]] <- unique(pairs[, c("site", "cell")])
  }
  res <- do.call(rbind, out)
  attr(res, "thr_raw") <- thr_raw
  res
})
stopifnot(isTRUE(all.equal(attr(site_cov, "thr_raw"), thr_raw)))

greedy_sites <- function(top_n = 150, min_gain = 1000) {
  covered <- covered_op
  el  <- site_cov[!covered[site_cov$cell], ]
  cum <- sum(w[covered]) / sum(w) * 100
  rows <- list()
  for (r in seq_len(top_n)) {
    el <- el[!covered[el$cell], , drop = FALSE]
    if (!nrow(el)) break
    gain <- tapply(w[el$cell], el$site, sum)
    if (max(gain) < min_gain) break
    best  <- as.integer(names(gain)[which.max(gain)])
    covered[unique(el$cell[el$site == best])] <- TRUE
    cum <- cum + max(gain) / sum(w) * 100
    rows[[r]] <- tibble(rank = r, cell_id = grid$cell_id[best],
                        governorate = acc$Gov_EN[best], region = acc$region_en[best],
                        dens_class = as.character(acc$dens_class[best]),
                        newly_covered_pop = max(gain), cum_within_8_pct = cum)
  }
  bind_rows(rows) |> mutate(cum_gain_m = cumsum(newly_covered_pop) / 1e6)
}
sites <- greedy_sites()
gain_d1  <- (kpi$value[kpi$quantity == "Modelled KPI with residual delay reduced by 1 min"] -
             PAR$obs$resp8_share) / 100 * sum(w) / 1e6
gain_d05 <- (kpi$value[kpi$quantity == "Modelled KPI with residual delay reduced by 0.5 min"] -
             PAR$obs$resp8_share) / 100 * sum(w) / 1e6
sites <- sites |> mutate(
  matches_delta_0_5min = cum_gain_m >= gain_d05,
  matches_delta_1min   = cum_gain_m >= gain_d1)
write_csv(sites, file.path(tabs, "tableS_station_sites.csv"))
message(sprintf(paste0("New stations (8-min, operational scenario): 15 sites +%.2fM; ",
  "sites to match delta -0.5 min (+%.2fM): %s; delta -1 min (+%.2fM): %s"),
  sites$cum_gain_m[min(15, nrow(sites))], gain_d05,
  if (any(sites$matches_delta_0_5min)) which(sites$matches_delta_0_5min)[1] else ">150",
  gain_d1,
  if (any(sites$matches_delta_1min)) which(sites$matches_delta_1min)[1] else ">150"))

message("Done. Tables -> output/tables | Figures -> output/figures | ",
        sprintf("alpha = %.3f | trauma set '%s' (n=%d)%s",
                alpha, PAR$trauma_set, nrow(trauma_sf),
                if (has_t2) "" else " | stroke Tier 2 pending clinician list"))
