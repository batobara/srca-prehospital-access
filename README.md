# Access to Trauma, Myocardial Infarction, and Stroke Care Through Saudi Arabia's Ambulance Network: A Modeling Study

Replication package for the manuscript (Batobara A, Alfaraidhy M, Lim C; under review). One reproducible R script models the full prehospital chain, call receipt to hospital door, for every populated square kilometer of Saudi Arabia, calibrated to published SRCA 2025 operational moments, and ranks capability upgrades across the MoH (Health Holding Company) hospital network and new ambulance-station sites. Please cite the manuscript when using this code or the tables.

## Contents (three files)

| File | What it is |
|---|---|
| `analysis.R` | The single analysis pipeline (Sections 0-12: supply, demand, routing, calibration, outcomes, KPI decomposition, bypass, inequality, spatial clustering, sensitivity and trauma capability tiers, figures, capability upgrades, rank stability, and new-station siting). |
| `derived_tables.zip` | The manuscript's 19 derived tables exactly as produced by the deposited run; every number in the paper traces to these CSVs (unzip to `output/tables/` to compare against a rerun). |
| `README.md` | This file, including the license. |

## Data files

No input data ship with this release. The facility coordinate files (509 SRCA ambulance stations; the 73 catheterization centers and 71 stroke receiving hospitals on the SRCA lists; the MoH provider-registry compilation) were geocoded facility by facility by the authors, each point verified against satellite imagery and commercial map services, and carry author-added attributes. They are deposited in a separate restricted-access Zenodo record; access requires the authors' permission. Their sources are public and cited in the manuscript (SRCA Annual Report 2025 facility listings; the MoH facility registry), so equivalent files can also be compiled independently using the column schema below.

## Inputs a full rerun needs (place in `raw/`)

| File expected by the script | Source |
|---|---|
| `raw/gcc-states-latest.osm.pbf` | Geofabrik GCC extract, https://download.geofabrik.de/asia/gcc-states.html (~250 MB; OpenStreetMap contributors, ODbL) |
| `raw/sau_ppp_2020_constrained.tif` | WorldPop 2020 constrained 100-m population, doi:10.5258/SOTON/WP00685 |
| `raw/PopulationbyNationalitybyRegionGovernorateCityandNationalityARCSV.csv` | GASTAT Saudi Census 2022, population by nationality/region/governorate/city (stats.gov.sa) |
| `raw/governorate/Governorate.gpkg` | Governorate boundaries (150 units; official boundary release) |
| `raw/Regions/Regions.shp` | The 13 administrative regions (official boundary release) |
| `raw/srca_report_regional_2025.csv` | Transcribe from SRCA Annual Report 2025 Tables 22/2 (missions), 23/2 (transports), 28/2 (launch points): columns `region_en, missions_2025, transports_2025, launch_points_2025`; the script asserts the printed totals (1,416,037; 568,827; 520). |
| `raw/SRCA Centers/srca_centers.csv` | One row per distinct station: `lat, lon` (WGS84). |
| `raw/Cath Centers/cath_centers.csv` | One row per center: `lat, lon, moh_facility_id` (MoH registry ID, blank if not MoH), `sector` (MoH, Private, University, Military, National Guard, Other government). |
| `raw/Stroke Centers/stroke_centers.csv` | As for catheterization centers; an optional `reperfusion_capable` column (TRUE/FALSE) activates the Tier-2 stroke arm. Both lists also need `name_en`: the trauma capability tiers match the two Level I-equivalent centers by name (`King Saud Medical City`, `King Abdulaziz Medical City Riyadh`). |
| `raw/list_of_healthcare_providers.csv` | MoH secondary and tertiary hospitals: `facility_id, name_en, type` (2ry or 3ry), `scope, bed_capacity, lat, lon, governorate_en, region_en, region_ar, cluster_en`. The broad trauma tier takes type-3ry rows whose `scope` is `Medical City` or `Specialist`, plus public stroke receiving hospitals in the same governorate and the same town (one connected cluster of urban or peri-urban cells) as a public catheterization center. |

## Running

```bash
Rscript analysis.R
```

Dependencies install automatically from CRAN on first run (sf, terra, dplyr, tidyr, readr, stringr, dodgr, osmextract, units, RANN, curl, ggplot2, patchwork, viridis, jsonlite, spdep). R 4.6 was used for the deposited run; the permutation seed is fixed (`set.seed(2026)`). First full run: several hours and 32 GB RAM (road-graph build and routing); cached reruns complete in minutes. The script stops at Section 1 until the input files above are in place, by design.

## Script-to-manuscript map

Tables 1-2 are assembled from `table1_national.csv`, `tableS_urban_rural.csv`, and the model parameters. Supplement items correspond to these files (all inside `derived_tables.zip`):

| Item | Source table |
|---|---|
| eTable 1 | `table1_national.csv` (cumulative columns) |
| eTable 2 | `tableS_transport_legs.csv` |
| eTable 3 | `tableS_bypass.csv` |
| eTable 4 | `tableS_kpi_decomposition.csv`, `validation_moments.csv` |
| eTable 5 | `tableS_regional.csv` |
| eTable 6 | `table2_governorate.csv` |
| eTable 7 | `tableS_inequality.csv` |
| eTable 8 | `tableS_spatial_moran.csv` |
| eTable 9 | `tableS_lisa_lowlow.csv` |
| eTable 10 | `tableS_sensitivity.csv` |
| eTable 11 | `tableS_trauma_tiers.csv` |
| eTable 12 | `tableS_upgrade_sites_stemi.csv` |
| eTable 13 | `tableS_upgrade_sites_stroke.csv` |
| eTable 14 | `tableS_upgrade_sites_trauma.csv` |
| eTable 15 | `tableS_station_sites.csv` |
| eMethods 1 (broad trauma tier: same-town hospital pairs, straight-line and road km) | `tableS_broad_tier_pairs.csv` |
| eMethods 3 (rank stability) | `tableS_upgrade_stability.csv` |

## Notes

- The stroke Tier-2 (reperfusion-capable) arm activates automatically when a clinician-verified list is supplied (Section 1); the deposited run used Tier-1 receiving designation only, as reported.
- Facility sector labels carry designation as published, not verified clinical capability.
- The calibration factor is fitted to one observed moment (mean travel-to-scene, 7 min 46 s). The residual delay delta is then solved so the modeled 8-minute share equals the observed 45.43%; it is an unexplained residual (operational and model-related), not a measured delay. The reported mean response, the stroke- and STEMI-pathway responses, and the STEMI pathway's mean time to reach a cardiac hospital (37 minutes; compared with modeled medians) serve as consistency checks.
- STEMI access is measured from first medical contact (ambulance arrival) to device, the ACC/AHA metric (90 minutes; 120 minutes as the fibrinolysis threshold): 15 minutes on scene + calibrated transport + an assumed 30-minute door-to-device time (`PAR$door_device_min`; 45 and 63 minutes in sensitivity rows). SRCA's 37-minute STEMI-pathway figure is the time to reach a cardiac hospital and is used only as a consistency check. Call-to-door times are also reported.
- Origins and destinations snap only to the strongly connected core of each road component (`graph_core_vertices.rds`); `data/processed/reachability_check.csv` records that every facility and station receives routed residents or shares a campus with one that does.
- Candidate sites for new stations are populated 1-km cells, treated as locations only.

## License

MIT License

Copyright (c) 2026 Ahmed Batobara, Maha Alfaraidhy, Chris Lim

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
