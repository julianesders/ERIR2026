# ERIR2026 — Project Context

Seminar paper: "Communal Projects and EV Charging Inequality." Author: Julian
Esders. Estimates the effect of German federal mobility-concept (EMK) funding
on Battery-Electric-Vehicle (BEV) adoption at the German municipal (Gemeinde)
level, and probes whether the funded subsidy reaches under-resourced
Gemeinden or reinforces existing capacity gaps.

For the long-form write-ups see `documents/pipeline_summary.md` (data prep,
step by step) and `documents/econometric_specification.md` (every regression
in detail). This file is the navigational summary that's always-in-context.

## Repo layout

```
ERIR2026/
├── wd/
│   ├── 01_data/
│   │   ├── 01_raw/                   # raw drops from BBSR, KBA, INKAR, BNetzA
│   │   ├── 02_clean/                 # cleaned per-source
│   │   └── 03_final/                 # merged panel + estimation frames (CSV)
│   ├── 02_code/
│   │   ├── 00_scraping/              # one-off scrapers (ADAC, etc.)
│   │   ├── 01_cleaning/              # per-source Python cleaners
│   │   ├── 02_merging/               # build the AGS8 panel
│   │   └── 03_analysis/              # R analysis pipeline (see its CLAUDE.md)
│   └── 04_results/                   # CSV+TeX twins, figures, by script stem
├── documents/                        # paper drafts, slides, the two long docs
├── materials/                        # source material (PDFs, slides, etc.)
└── README.md
```

## Pipeline at a glance

Python (01_cleaning + 02_merging) → CSV → R (03_analysis) → CSV/TeX/PDF.

```
raw data sources
  ↓ (01_cleaning/*.py — per source: INKAR, KBA, BNetzA, BBSR, elections, EMK)
01_data/02_clean/
  ↓ (02_merging/02_merge_emk_panel_ags8.py — build the unified AGS8 panel)
01_data/03_final/emk_inkar_panel_ags8.csv
  ↓ (03_analysis/00_prep_analysis.R — build estimation frames)
01_data/03_final/frame_hazard.csv          (DIRECT-onset risk set)
01_data/03_final/frame_hazard_cov.csv      (BROAD-onset risk set)
01_data/03_final/frame_logit_full.csv      (full panel, direct, no censoring)
01_data/03_final/frame_logit_cov_full.csv  (full panel, broad, no censoring)
01_data/03_final/frame_did_broad.csv       (panel for staggered DiD, broad)
01_data/03_final/frame_did_direct.csv      (panel for staggered DiD, direct)
  ↓ (03_analysis/{02_hazard, 03_did_main, ...}.R)
04_results/<script_stem>/{*.csv, *.tex, *.pdf}
```

## Conventions

- **Geography.** AGS8 = Gemeinde (~11k units), AGS5 = Kreis (~400 units, used
  for SE clustering), AGS2 = Bundesland (16 units, used for state FEs and
  Stadtstaat exclusion).
- **String dtype on AGS codes.** Always passed as `character` with leading
  zeros preserved (`sprintf("%08d", ...)`); never read as integer.
- **Stadtstaaten** (AGS2 ∈ {02, 04, 11}) are excluded only from hazard
  specifications that include municipal `pers_z` (personnel), because the
  city-state structure conflates municipal/Länder roles. All other hazard
  specs and the DiD frames keep Hamburg/Bremen/Berlin. List lives in
  `wd/02_code/03_analysis/_dict.R::STADTSTAATEN`.
- **Per-100k rates as primary outcomes.** Each is built in the merge as
  `count / xbev * 100_000` with `xbev` = midyear population. log1p twins are
  computed in `00_prep_analysis.R` (`log1p_<rate>`) for proportional-effect
  specifications.
- **Outputs** live in `wd/04_results/<script_stem>/`. Every CSV-of-estimates
  has a TeX twin with the same numbers.
- **Path resolution** boilerplate at the top of every R script reads
  `commandArgs(trailingOnly = FALSE)` to find its own location, so scripts
  can be run from `Rscript` or RStudio without setting `wd`.
- **Z-scoring on the estimation sample.** The `z()` helper in `_dict.R` is
  applied *after* every sample restriction, so each model's regressors are
  standardised on its own estimating sample (NOT on the full panel).
- **Baseline snapshots.** Time-invariant covariates for CS come from the
  per-AGS8 earliest available year in `BASE_WINDOW = 2014:2016`. Built once
  in `00_prep_analysis.R::base_dt`.

## Data sources

- **INKAR (BBSR):** population, density, purchasing power (Kaufkraft),
  fiscal capacity (Steuerkraft), demographic structure.
- **KBA:** vehicle stock and new-registration counts at AGS8, by drivetrain
  (BEV vs ICE) and holder type (corporate vs private). Available 2014-2023
  flows / 2012-2023 stock.
- **BNetzA Ladesäulenregister:** public charging stations + chargepoints at
  AGS8. Reliable from 2017 (Ladesäulenverordnung); pre-2017 under-coverage
  flagged by the merge.
- **BBSR EMK list:** Bewilligungen (project starts) of the federal mobility-
  concept ("Elektromobilitätskonzepte") program. Built into two onset flags:
  `first_treat_direct` (Gemeinde-level project) and `first_treat_broad`
  (covered by any project at Gemeinde *or* Kreis level).
- **German federal/state/municipal elections:** Grüne vote share by AGS8.
- **BBSR territorial reform crosswalk:** AGS8 changes 1995-2023, applied at
  merge time to back-fill consistent AGS8 IDs.

## Run order

For a from-scratch rebuild:

1. `01_cleaning/*.py` — once per data update.
2. `02_merging/02_merge_emk_panel_ags8.py` — outputs `emk_inkar_panel_ags8.csv`.
3. `03_analysis/01_spatial_weights_ags8.py` — outputs `spatial_neighbors_ags8.csv`.
4. `03_analysis/00_prep_analysis.R` — outputs the six estimation frames.
5. `03_analysis/{01,02,03,04,05,06,07}_*.R` — run any in any order; 07
   aggregates the rest into a manifest.

## Key conventions

- Hazard / DiD never share a frame: hazard uses `frame_hazard.csv` (one row
  per AGS8-year in the risk set), DiD uses `frame_did_broad.csv` /
  `frame_did_direct.csv` (full panel with `gname_cs`).
- The DiD frames are CSV, not RDS. All data flow is CSV.
- Charge-point outcomes are excluded from DiD because BNetzA pre-2017
  under-coverage makes pre-treatment trends unreliable; charge points stay
  as a covariate.
- DiD estimator is Callaway-Sant'Anna (CS) only; BJS was dropped.

## Where to look

- `wd/02_code/03_analysis/CLAUDE.md` — analysis pipeline detail, every R
  script's purpose, econometric choices.
- `documents/pipeline_summary.md` — long-form data prep documentation.
- `documents/econometric_specification.md` — long-form spec for every
  regression in the paper.
- `wd/02_code/03_analysis/_dict.R` — shared constants and table-writer
  helpers (`z`, `wq`, `ES_MIN/MAX`, `XFORMLA_CS`, `STADTSTAATEN`, `COHORT_MIN`,
  `MIN_COHORTS_PER_E`, `OUTCOME_LABELS`, `write_longtblr`, `write_coef_longtblr`,
  `write_estimates_csv`, `resolve_root`, `results_dir`).
- `wd/02_code/03_analysis/_did_helpers.R` — shared CS estimation + output
  machinery (`run_spec`, `run_cs`, `cs_pre_test`, `es_graph`,
  `write_es_longtblr`, `emit_section`, `.att_row`); sourced by `03_did_main.R`
  and `06_spillovers.R` so both produce identical-format figures/tables/CSVs.