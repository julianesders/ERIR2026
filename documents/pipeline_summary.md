# ERIR2026 — Pipeline Summary

Schematic walkthrough of every step from raw delivery files to final estimate
tables. Internal handoff style: file paths, input/output, key choices.

Last update: 2026-06-12. Reflects the analysis-plan v2 rewrite.

---

## 0. End-to-end flow

```
RAW SOURCES                  CLEANING                MERGE              ANALYSIS
─────────────                ──────────              ─────              ────────
EMK web (scrape)        ┐
Anschriftenverz. xlsx   │
INKAR 2025 (63M rows)   │    01_cleaning/*.py        02_merging/*.py    03_analysis/*.R
Ladestationenregister   │    (12 scripts)            (2 scripts)
GERDA elections         │      │                       │                  │
KBA delivery txt        ├──►  │   ──►  emk_inkar_panel_ags8.csv  ──►   00_prep_analysis.R
ADAC scrape             │      │                                          │
Genesis area 11111      │      │                                          │
VG250 geometries        │      │                                          ▼
INKAR personnel 74111   ┘      │                                  frame_{hazard,did_*}.csv
                               │                                          │
                                                                          ▼
                                                          01–06 (descriptives, hazard,
                                                          DiD main, DiD robust, hetero,
                                                          spillover) → 07_assemble
```

Runner: `wd/main.py` orchestrates steps 0–14. R analysis is sourced manually
(`Rscript NN_*.R` or RStudio source).

---

## 1. Cleaning (Python, `wd/02_code/01_cleaning/`)

One row per script. Each one runs from `main.py` in numbered order; every
script writes to `01_data/02_intermediate/<topic>/`.

| Step | Script | Purpose (1-liner) | Key output |
|------|--------|--------------------|------------|
| 0  | `00_scraping/scrape_emkonzepte.py`        | Scrape EMK project metadata from Förderkatalog | `01_data/01_raw/scraping_raw_60/` |
| 1  | `01_clean_anschriftenverzeichnis.py`      | Clean Bundesländer's official address list → AGS8 universe | `anschriftenverzeichnis/anschriftenverzeichnis_ags8.csv` |
| 2  | `02_clean_emk.py`                         | Parse EMK metadata; create `tag_*` / `space_*` dummies; standardise start/end dates | `emk/emk_all.csv` |
| 3  | `00_inkar_extract.py`                     | Chunked extract of INKAR indicators (Kuerzel filter) from the 63M-row CSV | `inkar/inkar_joint_panel.csv` |
| 4  | `03_clean_inkar_ags8.py`                  | Pivot wide; broadcast Kreis vars to AGS8 via address xwalk | `inkar/inkar_ags8_panel.csv` |
| 5  | `04_clean_ladestationen_ags8.py`          | Spatial join Ladestationenregister → VG250 → AGS8 counts/year | `ladestationen/ladestationen_ags8_panel.csv` |
| 6  | `05_clean_elections_ags8.py`              | GERDA national/state/muni election results, harmonised to 2023 AGS8 via 2021→23 BBSR crosswalk | `elections/elections_ags8_panel.csv` |
| 7  | `06_build_ags8_base.py`                   | Validate the AGS8 spine vs INKAR | (validation prints only) |
| 8  | `07_clean_personal_ags5.py`               | Munical personnel VZE (Genesis 74111) by AGS5 | `inkar/personal_ags5_panel.csv` |
| 9  | `08_clean_kba_ags8.py`                    | Unpack KBA delivery zips; FWF parse Bestand + Neuzulassung at AGS8 × fuel × ownership | `kba/kba_panel.csv`, `kba/kba_neuz_model_panel.csv` |
| 10 | `09_aggregate_kba_vars.py`                | Reshape KBA to wide AGS8 × year with `{B,N,A}_<fuel>_<own>` and ev-share columns | `kba/kba_ags8_panel.csv` |
| 11 | `10_scrape_adac_prices.py`                | (Optional) scrape ADAC list prices for BEV models | `adac/...` |
| 12 | `11_clean_area_ags8.py`                   | Parse Genesis 11111 Gebietsfläche; reindex to full year range; ffill+bfill | `area/area_ags8_panel.csv` |

KBA boundary harmonisation (P0.2 in the plan) is deferred — BBSR annual files
2013–2020 not yet available.

---

## 2. Merging (Python, `wd/02_code/02_merging/`)

### 2.1 `01_match_emk_ags.py`

Fuzzy-match EMK project metadata to AGS8 / AGS5 (where the Förderkatalog
provides only a Stadt/Kreis name).

| Output | Content |
|---|---|
| `emk/emk_ags_matched.csv` | One row per project; `AGS8` filled when match is unambiguous, else NA; `AGS5` filled for all rows |

Direct = 187 projects with non-NA `AGS8`. Kreis-broadcast = 54 projects with
NA `AGS8` but matched `AGS5`.

### 2.2 `02_merge_emk_panel_ags8.py`

Build the canonical panel. Output: `01_data/03_final/emk_inkar_panel_ags8.csv`.

#### Inputs and joins

| Side | Source | Key |
|---|---|---|
| Spine | INKAR panel (AGS8 × year) | `(AGS8, year)` |
| Area  | `area/area_ags8_panel.csv` | `(AGS8, year)` |
| Ladestationen | `ladestationen/...` | `(AGS8, year)` |
| KBA | `kba/kba_ags8_panel.csv` (count cols + shares + N_benzin_diesel_overall) | `(AGS8, year)` |
| Elections | `elections/elections_ags8_panel.csv` | `(AGS8, year)` |
| Personnel | `inkar/personal_ags5_panel.csv`, divided by AGS5 population (per-100k) | `(AGS5, year)` |
| Treatment map | derived from `emk_ags_matched.csv` (see below) | `AGS8` |

#### Treatment derivation

```
first_treat_direct  := min(start_year) over projects with non-NA AGS8 matching unit
first_treat_broad   := min(start_year) over (direct ∪ Kreis-broadcast) matching unit
kreis_funded_year   := min(start_year) over AGS5-level projects in own Kreis
treat_type          := "direct"          if first_treat_direct notna
                    := "broadcast_only"  elif first_treat_broad  notna
                    := "never"           otherwise
```

#### Derived columns

| Column | Definition |
|---|---|
| `log_pop_dens` | `log(xbev / area_qkm)` |
| `log_steuerkraft` | `log1p(clip0(q_gest_bev))` (negatives clipped to 0) |
| `log_kaufkraft` | `log(q_kaufkraft)` (fallback log1p if any non-positive) |
| `bev_{stock,neuzulassungen,corporate,private}_p100k` | KBA counts ÷ xbev × 100,000 |
| `ice_neuzulassungen_p100k` | `N_benzin_diesel_overall ÷ xbev × 100,000` (placebo) |
| `eco_index` | PC1 of log1p(`bev_stock_p100k`, `ev_chargepoints_p100k`), fit on year≥2017 only (BNetzA registry pre-2017 under-coverage), sign-flipped to "more is more EV ecosystem"; scored for every row where both inputs are present |
| `*_imp` | Boolean flag for KBA cells that were linearly interpolated (`limit_area="inside"`) |

#### Lag block

Year-based merge to produce L1/L2/L3 of: `log_steuerkraft`, `log_kaufkraft`,
`eco_index`, `bev_stock_p100k`, `ev_chargepoints_p100k`,
`{fed,state,muni}_gruene`, `n_vze_personal`, `N_elektro_{overall,private,corporate}`,
`N_benzin_diesel_overall`, `N_ev_share_*`, plus raw `q_gest_bev`.

---

## 3. Spatial structure (`wd/02_code/03_analysis/01_spatial_weights_ags8.py`)

Builds neighbour indicators for the spillover analyses.

| Layer | Geometry | Treatment def. | Output columns |
|---|---|---|---|
| Granular  | Queen contiguity at **AGS8** (Gemeinde) | DIRECT (`first_treat_direct ≤ year`) | `n_nbrs_gem_{1,2,3}`, `direct_treated_nbrs_gem_{1,2,3}`, `direct_treated_any_nbrs_gem_{1,2,3}` |
| Aggregated | Queen contiguity at **AGS5** (Kreise, after dissolve) | BROAD per Kreis (any AGS8 with `first_treat_broad ≤ year`) | `n_nbrs_kreis`, `broad_treated_nbrs_kreis`, `broad_treated_any_nbrs_kreis` |

Output: `01_data/03_final/spatial_neighbors_ags8.csv`. One row per
`(AGS8, year)`.

---

## 4. Analysis (R, `wd/02_code/03_analysis/`)

All scripts source `_dict.R` (label map, ES horizons, `z()`, `wq()`, output dir
helpers). All write to `04_results/<NN_scriptname>/`. CSV+TeX twins for every
estimate table.

### 4.1 `00_prep_analysis.R` — Estimation frames

| Step | Detail |
|---|---|
| Load | `emk_inkar_panel_ags8.csv` + `spatial_neighbors_ags8.csv` |
| Baseline snapshot | per AGS8, earliest year in {2014,2015,2016} with `log_kaufkraft / log_steuerkraft / log_pop_dens / bev_stock_p100k / ev_chargepoints_p100k / muni_gruene / xbev` observed |
| Tercile / quintile cuts | unweighted + population-weighted on `kk_base`, `sk_base` |
| z-scoring | for `_z` suffix on baseline snapshot (used as time-invariant covs in CS-dr) |
| Frames | see table below |
| Cohort table | onsets per year × {broad, direct}; pre/post-period availability |

#### Frames written to `01_data/03_final/`

| File | Risk-set / sample | Event variable | Notes |
|---|---|---|---|
| `frame_hazard.csv`       | `year ≥ 2015`, drop unit-years after `first_treat_direct` | `onset_direct` | Broadcast-only units stay in risk set with `kreis_funded` switching on |
| `frame_hazard_cov.csv`   | `year ≥ 2015`, drop unit-years after `first_treat_broad`  | `onset_broad`  | Appendix |
| `frame_logit_full.csv`   | full panel (no risk-set censoring); post-onset rows have mechanical-zero `onset_direct` | `onset_direct` | For absorbing-treatment logit in `02b_logit_uncensored.R` |
| `frame_logit_cov_full.csv` | full panel; post-onset rows have mechanical-zero `onset_broad` | `onset_broad` | Written by `00_prep_analysis.R` (currently unused) |
| `frame_did_broad.csv`    | full panel | `gname_cs` (0 = never) / `gname_bjs` (0 = never) from `first_treat_broad` | Carries baseline snapshot + spatial cols |
| `frame_did_direct.csv`   | full panel, broadcast-only units DROPPED | g-names from `first_treat_direct` | Sharp identification |

`MAX_COHORT = NA` for first pass (no policy-regime cap).

### 4.2 `01_descriptives.R` — Inequality storyline (access, not euros)

| Output | Content |
|---|---|
| `fig_concentration_coverage_sk.pdf`  | Lorenz-style curve, AGS8 ranked by `sk_base` (Steuerkraft) |
| `fig_concentration_multi.pdf`        | Multi-ranking Lorenz curves (direct treatment, 6 ranking variables) |
| `tab_concentration_index.csv`        | Concentration indices (sk broad/direct + multi-ranking direct) |
| `tab_balance.{tex,csv}`              | Means by `treat_type` + normalized differences vs "never" |
| `tab_desc_means.{tex,csv}`           | Direct vs Never: unit means + t-test p-values (2015–latest) |
| `tab_desc_medians.{tex,csv}`         | Direct vs Never: unit medians + Mood's median test p-values |
| `fig_map_treatment.pdf`              | VG250 choropleth of `treat_type` (optional; needs sf + gpkg) |
| `tab_corr.csv`                       | Pairwise correlations on hazard-frame channels |
| `tab_vif.{tex,csv}`                  | VIFs on hazard-frame channels |
| `fig_map_bev_dual_2023.png`          | BEV stock p100k + BEV share of stock, 2023 (optional) |

### 4.3 `02_hazard.R` — Part (i): drivers of EMK onset

Estimator: `fixest::feglm` with `family = binomial("cloglog")`, FE `| year + AGS2`,
SEs clustered at AGS2. Logit twin for robustness. All channels lag-1 (z'd on
estimation sample). `log_dens_z` is a universal control in all columns.

#### Columns

| # | Spec | Sample |
|---|---|---|
| (1) | `~ sk_z + bev_z + chg_z + state_gruene_z + log_dens_z` | full |
| (2) | (1) + `pers_z` | no-Stadtstaaten (drop AGS2 ∈ {02,04,11}) |
| (3) | `sk_z + eco_z + pers_z + state_gruene_z + log_dens_z` | full |
| (4) | (1) replacing `state_gruene_z` with `fed_gruene_z` | full |
| (5) | (1) replacing `state_gruene_z` with `kreis_funded` | full |

#### Outputs

| File | Content |
|---|---|
| `tab_hazard_coef.{tex,csv}`  | Headline cloglog coefficients (with clustered SEs) |
| `tab_hazard_logit.{tex,csv}` | Logit twin |
| `tab_hazard_ame.{tex,csv}`   | Point AMEs via `avg_slopes(vcov=FALSE)` (SEs only on the linear-index coefs) |
| `tab_hazard_robust.csv`      | brglm2 penalized logit + LPM + cloglog complete-case |
| `tab_hazard_diag.csv`        | Events per AGS2 / per year |
| `tab_hazard_cov.{tex,csv}`   | Appendix: coverage-event variant on `frame_hazard_cov.csv` |
| `fig_coef_stability.pdf`     | (1) vs (3) coefficient comparison |

### 4.4 `03_did_main.R` — Part (ii): adoption response

Two co-primary staggered-DiD estimators run for every (outcome × design) cell:

- **CS** (`did::att_gt` → `aggte`), `est_method = "dr"`, `xformla = NULL`
  (unconditional parallel trends for the headline; conditional variant is in `04`),
  `bstrap=TRUE, biters=2000, clustervars="AGS5", allow_unbalanced_panel=TRUE`.
- **BJS** (`didimputation::did_imputation`), defaults; `cluster_var="AGS5"`.

#### Design grid

| Role | Frame | CS control | Anticipation |
|---|---|---|---|
| Main      | broad  | `nevertreated`   | 0 |
| Sharp     | direct | `nevertreated`   | 0 |
| Selection | broad  | `notyettreated`  | 1 |

#### Outcomes

`bev_neuzulassungen_p100k` (headline), `bev_corporate_p100k`,
`bev_private_p100k`, `bev_stock_p100k` (secondary, cumulative),
`ice_neuzulassungen_p100k` (placebo). Charging-points is *not* an outcome in
this paper; it is a covariate/baseline control only.

#### Outputs

| File | Content |
|---|---|
| `est_att_main.csv`         | Long form: outcome × design × estimator (BJS, CS-dr) |
| `tab_att_main.{tex,csv}`   | Wide form for main + sharp roles |
| `pretests.csv`             | CS Wald pre-test W, p-value (BJS pretrends in the plot leads) |
| `es_<outcome>_<frame>_<ctrl>.pdf` | Event studies (BJS + CS overlay, ribbon = 95% CI, dashed line at -0.5) |

### 4.5 `04_did_robustness.R` — Sensitivity

Loop over variants × outcomes × estimators.

#### Variants

| Name | Filter | Anticipation | xformla | Control |
|---|---|---|---|---|
| `baseline`           | none                                            | 0 | NULL         | nevertreated |
| `conditional`        | none                                            | 0 | `XFORMLA_CS` | nevertreated |
| `anticipation_1`     | none                                            | 1 | NULL         | nevertreated |
| `drop_2016`          | drop cohort 2016 (1 pre-period)                 | 0 | NULL         | nevertreated |
| `drop_covid`         | drop `year ∈ {2020, 2021}` (short-run only)     | 0 | NULL         | nevertreated |
| `complete_case`      | drop rows where `N_elektro_overall_imp == TRUE` | 0 | NULL         | nevertreated |
| `notyet_evertreated` | ever-treated only                               | 0 | NULL         | notyettreated |

#### Estimators

Each variant runs BJS (static, `horizon=NULL`) and CS (simple ATT, `biters=999`).

Outputs: `est_att_robust.csv`, `fig_robust_grid.pdf` (forest of ATTs).

### 4.6 `05_heterogeneity.R` — Inequality payoff

| Block | Implementation | Output |
|---|---|---|
| Tercile-stratified CS | Split sample on `kk_base_terc` / `sk_base_terc`; rerun CS-dr (never-ctrl) per tercile; report dynamic and simple aggregations | `tab_att_terciles.{tex,csv}`, `fig_att_by_tercile.pdf` |
| Top–bottom diff (KK)  | Normal-approx CI on `ATT(top) − ATT(bottom)` | `tab_diff_top_bottom_kk.csv` |
| Treat-type heterogeneity | CS on broad frame restricted to `treat_type ∈ {direct, never}` vs `{broadcast_only, never}` | `tab_att_treat_type.{tex,csv}` |

### 4.7 `06_spillovers.R` — Spatial spillovers

| Block | Detail | Output |
|---|---|---|
| Donut robustness | Drop never-treated controls with `direct_treated_any_nbrs_gem_1 == 1` OR `broad_treated_any_nbrs_kreis == 1`. Rerun CS-dr main spec | `est_att_donut.csv` |
| Descriptive ES   | Restrict to never-treated AGS8; pseudo-treatment = first year `direct_treated_any_nbrs_gem_1 == 1`; BJS event study | `es_spillover.pdf`, `es_spillover.csv` |

### 4.8 `07_assemble.R` — Manifest

Walks `04_results/` and writes `manifest.csv`: file, script_dir, ext, mtime,
size, content type. No re-computation.

---

## 5. Treatment-frame matrix

| Frame | Sample | Treatment | Used by |
|---|---|---|---|
| `frame_hazard`      | year ≥ 2015, pre-onset only | onset of `first_treat_direct` | 02_hazard.R |
| `frame_hazard_cov`  | year ≥ 2015, pre-onset only | onset of `first_treat_broad`  | 02_hazard.R (appendix) |
| `frame_did_broad`   | full panel                  | g-names from `first_treat_broad`  | 03–06 (main, robust, hetero, spillover) |
| `frame_did_direct`  | drops broadcast-only units  | g-names from `first_treat_direct` | 03_did_main.R (sharp role) |

---

## 6. Lag / baseline policy (single source of truth)

| Variable | Timing | Transform | Scaled |
|---|---|---|---|
| `log_steuerkraft` | t–1 | `log1p(clip0)` (in merge) | z on est. sample |
| `log_kaufkraft`   | t–1 | `log`                     | z |
| `bev_stock_p100k` | t–1 | `log1p`                   | z |
| `ev_chargepoints_p100k` | t–1 | `log1p`             | z |
| `n_vze_personal`  | t–1 | `log1p` (after p100k)     | z (no-Stadtstaaten sample) |
| `muni/state/fed_gruene` | t–1 (LOCF election) | none | z |
| `eco_index`       | t–1 | PCA output                | as-is |
| `log_pop_dens`    | t   | `log`                     | z |
| `kreis_funded`    | strict past: `1{kreis_funded_year < t}` | dummy | — |
| DiD baseline covs | snapshot 2014–2016 (per-var earliest) | log1p where applicable | z (cross-section) |

---

## 7. Conventions

- AGS as zero-padded strings: AGS8 (8), AGS5 (5), AGS2 (2). Always
  `colClasses=list(character=c("AGS8","AGS5","AGS2"))` in R; `dtype=str` +
  `.str.zfill()` in Python.
- Output dir: `04_results/<NN_scriptname>/`.
- CSV+TeX twin for every estimate table.
- z-scoring helper: `z <- function(x) as.numeric(scale(x))`, applied AFTER
  sample restrictions.
- Every script prints sample sizes and event counts.
- R path resolution: `argv/self_flag` block (works both `Rscript` and
  RStudio source).

---

## 8. Open / parked

| Item | Status |
|---|---|
| KBA boundary harmonisation to 2023 (P0.2) | Parked — need BBSR ref-gemeinden 2013–2020 |
| `modellregion_pre2015` Kreis dummy | TODO placeholder column (`= 0`) in merge; will plug in once data arrives |
| `MAX_COHORT` policy-regime cap | `NA` (no cap) for first pass — toggle in `00_prep_analysis.R` |
| BJS joint pretrend p-value | Not in the package; using CS Wald + per-term BJS leads in plots |
| `kreis_funded` in hazard col (4) | Implemented as strict-past `1{kreis_funded_year < year}` |
