# 03_analysis — Analysis Pipeline

R scripts that produce every estimate in the paper.

## Label policy

**All output table and figure labels must be in English.** German economic
terminology is translated as follows (enforced in `_dict.R` and every
`ET_DICT` override):

| German term | English label |
|---|---|
| Steuerkraft | Tax capacity |
| Kaufkraft | Purchasing power |
| Grüne (vote share) | Green (vote share) |
| Kreis / Kreise | County / Counties |
| Gemeinde / Gemeinden | Municipality / Municipalities |
| Bundesland | State |
| Stadtstaaten | City-states |
| VZE (Vollzeitäquivalente) | FTE (Full-Time Equivalents) |
| kreis_funded | County-funded | Read the root `CLAUDE.md`
for project-level context (geography, conventions, run order). This file
covers the analysis layer specifically — what each script does, what frame
it consumes/produces, and the econometric choices it makes.

For long-form rationale: `documents/econometric_specification.md`.

## Run order

1. `01_spatial_weights_ags8.py` (Python — produces `spatial_neighbors_ags8.csv`)
2. `00_prep_analysis.R` — builds the six estimation frames
3. Then any of `01_descriptives.R`, `02_hazard.R`, `02b_logit_uncensored.R`,
   `03_did_main.R`, `05_heterogeneity.R`, `06_spillovers.R`
   (`04_did_robustness.R` was deleted; robustness lives in `03_did_main.R`)
4. `07_assemble.R` — gathers everything into a manifest

`_dict.R` is sourced from every script; never run on its own.

## Frames built by `00_prep_analysis.R`

| File | Purpose | Restriction |
|---|---|---|
| `frame_hazard.csv` | Hazard model on direct (Gemeinde-level) EMK onsets | risk set drops post-onset unit-years; broadcast-only units stay in; `kreis_funded` switches on after the AGS5's Kreis-level project starts |
| `frame_hazard_cov.csv` | Hazard model on broad (coverage) onsets | risk set drops post-broad-onset unit-years |
| `frame_logit_full.csv` | Absorbing-logit counterfactual to `frame_hazard.csv` | full panel, no risk-set censoring; post-onset AGS8-years retained |
| `frame_logit_cov_full.csv` | Absorbing-logit counterfactual to `frame_hazard_cov.csv` | full panel; post-onset retained |
| `frame_did_broad.csv` | Staggered DiD with broad treatment definition | full panel with `gname_cs`; all AGS8 retained |
| `frame_did_direct.csv` | Staggered DiD with direct treatment definition | drops `treat_type == "broadcast_only"`; never-treated pool same as broad |

Constants set in `00_prep_analysis.R`: `MIN_YEAR = 2010L`, `BASE_WINDOW =
2014:2016`, `START_YEAR = 2015L` (hazard risk-set start), `POP_MIN = 500L`
(unit-level population filter), `WINSOR_Q = 0.99` (cap on per-100k outcomes).

## Scripts

### `_dict.R` — shared constants and table writers

**Constants:**
- `STADTSTAATEN`: `c("02", "04", "11")`. Dropped from spec (2) of the hazard
  model (the one with `pers_z`) because the city-state structure conflates
  municipal/Länder roles. All other specs keep them.
- `ES_MIN = -4L`, `ES_MAX = 7L`; `es_max_data_driven(dat, yname)` gives the
  per-frame cap = last outcome year minus earliest reachable cohort.
- `XFORMLA_CS = ~ sk_base_z + state_green_base_z + dens_base_z`. Three
  baseline z-covariates for the conditional CS variant. Order matches output
  tables: tax capacity | state Green vote share | log pop. density. Baseline
  BEV/charging covariates excluded (mass-zero across 2014–16 cross-section
  → singular DR design matrix).
- `z(x)`: scale to mean-0 SD-1 on the input sample. Apply *after* every
  sample restriction.
- `wq(x, w, probs)`: population-weighted quantiles.
- `resolve_root()`, `results_dir(root, stem)`: path helpers.

**Table writers** (all regression and descriptive tables use these):
- `write_longtblr(stem, caption, label, note, colspec, header_rows, body_rows,
  footer_rows, rowsep="-3pt")`: core TeX builder using tabularray `longtblr`.
  Emits `label = {}` inside the `[...]` options block.
- `write_coef_longtblr(models, stem, caption, label, note, groups, var_labels,
  event_counts, show_pr2)`: builds a raw-coefficient table from a named list
  of fixest models. Groups are a named list; `NULL` entries insert `\hline`;
  non-NULL entries are separated by blank lines (no `\SetCell` group headers).
- `write_estimates_csv(models, file)`: writes a long-form CSV of
  estimate/SE/CI for every model in the list.

### `00_prep_analysis.R`

- Loads `emk_inkar_panel_ags8.csv` + spatial neighbours.
- Reclassifies pre-`MIN_YEAR` treated units as untreated controls.
- Unit-level POP_MIN filter (drops AGS8 whose min panel-period population
  is below 500).
- Winsorises per-100k outcomes at the 99th percentile.
- Adds `log1p_<rate>` twins for every per-100k outcome.
- Builds `base_dt` (per-AGS8 baseline snapshot, earliest year in 2014-2016);
  z-scores the columns; falls back to `state_gruene` baseline where
  `muni_gruene` is missing.
- Writes the six estimation frames and a cohort table.

### `01_descriptives.R`

Descriptive plots and tables: concentration curves (Lorenz-style), treated vs
control balance, choropleth maps, and VIF/correlation diagnostics. Outputs via
`write_longtblr()`:

| Stem | Content |
|---|---|
| `desc_balance` | Balance table (means + tests, treated vs control) |
| `desc_means` | Descriptive means table |
| `desc_medians` | Descriptive medians table |
| `desc_vif` | VIF table for hazard-frame channel set |

Note: VIF section includes `kk_z` and `muni_gruene_z` as diagnostics even
though neither appears in the main hazard specifications.

### `01_spatial_weights_ags8.py`

Builds AGS8 queen-contiguity neighbours (1st through 3rd order) and computes
per-AGS8 per-year flags for whether any neighbour had an absorbing EMK project.
Output: `spatial_neighbors_ags8.csv` with columns
`emk_absorbing_any_nbrs_{1,2,3}`. Merged into the panel in
`00_prep_analysis.R`.

### `02_hazard.R`

**Discrete-time hazard model** of EMK onset on the `frame_hazard.csv`
(direct-onset risk set). Identification: cross-section + time variation in
baseline covariates predicts who applies/wins funding.

**Specification:** cloglog link (primary) and logit link (appendix robustness),
year + AGS2 fixed effects, SEs clustered on AGS5. All channels lagged one year
(L1) and z-scored on the estimation sample. City-states dropped only in
spec (2), which includes `pers_z`.

**Three-spec grid** (same across cloglog and logit runs):
- (1) Baseline: `sk_z + bev_z + chg_z + state_gruene_z + log_dens_z`
- (2) + Personnel: adds `pers_z`; drops city-states
- (3) + County coverage: adds `kreis_funded`

**Election-variant table** uses spec (3) as the base, varying the Green
party electoral level: state (col. 1), federal (col. 2), municipal (col. 3).
Municipal Green (`muni_gruene_z`) has NAs where no municipal election was held,
so col. 3 runs on a smaller sample.

**Outputs** in `04_results/02_hazard/`:

| Stem | Content | Where |
|---|---|---|
| `hazard_hr_ame` | cloglog HR + AME (7-col, 3 specs) | Main text |
| `hazard_coef` | cloglog raw log-hazard coef (appendix) | Appendix |
| `hazard_logit_hr_ame` | logit OR + AME (7-col) | Appendix |
| `hazard_logit_coef` | logit raw log-odds coef | Appendix |
| `hazard_election_hr_ame` | election-variant HR + AME | Appendix |
| `hazard_election_coef` | election-variant raw coef | Appendix |
| `hazard_broad_hr_ame` | broad-treatment cloglog HR + AME (5-col, 2 specs) | Appendix |
| `hazard_broad_coef` | broad-treatment cloglog raw coef | Appendix |

Note: penalised logit (brglm2) was removed from the pipeline.

HR = `exp(β̂)`, 95% CI = `[exp(β̂ − 1.96·SE), exp(β̂ + 1.96·SE)]`. OR
(logit) uses the same formula. AMEs computed via `avg_slopes(m, vcov = ~AGS5)`
from `marginaleffects`.

### `02b_logit_uncensored.R`

**Absorbing-treatment logit** on `frame_logit_full.csv` (full panel, no
risk-set censoring). Different outcome and sample logic from `02_hazard.R`.
Outcome: `treat_absorb = 1{year >= first_treat}` (absorbing dummy).
Year + AGS2 FE, SEs clustered on county. State Green vote share only (no
federal or municipal election variants).

Specification grid:
- Direct (3 specs): (D1) baseline, (D2) +pers, (D3) +county_funded
- Broad (2 specs): (B1) baseline, (B2) +pers; county-funded excluded
  (near-collinear with broad absorbing dummy)

Tables combine Direct and Broad into a single 5-column layout with spanning
headers indicating the treatment definition:

Outputs in `04_results/02b_logit_uncensored/`:
- `logit_coef.{tex,csv}` — raw coefficients, Direct (1)-(3) | Broad (1)-(2)
- `logit_ame.{tex,csv}` — average marginal effects, same layout
- `logit_diag_{direct,broad}.csv` — treated obs per state/year

### `03_did_main.R` — all DiD estimates

**Callaway-Sant'Anna (CS)** estimator (`did::att_gt`, `est_method = "dr"`).
Bootstrap: multiplier, `CS_BITERS = 2000L`, clustered on AGS5.
Event-study simultaneous bands: `cband = TRUE` in `aggte()`.
`allow_unbalanced_panel = TRUE`.

**Primary outcome:** `bev_neuzulassungen_p100k` (level only; winsorised).
**Two specs per frame:** unconditional (`xformla = NULL`) and conditional
(`xformla = XFORMLA_CS`), producing side-by-side columns in each table.

**Sections:**

| Section | Frame | Control | Outcomes | Files |
|---|---|---|---|---|
| A Main | direct | nevertreated | BEV overall | `es_main_direct.*` |
| B Main | broad | nevertreated | BEV overall | `es_main_broad.*` |
| C Robustness | direct (ever-treated only) | notyettreated | BEV overall | `es_robust_notyet.*` |
| D | direct | nevertreated | Corporate BEV | `es_corp.*` |
| E | direct | nevertreated | Private BEV | `es_priv.*` |
| F | direct | nevertreated | ICE (placebo) | `es_ice.*` |

Each section produces: `<stem>.png` (event-study graph, simultaneous CI, 300 dpi via ragg),
`<stem>.tex` (longtblr event-study table, appendix), `<stem>.csv` (twin).
Plus `est_att_main.csv` (all overall ATTs).

**Pre-treatment test:** Wald χ² on joint pre-period ATTs, covariance from
bootstrap influence functions, reported at the bottom of each table.

**A1 guard:** verifies never-treated sample > 50% of rows and > 1000 units.

**`04_did_robustness.R` has been deleted** — its content is now fully
contained in section C of `03_did_main.R`.

### `05_heterogeneity.R`

Tercile-stratified CS on `bev_neuzulassungen_p100k`. Splits treated and
control units by `kk_base_terc` (Kaufkraft, headline) and `sk_base_terc`
(Steuerkraft, appendix). Plus treat-type heterogeneity on the broad frame.
Uses unconditional CS to match the headline spec.

Outputs via `write_longtblr()`:
- `did_heterogeneity_terciles.{tex,csv}`
- `did_heterogeneity_treat_type.{tex,csv}`

### `06_spillovers.R`

Two SUTVA probes on the broad frame:

1. **Donut:** rerun the headline CS after dropping never-treated controls
   flagged by `emk_absorbing_any_nbrs_1 == 1L`. Survival of the headline ATT
   indicates control contamination isn't driving it.
2. **Descriptive spillover ES:** among never-treated only, pseudo-treatment
   = first year a 1st-order neighbour appears as treated; BJS event study.
   Descriptive only — selection-on-geography caveat; never reported as causal.

### `07_assemble.R`

Globs `04_results/**/est_att_*.csv` and `04_results/**/tab_*.csv` into one
manifest. Cheap; run last; no estimation.

## Economic identifying logic

The paper is structured in two parts.

**Part (i) — who gets funded (selection)?** Hazard model (`02_hazard.R`).
Identification is cross-sectional + time variation in the channels that predict
EMK application/award; documents which Gemeinde characteristics raise the
per-year onset probability. Used as evidence for selection on observables.

**Part (ii) — what does funding do (adoption response)?** Staggered DiD
(`03_did_main.R`, heterogeneity in 05, SUTVA in 06). Treatment is the year
a Gemeinde is first covered by an EMK project. Outcome is BEV new
registrations per 100k (level). Primary estimator: CS unconditional
(`xformla = NULL`); conditional variant (`xformla = XFORMLA_CS`) shown as
a second column in every table for comparison. Robustness (not-yet-treated
control group) is section C of `03_did_main.R`.

**Heterogeneity** (`05`) reports the inequality payoff: does the funded ATT
differ by baseline Kaufkraft tercile?

**Spillovers** (`06`) probe SUTVA: do never-treated controls with treated
neighbours behave differently from clean controls?

## When iterating

- **Console floods kill RStudio responsiveness.** Keep
  `print_details = FALSE` in every `att_gt` call.
- **Bootstrap is the time sink** in CS. Asymptotic-only (`bstrap = FALSE`)
  is fine for iteration; flip back for the final run.
- **`notyettreated` is more expensive than `nevertreated`** because the
  comparison pool changes per period; section C of `03_did_main.R` runs it
  on the ever-treated subsample only. Use `nevertreated` while iterating.