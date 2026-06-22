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

`_dict.R` is sourced from every script; never run on its own. `_did_helpers.R`
(shared CS estimation + output machinery) is sourced by `03_did_main.R` and
`06_spillovers.R`; also never run on its own.

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
  per-frame display cap = the largest event time `e` for which at least
  `MIN_COHORTS_PER_E` eligible cohorts still contribute an `ATT(g, g+e)`,
  computed on the post-drop frame (no longer `last_year − earliest_cohort`).
- `COHORT_MIN = 5L`: cohorts with fewer than 5 treated units are dropped
  *entirely* from the DiD frames in `00_prep_analysis.R` — singleton /
  near-singleton cohorts give rank-deficient within-cohort covariance and
  unstable DR propensity fits.
- `MIN_COHORTS_PER_E = 3L`: cohort-count floor used by `es_max_data_driven()`
  to cap the dynamic event-study horizon (lower to 2 for a more permissive
  right tail).
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

### `_did_helpers.R` — shared CS estimation + output machinery

Sourced by **both** `03_did_main.R` and `06_spillovers.R` (after `_dict.R`), so
every CS event study — main or spillover — produces the same-format figure,
table and CSV; only the estimation sample differs. Definitions only, no
top-level estimation. Contents:

- Constants: `CS_BITERS = 2000L`, `OUTCOME_{BEV,CORP,PRIV,ICE}`, `COL_LABELS`
  (column key → display label, incl. `donut` / `spillover`), plot theme
  (`PLOT_BLUE`, `PLOT_FILL`, `theme_es`), `.read_frame`.
- `run_cs(yname, dat, xformla, control_group, est_method, gname)`: the
  `att_gt` call. `gname` defaults to `"gname_cs"`; pass a different cohort
  column (e.g. a pseudo-treatment) to reuse the machinery on another design.
- `cs_es_agg` / `cs_att_agg` / `cs_pre_test`: dynamic + simple aggregations and
  the joint-Wald pre-test (`cs$W`/`cs$Wpval`).
- `run_spec(dat, yname, xformla, control_group, est_method, label, gname)`:
  one estimation cell → `{es_agg, att_agg, pre_tst, n_treated, n_control}`.
- `es_graph` (ragg PNG, simultaneous bands), `write_es_longtblr` (longtblr +
  CSV twin), `.att_row` (summary row), `emit_section` (figure + main table +
  optional dr appendix twin).

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
- **Drops small cohorts (`< COHORT_MIN`) from each DiD frame independently**
  (direct and broad have different cohort sizes). Units in dropped cohorts
  are removed entirely (their rows have `gname_cs > 0`, so the never-treated
  pool is untouched; they are NOT reclassified as controls). Logs the dropped
  cohort years/units/rows and asserts `gname_bjs` stays in sync with
  `gname_cs`. For the current data this removes direct cohorts
  {2018, 2020, 2021, 2023} (survivors {2016, 2017, 2019, 2022}) and broad
  singletons {2020, 2023}.
- Writes the six estimation frames and a cohort table. The cohort table is
  built *before* the drop (documents the full structure) and carries a
  `dropped` boolean flagging cohorts below `COHORT_MIN`.

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

**Callaway-Sant'Anna (CS)** estimator (`did::att_gt`); `est_method` is
threaded per spec. Bootstrap: multiplier, `CS_BITERS = 2000L`, clustered on
AGS5. Event-study simultaneous bands: `cband = TRUE` in `aggte()` (pointwise
1.96 fallback when the simultaneous critical value is unavailable).
`allow_unbalanced_panel = TRUE`.

**Primary outcome:** `bev_neuzulassungen_p100k` (level only; winsorised).
**Three specs per main frame:** unconditional (`xformla = NULL`,
estimator-invariant — reg ≡ dr, reported once) and conditional
(`xformla = XFORMLA_CS`) run with **both** `est_method = "reg"` (outcome
regression) and `est_method = "dr"` (doubly robust). The **main table** shows
*Unconditional + Conditional (reg)*; the **appendix twin** (`*_dr.{tex,csv}`)
shows *Conditional (dr)*. The event-study **figure** shows all three columns.
The main `.csv` is the union of all estimated columns (no estimate lost to the
display split). A `dr` overlap diagnostic prints treated-N per surviving cohort
for every `dr` conditional cell and warns when a contributing cohort has `< 10`
treated.

**Sections:**

| Section | Frame | Control | Columns | Files |
|---|---|---|---|---|
| A Main | direct | nevertreated | uncond, cond(reg), cond(dr) | `es_main_direct.*`, `es_main_direct_dr.*` |
| B Main | broad | nevertreated | uncond, cond(reg), cond(dr) | `es_main_broad.*`, `es_main_broad_dr.*` |
| C Robustness | direct (ever-treated only) | notyettreated | uncond, reg, dr | `es_robust_notyet.*`, `es_robust_notyet_dr.*` |
| D | direct | nevertreated | Corporate BEV (uncond) | `es_corp.*` |
| E | direct | nevertreated | Private BEV (uncond) | `es_priv.*` |
| F | direct | nevertreated | ICE placebo (uncond) | `es_ice.*` |

Each section produces: `<stem>.png` (event-study graph, all estimable columns,
simultaneous CI, 300 dpi via ragg), `<stem>.tex` (longtblr event-study table),
`<stem>.csv` (twin). Conditional sections also produce the `<stem>_dr.{tex,csv}`
appendix twin. Plus `est_att_main.csv` — one row per estimated cell with an
`est_method` column (uncond rows appear once, not duplicated reg/dr).

**Pre-treatment test:** the `did` package's own joint Wald χ² of H₀ that all
pre-treatment ATT(g,t) = 0 (`cs$W` / `cs$Wpval`, AGS5-clustered multiplier
bootstrap; df = number of pre-treatment (g,t) cells). Reported as `--` when
unavailable — **no** hand-rolled statistic from the stored event-study
influence function (its raw covariance is ~2 orders of magnitude too small to
reproduce the clustered SEs) and **no** diagonal-sum-of-squares fallback. The
CSV records `pre_test_method` (`joint_wald` / `unavailable` / `no_pre_periods`)
so each row's source is auditable.

**Display horizon:** `es_max_data_driven()` caps `max_e` at the largest `e`
supported by `≥ MIN_COHORTS_PER_E` cohorts (evaluates to 4 for the direct
frame post-drop at the default `MIN_COHORTS_PER_E = 3`).

**A1 guard:** verifies the (post-drop) estimation sample has never-treated
> 50% of rows and > 1000 units.

**`04_did_robustness.R` has been deleted** — its content is now fully
contained in section C of `03_did_main.R`.

### `05_heterogeneity.R`

Tercile-stratified CS on `bev_neuzulassungen_p100k`. Splits treated and
control units by `sk_base_terc` (Steuerkraft / tax capacity) — the **only**
wealth ranking; Kaufkraft is deliberately excluded from the DiD heterogeneity.
Plus treat-type heterogeneity on the broad frame. Uses unconditional CS to
match the headline spec.

Outputs via `write_longtblr()`:
- `did_heterogeneity_terciles.{tex,csv}`
- `tab_diff_top_bottom_sk.csv` (top–bottom Steuerkraft tercile difference)
- `did_heterogeneity_treat_type.{tex,csv}`

### `06_spillovers.R`

Two SUTVA probes on the **direct frame** (sharp Gemeinde-level treatment —
adjacency-based spillover is mechanically coherent there; the broad frame is
not used for spillovers because broad treatment is Kreis-level coverage and
does not map onto Gemeinde adjacency):

1. **Donut:** rerun the headline direct CS (direct frame, never-treated control,
   unconditional — replicates `03_did_main.R` section A exactly on the full
   frame) after dropping never-treated controls whose 1st-order Gemeinde
   neighbour was directly treated (`direct_treated_any_nbrs_gem_1 == 1L`).
   Survival of the headline ATT indicates control contamination isn't driving
   it. (Realised: drops ≈1,028 controls, keeps ≈6,064; simple ATT 30.6 → 34.3
   per 100k — effect survives.)
2. **Descriptive spillover ES:** among never-direct-treated only,
   pseudo-treatment = first year a 1st-order Gemeinde neighbour appears as
   directly treated (`direct_treated_any_nbrs_gem_1`); **CS** event study
   (CS only — BJS was dropped pipeline-wide), dynamic aggregation with
   simultaneous bands. Descriptive only — selection-on-geography caveat;
   never reported as causal.

Both probes run through the **shared CS machinery** in `_did_helpers.R`
(`run_spec` → `emit_section`), so their outputs are byte-identical in format to
the main DiD sections — only the sample differs. Outputs in
`04_results/06_spillovers/`:

| Stem | Content |
|---|---|
| `es_donut.{png,tex,csv}`     | Donut event study (single column, main-spec format) |
| `es_spillover.{png,tex,csv}` | Descriptive spillover event study (same format) |
| `est_att_spillover.csv`      | Overall ATTs (same schema as `est_att_main.csv`) |

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