# ERIR2026 — Econometric Specification

Detailed methods document. Pairs with `pipeline_summary.md` (which is the
implementation walkthrough). Last update: 2026-06-12.

---

## 0. Notation

Indices

| Symbol | Meaning |
|---|---|
| $i$ | Gemeinde, identified by AGS8 |
| $k(i)$ | Kreis containing Gemeinde $i$, identified by AGS5 |
| $s(i)$ | Bundesland of $i$, identified by AGS2 |
| $t$ | Calendar year |
| $g_i$ | First-treatment year of $i$ ("cohort") |
| $e = t - g_i$ | Event time relative to first treatment |

Treatment variants

| Symbol | Definition |
|---|---|
| $g_i^{\text{dir}}$ | First year $i$ has a *direct* AGS8-matched EMK project; $\infty$ (or NA) if never |
| $g_i^{\text{broad}}$ | First year $i$ is under *any* EMK coverage (direct $\cup$ Kreis-broadcast); $\infty$ if never |
| $\text{treat\_type}_i \in \{\text{direct, broadcast\_only, never}\}$ | Classification used as a grouping variable |
| $\text{kreis\_funded}_{i,t}$ | $\mathbf{1}\{$ Kreis $k(i)$ has a Kreis-level project with start $< t\}$ |

Outcomes

| Symbol | Definition |
|---|---|
| $Y_{it}^{\text{BEV}}$ | BEV new registrations per 100,000 inhabitants in Gemeinde $i$, year $t$ |
| $Y_{it}^{\text{corp}}$, $Y_{it}^{\text{priv}}$ | Sub-decomposition by holder type |
| $Y_{it}^{\text{stock}}$ | BEV stock per 100,000 |
| $Y_{it}^{\text{charge}}$ | Public charge points per 100,000 (intermediate channel) |
| $Y_{it}^{\text{ICE}}$ | ICE new registrations per 100,000 (placebo) |

Time-varying covariates $X_{it}$, time-invariant baseline covariates $W_i$
(snapshot at earliest 2014–2016 year per variable) — defined in §3.4.

---

## 1. Identification and ex-ante framing

The paper splits into two complementary questions:

- **Part (i) — Selection.** Conditional on remaining at risk, what observable
  Gemeinde features predict EMK direct onset? This is a discrete-time
  proportional-hazards question; the parameter of interest is the
  conditional hazard ratio of each channel.

- **Part (ii) — Adoption response.** Among funded Gemeinden, what is the
  causal effect of EMK on BEV adoption (and on the charging mechanism that
  rationalises any effect)? This is a staggered DiD; the parameter of
  interest is $\text{ATT}(g, t)$ and aggregations thereof.

The two questions are linked by the inequality framing: if selection in
Part (i) is monotone in baseline wealth (high-Kaufkraft / high-Steuerkraft
Gemeinden are more likely treated) AND the Part-(ii) effect is positive
across the wealth distribution, then EMK funding amplifies the rich-vs-poor
gap in BEV adoption. The two-by-two analysis (selection $\times$ effect by
tercile) is the inequality payoff.

---

## 2. Part (i) — Discrete-time hazard of EMK onset

### 2.1 Risk set and sample (`frame_hazard.csv`)

Indexed unit-years $(i, t)$ surviving to the start of year $t$.

| Filter | Rule |
|---|---|
| Period          | $t \ge 2015$ |
| Risk set         | drop $(i,t)$ with $t > g_i^{\text{dir}}$ (absorbing exit at onset) |
| Broadcast units  | KEPT in the risk set: they never have a direct event, but $\text{kreis\_funded}_{i,t}$ may turn on |
| Complete-cases   | drop $(i,t)$ missing any L1 channel covariate |

Notes:
- Stadtstaaten (AGS2 $\in \{02, 04, 11\}$) are flagged via `ns_flag` but kept
  in the main sample; column (2) below uses the complement.
- Realised sample (current data): $N \approx 85{,}302$ unit-years from
  $\approx 9{,}607$ AGS8, with $\approx 166$ direct onsets (rate $\approx 0.19\%$).

A second frame `frame_hazard_cov.csv` replaces $g_i^{\text{dir}}$ with
$g_i^{\text{broad}}$; used in an appendix table on "determinants of any
coverage."

### 2.2 Estimating equation

Let $h_{it}$ be the conditional hazard of direct onset in year $t$ for unit
$i$ surviving to $t$. Primary specification: complementary log-log link.

$$
\log\bigl(-\log(1 - h_{it})\bigr) = \alpha_t + \delta_{s(i)} + X_{it}'\beta
$$

with year FE $\alpha_t$ and Bundesland FE $\delta_{s(i)}$.

Equivalently, $h_{it} = 1 - \exp\bigl(-\exp(\alpha_t + \delta_{s(i)} + X_{it}'\beta)\bigr)$.

The cloglog link is the discrete-time analogue of the Cox proportional-hazards
model: under it, $\exp(\beta_k)$ is the hazard ratio for a one-unit increase
in $X^{(k)}$. Logit is reported as a robustness twin
($\log\frac{h_{it}}{1-h_{it}} = \alpha_t + \delta_{s(i)} + X_{it}'\beta$).
Year FE supply the non-parametric baseline hazard; AGS2 FE absorb
time-invariant state-level heterogeneity (Förderkulturen, party landscape,
demographics).

#### Channel set $X_{it}$ (all lag-1, $z$-scored on the estimation sample)

| Symbol | Channel | Construction |
|---|---|---|
| `log_dens_z`     | (contemporaneous) log pop. density | $z\{\log(\text{xbev}_{it} / \text{area}_i)\}$ |
| `sk_z`           | Steuerkraft per capita (L1) | $z\{\log1p(\text{q\_gest\_bev}_{i,t-1})\}$ |
| `bev_z`          | BEV stock per 100k (L1) | $z\{\log1p(\text{bev\_stock\_p100k}_{i,t-1})\}$ |
| `chg_z`          | Charging points per 100k (L1) | $z\{\log1p(\text{ev\_chargepoints\_p100k}_{i,t-1})\}$ |
| `state_gruene_z` | State Grüne vote share (L1, LOCF) | $z\{\text{state\_gruene}_{i,t-1}\}$ |
| `fed_gruene_z`   | Federal Grüne vote share (L1, LOCF) | $z\{\text{fed\_gruene}_{i,t-1}\}$ |
| `pers_z`         | Personnel VZE per 100k (L1) | $z\{\log1p(\text{n\_vze\_personal}_{i,t-1})\}$ |
| `eco_z`          | PCA composite (L1) | First PC of $(\log1p \text{bev}, \log1p\text{chg})$, fit on year $\ge 2017$ (BNetzA registry pre-2017 under-coverage), sign-flipped to positive |
| `kreis_funded`   | Strict-past Kreis funding dummy | $\mathbf{1}\{\text{kreis\_funded\_year}_i < t\}$ |

### 2.3 Specification grid

`log_dens_z` is a universal control present in all columns.

| Col | Equation | Sample |
|---|---|---|
| (1) | $X = \{\text{sk\_z, bev\_z, chg\_z, state\_gruene\_z, log\_dens\_z}\}$ | full |
| (2) | (1) $+$ `pers_z` | drop Stadtstaaten (AGS2 $\in \{02,04,11\}$) |
| (3) | `sk_z + eco_z + pers_z + state_gruene_z + log_dens_z` | full |
| (4) | (1) replacing `state_gruene_z` with `fed_gruene_z` | full |
| (5) | (1) replacing `state_gruene_z` with `kreis_funded` | full |

### 2.4 Inference

Robust score-type standard errors clustered at AGS2 (Bundesland)
(`cluster = ~AGS2` in `fixest::feglm`). Treatment at AGS8 level with AGS2
FE means no collinearity issue. Bundesland-level clustering absorbs the
Förderkultur and correlated application behaviour within states.

### 2.5 Average marginal effects (AMEs)

Because the cloglog link is non-linear, $\beta_k$ is the partial effect on
the linear index, not on $h$. The AME on the hazard scale is

$$
\overline{\text{AME}}_k = \frac{1}{N}\sum_{i,t} \beta_k \cdot \exp\bigl(\eta_{it}\bigr) \cdot \exp\bigl(-\exp(\eta_{it})\bigr),
$$

with $\eta_{it} = \alpha_t + \delta_{s(i)} + X_{it}'\beta$. Computed via
`marginaleffects::avg_slopes(model, vcov = FALSE)`. The `vcov = FALSE`
argument is necessary because `marginaleffects` does not propagate
fixest's clustered VCOV through the delta method on FE models; we report
the point AME and rely on the clustered coefficient SEs in
`tab_hazard_coef` for inference.

### 2.6 Robustness

| Variant | Engine | What it tests |
|---|---|---|
| Logit twin                 | `feglm(..., family=logit)`            | Link sensitivity |
| Penalised logit            | `brglm2::brglm_fit` with explicit year+AGS2 dummies | Rare-event bias (Firth-type penalty) |
| LPM                        | `feols(...)`                          | Linearity / specification stability |
| Complete-case              | rows with `N_elektro_overall_imp == FALSE` | Sensitivity to KBA interpolation |
| AGS5-level appendix hazard | manual, on $\approx 50$ Kreis events | Different unit of treatment |

### 2.7 Identification assumptions (Part i)

The hazard interpretation requires:
1. **Discrete-time proportional hazards**: $\beta_k$ is invariant in $t$ and
   $s$. Logit/LPM twins probe the link choice.
2. **No selection conditional on observables**: given $X_{it}$, the
   unobservable determinants of $g_i^{\text{dir}}$ are independent of $t$
   beyond what year FE absorb. Implicitly, the Förderkultur and unobserved
   policy entrepreneur networks are state-constant (captured by AGS2 FE).
3. **Correct timing**: all channels are L1 (year $t-1$) except
   `log_pop_dens` (contemporaneous; treated as predetermined within-year)
   and `kreis_funded` (strict past).

Part (i) is descriptive of conditional selection — not a causal estimate of
"what raises a Gemeinde's chance of being funded." The take-away is which
covariates correlate with onset, and whether the correlations rationalise
the patterns in Part (ii).

---

## 3. Part (ii) — Staggered Difference-in-Differences

### 3.1 Frames

| Frame | Rule | Used by |
|---|---|---|
| `frame_did_broad`  | Full panel; $g_i = g_i^{\text{broad}}$; broadcast units count as treated | Main, robust, hetero, spillover |
| `frame_did_direct` | Drop broadcast-only units entirely; $g_i = g_i^{\text{dir}}$ | Sharp design |

Both frames carry the time-invariant baseline snapshot $W_i$ (see §3.4) and
the spatial-neighbour columns from `01_spatial_weights_ags8.py`.

**Two data-quality filters are applied in `00_prep_analysis.R` before the
frames are saved:**

1. **Drop `xbev < 500`.** KBA registers vehicles at the holder's HQ AGS8;
   *Großkunden-Halter* (leasing companies headquartered in tiny Gemeinden)
   produce extreme per-100k spikes (e.g. AGS8 01054108: population 304, 143
   BEVs registered in 2022 → 47,000 per 100k). Threshold of 500 inhabitants
   catches the worst artifacts and loses essentially no treated AGS8 (large
   Gemeinden are over-represented among EMK recipients anyway).
2. **Winsorize the per-100k outcomes at the 99th percentile** (`WINSOR_Q = 0.99`). Caps
   remaining fleet-registration artifacts at larger Gemeinden while
   preserving substantively important variation in the body of the
   distribution. Applied to all `bev_*_p100k` and `ice_neuzulassungen_p100k`.

For the CS package, $g_i$ is coded as 0 for never-treated. For BJS
(didimputation v0.5.1), $g_i$ is also coded as 0 for never-treated. Unit ID
is `ags8_id` (integer .GRP of AGS8).

### 3.2 Two co-primary estimators

#### 3.2.1 Callaway-Sant'Anna (CS), doubly robust

Group-time average treatment effects on the treated:

$$
\text{ATT}(g, t) = \mathbb{E}[Y_{it}(g) - Y_{it}(0) \mid g_i = g],
$$

estimated, for each $(g, t)$ with $g \le t$, as a 2×2 DiD between cohort $g$
and a control group ($\mathcal{N}$ never-treated, or $\mathcal{N}_t$
not-yet-treated). The doubly-robust estimator combines an outcome regression
on the control group and a propensity score $p(W) = \Pr(g_i = g \mid W_i)$:

$$
\widehat{\text{ATT}}(g, t)
  = \mathbb{E}\!\left[ \left(\frac{D_g}{\bar D_g} - \frac{p(W_i)(1 - D_g)}{(1 - p(W_i))\bar p}\right)
  \big(\Delta Y_{i,t,g} - \widehat{m}(W_i)\bigr)\right],
$$

where $D_g = \mathbf{1}\{g_i = g\}$, $\Delta Y_{i,t,g} = Y_{it} - Y_{i,g-1}$
(varying base period), and $\widehat{m}$ is the regression of $\Delta Y$ on
$W$ in the control group. The "doubly robust" property: consistent if EITHER
the propensity score OR the outcome regression is correctly specified.

#### 3.2.2 Borusyak-Jaravel-Spiess (BJS) imputation

Posit $Y_{it}(0) = \alpha_i + \lambda_t + u_{it}$ (TWFE in untreated
potential outcome). Estimate $\hat\alpha_i, \hat\lambda_t$ on the
*untreated* observations only (units that are not yet treated by year $t$,
plus all never-treated). Impute $\widehat{Y_{it}(0)} = \hat\alpha_i +
\hat\lambda_t$ for treated cells and compute

$$
\widehat{\tau}_{it} = Y_{it} - \widehat{Y_{it}(0)}.
$$

Aggregate to overall ATT or to event-time effects $\widehat{\tau}_e =
\mathbb{E}[\widehat\tau_{it} \mid t - g_i = e]$.

This estimator absorbs **time-invariant** unit heterogeneity into
$\alpha_i$, so the baseline covariates $W_i$ are absorbed automatically.
This is why CS-dr is the one that carries explicit baseline covariates —
the two estimators end up conditioning on the same time-invariant
heterogeneity through different mechanisms (FE vs `xformla`), making the
estimates comparable.

### 3.3 Design grid (in `03_did_main.R`)

| Role | Frame | CS control | Anticipation |
|---|---|---|---|
| **Main**       | broad  | never-treated  | 0 |
| **Sharp**      | direct | never-treated  | 0 |
| **Selection check** | broad  | not-yet-treated | 1 |

The selection check + anticipation = 1 jointly tests (a) whether using
future-treated units as controls changes the headline (Roth-Sant'Anna
selection concern) and (b) whether onset year is itself a noisy proxy for
the start of any policy attention.

### 3.4 Baseline covariates $W_i$ (CS-dr only)

For each AGS8, take the earliest year $\bar t_i \in \{2014, 2015, 2016\}$
with the variable observed:

| Component | Source variable | Transform |
|---|---|---|
| `kk_base_z`    | `log_kaufkraft` at $\bar t_i$         | cross-sectional $z$ |
| `sk_base_z`    | `log_steuerkraft` at $\bar t_i$       | $z$ |
| `dens_base_z`  | `log_pop_dens` at $\bar t_i$          | $z$ |
| `green_base_z` | `muni_gruene` at $\bar t_i$           | $z$ |
CS-dr formula: `xformla = ~ kk_base_z + sk_base_z + dens_base_z + green_base_z`.

Note: `bev_base_z` and `chg_base_z` are excluded — baseline 2014–2016
BEV/charging counts are mass-zero across the AGS8 cross-section and collapse
the DR design matrix to singular inside small $(g, t)$ cells.

These are **time-invariant** by construction. Rationale (in three words):
no-bad-controls; CS-has-no-unit-FE; absorbed-by-BJS-anyway. See the
`pipeline_summary.md` Q&A for the long version.

### 3.5 Outcomes

| Symbol | Column | Role |
|---|---|---|
| $Y^{\text{BEV}}_{it}$    | `bev_neuzulassungen_p100k` | headline |
| $Y^{\text{corp}}_{it}$   | `bev_corporate_p100k`      | holder-type decomposition |
| $Y^{\text{priv}}_{it}$   | `bev_private_p100k`        | holder-type decomposition |
| $Y^{\text{stock}}_{it}$  | `bev_stock_p100k`          | cumulative, secondary |
| $Y^{\text{ICE}}_{it}$    | `ice_neuzulassungen_p100k` | placebo (constructed from `N_benzin_diesel_overall` counts) |

Charging-point density is *not* estimated as an outcome. It enters the
analysis only as a covariate: as `chg_z` in the hazard channels (§2.2) and
as `chg_base_z` in the CS-dr baseline covariate vector (§3.4).

All outcomes are per-100k by construction (raw KBA counts divided by
year-specific population). No log transforms (per plan v2).

### 3.6 Aggregations

For CS:

$$
\text{ATT}^{\text{simple}}  = \frac{\sum_{(g,t): t \ge g} N_g \cdot \text{ATT}(g,t)}{\sum_{(g,t): t \ge g} N_g}
$$

$$
\text{ATT}^{\text{dynamic}}(e) = \mathbb{E}_g\bigl[\text{ATT}(g, g+e)\bigr], \quad e \in [-4, 4]
$$

Implemented as `aggte(type = "simple")` and `aggte(type = "dynamic",
min_e = -4, max_e = 4, balance_e = NULL, cband = TRUE)`. The dynamic
aggregation does **not** balance event time (composition changes across
$e$ are noted in captions).

For BJS, the event-time estimate $\widehat\tau_e$ is returned directly from
`did_imputation(horizon = TRUE, pretrends = TRUE)`. The static
ATT comes from `horizon = NULL` (label `"treat"` in the returned data
frame).

### 3.7 Inference

| Estimator | Variance | Cluster |
|---|---|---|
| CS-dr | Multiplier bootstrap, $B = 2000$, simultaneous bands via `cband = TRUE` | `clustervars = "AGS5"` |
| BJS   | Analytic, cluster-robust at the same level                                | `cluster_var = "AGS5"` |

AGS5 clustering matches the Kreis-broadcast structure: any single AGS5
project induces correlated treatment for all its AGS8, so within-Kreis
errors are dependent. Clustering at AGS8 would deflate SEs.

Bootstrap critical values for simultaneous bands ($\hat c$) are recovered
from the multiplier bootstrap percentiles; pointwise 95% bands use 1.96.

### 3.8 Pre-trend evidence

- **CS Wald test**: $W = \widehat{\theta}_{\text{pre}}'\,\widehat V_{\text{pre}}^{-1}\,\widehat{\theta}_{\text{pre}}$
  over pre-treatment $\text{ATT}(g,t)$ at $t < g$. Returned in
  `cs_obj$Wpval`. Written to `pretests.csv`.
- **BJS leads**: per-event-time pre-treatment estimates $\widehat\tau_{-4},
  \ldots, \widehat\tau_{-1}$ are displayed in the event-study plots
  (`es_<outcome>_<frame>_<ctrl>.pdf`). A joint pre-trend p-value is not
  exposed by didimputation v0.5.1.

### 3.9 Identification assumptions (Part ii)

For CS-dr with never-treated control and main / sharp roles:

1. **Conditional parallel trends in differences**: for each $g$,
   $\mathbb{E}[\Delta Y_{it}(0) \mid g_i = g, W_i] = \mathbb{E}[\Delta Y_{it}(0) \mid g_i = 0, W_i]$
   for $t \ge g$.
2. **No anticipation**: $Y_{it}(g) = Y_{it}(0)$ for $t < g$ (slightly
   weakened by `anticipation = 1` in the selection check).
3. **Overlap**: $0 < p(W_i) < 1$ within each $(g, t)$ cell; CS-dr's
   doubly-robust scoring tolerates moderate violations.
4. **SUTVA**: $Y_{it}(g_i)$ does not depend on $g_j, j \ne i$. The Kreis-
   broadcast mechanism creates within-Kreis dependence; the donut
   robustness drops never-treated controls with treated neighbours to probe
   this.

For BJS: TWFE structure on $Y_{it}(0)$ (separable unit + time effects on
the untreated branch).

### 3.10 Variant analyses (`04_did_robustness.R`)

| Variant | Sample modification | Anticipation | xformla | Control |
|---|---|---|---|---|
| `baseline`           | none                                              | 0 | NULL          | nevertreated |
| `conditional`        | none                                              | 0 | `XFORMLA_CS`  | nevertreated |
| `anticipation_1`     | none                                              | 1 | NULL          | nevertreated |
| `drop_2016`          | drop cohort 2016 (only one pre-period)            | 0 | NULL          | nevertreated |
| `drop_covid`         | drop $t \in \{2020, 2021\}$ (short-run only)      | 0 | NULL          | nevertreated |
| `complete_case`      | drop rows where `N_elektro_overall_imp == TRUE`   | 0 | NULL          | nevertreated |
| `notyet_evertreated` | ever-treated only                                 | 0 | NULL          | notyettreated |

Estimators per variant: BJS (static, `horizon = NULL`) and CS (simple ATT aggregation).

---

## 4. Heterogeneity (`05_heterogeneity.R`)

Sample-split CS-dr (broad frame, never-treated control) within unweighted
terciles of $W_i$. For each baseline rank $r \in \{\text{kk\_base}, \text{sk\_base}\}$
and each tercile $\tau \in \{1, 2, 3\}$, estimate

$$
\text{ATT}^{r,\tau} = \frac{\sum N_g\, \text{ATT}(g, t \mid \tau_i = \tau)}{\sum N_g}
$$

both as simple aggregation and dynamic event study. Difference of top vs
bottom tercile $\Delta = \text{ATT}^{r,3} - \text{ATT}^{r,1}$ is reported
with a normal-approximation CI (independence assumption across split-sample
bootstraps).

Treat-type heterogeneity: restrict the broad frame to
$\{$direct, never$\}$ and to $\{$broadcast-only, never$\}$ and compare the
two CS-dr ATTs.

---

## 5. Spillovers (`06_spillovers.R`)

Spatial structure from `01_spatial_weights_ags8.py`:

- **Granular**: AGS8 queen contiguity, 1st/2nd/3rd ring; treated-neighbour
  indicators use $g_i^{\text{dir}}$.
- **Aggregated**: AGS5 queen contiguity (after polygon dissolve);
  treated-neighbour indicator uses $g_i^{\text{broad}}$.

### 5.1 Donut robustness

Re-estimate the main CS-dr spec after dropping never-treated controls $i$
with $\text{direct\_treated\_any\_nbrs\_gem\_1}_{it} = 1$ OR
$\text{broad\_treated\_any\_nbrs\_kreis}_{it} = 1$. If the headline ATT
survives, contamination of "never-treated" controls by nearby treatment is
unlikely to drive the result.

### 5.2 Descriptive spillover event study

Restrict to never-treated AGS8 only. Pseudo-cohort
$\tilde g_i = \min\{t : \text{direct\_treated\_any\_nbrs\_gem\_1}_{it} = 1\}$.
BJS event study on the BEV outcome with this pseudo-treatment. Plotted as
$\widehat\tau_e$ in $e \in [-4, 4]$.

Caveat: the pseudo-cohort is endogenous to spatial-economic geography. The
event study is descriptive only — selection-on-geography rules out a
clean causal reading.

---

## 6. Standardisation policy (current implementation)

Hazard channels $X_{it}$ are $z$-scored on the estimation sample
(`zfit(ph)`), and re-scored on the no-Stadtstaaten subsample
(`zfit(ph_ns)`). Trade-off:

- Within-sample $z$: $\beta_k$ reads as "effect of 1 SD on this sample."
- Cross-spec $z$ comparability: coefficients across columns are NOT
  directly comparable in absolute magnitude because the SDs differ between
  (1) and (2).

Baseline covariates $W_i$ are $z$-scored on the AGS8 cross-section once
(in `00_prep_analysis.R`) and used unchanged across all CS-dr cells.

---

## 7. Sample sizes (current run; will refresh when KBA reruns)

| Sample | $N$ obs | $N$ AGS8 | Events |
|---|---|---|---|
| `frame_hazard` (year $\ge$ 2015, complete L1)      | 85,302 | 9,607 | 166 direct onsets |
| `frame_hazard` no-Stadtstaaten                     | 85,275 | 9,603 | 163 |
| `frame_did_broad`                                  | full panel | 10,775 | ever-treated: 1,827 |
| `frame_did_direct`                                 | drops broadcast-only | $\approx$ 8,948 + 167 | 167 |
| Cohorts (broad)                                    | 2015–2022 (varying counts) | — | — |

---

## 8. Open / parked methodological items

| Item | Status |
|---|---|
| Joint BJS pretrend p-value | Not exposed in didimputation v0.5.1; rely on CS Wald + BJS lead plots |
| `MAX_COHORT` policy-regime cap | NA for first pass; toggle if restricting to pre-2021 FRL regime |
| KBA boundary harmonisation to 2023 (Gebietsstand) | Awaiting BBSR annual ref-gemeinden files 2013–2020 |
| `modellregion_pre2015` Kreis dummy | Placeholder column; will enter selection equations once supplied |
| Switch hazard to common $z$-reference (full hazard frame) | Discussed; not yet implemented |
| Tercile difference: bootstrap vs normal approx | Currently normal approx; full split-sample bootstrap optional |
