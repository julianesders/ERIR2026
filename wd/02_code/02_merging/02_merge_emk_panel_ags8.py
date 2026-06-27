import numpy as np
import pandas as pd

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parents[2]))
from config import BASE_PATH
DATA_DIR_INTERMEDIATE = f"{BASE_PATH}/01_data/02_intermediate"
DATA_DIR_FINAL        = f"{BASE_PATH}/01_data/03_final"

# ── Load inputs ───────────────────────────────────────────────────────────────

# Spine: Gemeinden × all available years with INKAR variables
inkar = pd.read_csv(
    f"{DATA_DIR_INTERMEDIATE}/inkar/inkar_ags8_panel.csv",
    dtype={"AGS8": str, "AGS5": str},
)

personal = pd.read_csv(
    f"{DATA_DIR_INTERMEDIATE}/inkar/personal_ags5_panel.csv",
    dtype={"AGS5": str},
)
personal["AGS5"] = personal["AGS5"].str.zfill(5)

# Convert raw VZE headcount to per-100k using AGS5 population (sum of AGS8 xbev).
# Done here so that personal stays a clean raw-count file.
_pop_ags5 = (
    inkar.groupby(["AGS5", "year"], as_index=False)["xbev"]
    .sum()
    .rename(columns={"xbev": "_pop_ags5"})
)
personal = personal.merge(_pop_ags5, on=["AGS5", "year"], how="left")
personal["n_vze_personal"] = personal["n_vze_personal"] / personal["_pop_ags5"] * 100_000
personal = personal.drop(columns=["_pop_ags5"])

ladestationen = pd.read_csv(
    f"{DATA_DIR_INTERMEDIATE}/ladestationen/ladestationen_ags8_panel.csv",
    dtype={"AGS8": str},
)
ladestationen["AGS8"] = ladestationen["AGS8"].str.zfill(8)

elections = pd.read_csv(
    f"{DATA_DIR_INTERMEDIATE}/elections/elections_ags8_panel.csv",
    dtype={"AGS8": str},
)
elections["AGS8"] = elections["AGS8"].str.zfill(8)

kba = pd.read_csv(
    f"{DATA_DIR_INTERMEDIATE}/kba/kba_ags8_panel.csv",
    dtype={"AGS8": str},
    usecols=[
        "AGS8", "year",
        "B_elektro_overall",
        "N_elektro_overall", "N_elektro_private", "N_elektro_corporate",
        "N_benzin_diesel_overall",
        "N_ev_share_overall", "N_ev_share_private", "N_ev_share_corporate",
    ],
)
kba["AGS8"] = kba["AGS8"].str.zfill(8)

# Inf can arise from share denominators of 0 in the KBA aggregate; treat as NA.
for _c in ["N_ev_share_overall", "N_ev_share_private", "N_ev_share_corporate"]:
    kba[_c] = kba[_c].replace([np.inf, -np.inf], np.nan)

area = pd.read_csv(
    f"{DATA_DIR_INTERMEDIATE}/area/area_ags8_panel.csv",
    dtype={"AGS8": str},
)
area["AGS8"] = area["AGS8"].str.zfill(8)

emk = pd.read_csv(
    f"{DATA_DIR_INTERMEDIATE}/emk/emk_ags_matched.csv",
    dtype={"AGS8": str, "AGS5": str},
)
emk["AGS5"] = emk["AGS5"].str.zfill(5)
emk["AGS8"] = emk["AGS8"].where(emk["AGS8"].notna(), other=None)
emk.loc[emk["AGS8"].notna(), "AGS8"] = emk.loc[emk["AGS8"].notna(), "AGS8"].str.zfill(8)
emk["laufzeit_start"] = pd.to_datetime(emk["laufzeit_start"])
emk["laufzeit_end"]   = pd.to_datetime(emk["laufzeit_end"])
emk["start_year"]     = emk["laufzeit_start"].dt.year
emk["end_year"]       = emk["laufzeit_end"].dt.year


########################################
########   Treatment indicators  #######
########################################


# 187 projects matched directly to AGS8; 54 are Kreis-level (AGS8 = NaN).
# Kreis-level projects are broadcast to every Gemeinde in that Kreis so
# that treatment is also defined at the Gemeinde level. We track DIRECT and
# BROAD (direct ∪ Kreis-broadcast) separately throughout: the hazard risk set
# uses DIRECT events, the DiD main spec uses BROAD treatment.

ags_map = inkar[["AGS8", "AGS5"]].drop_duplicates()

emk_direct    = emk[emk["AGS8"].notna()][["AGS8", "AGS5", "start_year", "end_year"]]
emk_kreis     = emk[emk["AGS8"].isna()][["AGS5", "start_year", "end_year"]]
emk_kreis_exp = emk_kreis.merge(ags_map, on="AGS5", how="left")[
    ["AGS8", "AGS5", "start_year", "end_year"]
]
emk_all = pd.concat([emk_direct, emk_kreis_exp], ignore_index=True)

print(f"Direct AGS8 projects:            {len(emk_direct)} "
      f"({emk_direct['AGS8'].nunique()} unique Gemeinden)")
print(f"Kreis-level projects:            {len(emk_kreis)} "
      f"→ broadcast to {emk_kreis_exp['AGS8'].nunique()} Gemeinden")

panel_years = inkar[["AGS8", "year"]].drop_duplicates()

# Active / started indicators, broad coverage (direct ∪ Kreis broadcast)
activity = (
    panel_years
    .merge(emk_all[["AGS8", "start_year", "end_year"]], on="AGS8", how="left")
    .assign(
        active  =lambda df: (df["start_year"] <= df["year"]) & (df["year"] <= df["end_year"]),
        started =lambda df:  df["start_year"] <= df["year"],
    )
    .groupby(["AGS8", "year"], as_index=False)
    .agg(
        emk_active      =("active",  "any"),
        n_emk_active    =("active",  "sum"),
        emk_absorbing   =("started", "any"),
        emk_absorbing_n =("started", "sum"),
    )
)
for col in ["emk_active", "n_emk_active", "emk_absorbing", "emk_absorbing_n"]:
    activity[col] = activity[col].astype(int)

# Per-AGS8 first-treatment years for hazard / DiD frames
first_direct = (
    emk_direct.groupby("AGS8", as_index=False)["start_year"].min()
              .rename(columns={"start_year": "first_treat_direct"})
)
first_broad = (
    emk_all.groupby("AGS8", as_index=False)["start_year"].min()
           .rename(columns={"start_year": "first_treat_broad"})
)

# Kreis-funded year (AGS5-level): independent of own direct status.
# Used as a strict-past covariate in the hazard: 1{kreis_funded_year < t}.
first_kreis = (
    emk_kreis.groupby("AGS5", as_index=False)["start_year"].min()
             .rename(columns={"start_year": "kreis_funded_year"})
)

treat_map = (
    ags_map
    .merge(first_direct, on="AGS8", how="left")
    .merge(first_broad,  on="AGS8", how="left")
    .merge(first_kreis,  on="AGS5", how="left")
)
treat_map["treat_type"] = np.where(
    treat_map["first_treat_direct"].notna(), "direct",
    np.where(treat_map["first_treat_broad"].notna(), "broadcast_only", "never"),
)

print("Treat-type counts (AGS8):")
print(treat_map["treat_type"].value_counts().to_string())


########################################
########   Build final panel    ########
########################################


panel = inkar.merge(activity, on=["AGS8", "year"], how="left")
panel = panel.merge(treat_map[
    ["AGS8", "first_treat_direct", "first_treat_broad", "treat_type", "kreis_funded_year"]
], on="AGS8", how="left")
panel["treat_type"] = panel["treat_type"].fillna("never")
panel["AGS2"] = panel["AGS8"].str[:2]

# TODO: merge `modellregion_pre2015` (AGS5-level dummy) once supplied by user;
# join on AGS5. Placeholder so downstream scripts can reference the column.
panel["modellregion_pre2015"] = 0

# Drop rows where population is missing — xbev is the denominator for all
# per-capita variables; rows without it are uninformative.
_n_before = len(panel)
panel = panel[panel["xbev"].notna() & (panel["xbev"] > 0)].copy()
print(f"Dropped {_n_before - len(panel)} rows with missing or zero xbev "
      f"({panel['AGS8'].nunique()} AGS8 units remaining)")

# ── Gemeinde area (AGS8 × year) ───────────────────────────────────────────────

panel = panel.merge(area, on=["AGS8", "year"], how="left")

# ── Charging stations (AGS8 level) ────────────────────────────────────────────

panel = panel.merge(ladestationen, on=["AGS8", "year"], how="left")
panel["ev_stations"]     = panel["ev_stations"].fillna(0).astype(int)
panel["ev_chargepoints"] = panel["ev_chargepoints"].fillna(0).astype(int)
panel["ev_stations_p100k"]     = panel["ev_stations"]     / panel["xbev"] * 100_000
panel["ev_chargepoints_p100k"] = panel["ev_chargepoints"] / panel["xbev"] * 100_000

# ── BEV stock from KBA Bestand (AGS8 level) ───────────────────────────────────

panel = panel.merge(kba, on=["AGS8", "year"], how="left")

# Interpolate raw KBA counts before per-capita division so the denominator
# (population) stays year-specific. We record _imp flags so robustness can
# drop filled cells (complete-case).
#
# Interpolation policy (plan v2):
#   - N_* and N_ev_share_*  : interior-only linear, limit=2 (limit_area="inside")
#   - B_elektro_overall     : limit_direction="both", limit=2 (stock)

panel = panel.sort_values(["AGS8", "year"])

_n_interior_cols = [
    "N_elektro_overall", "N_elektro_private", "N_elektro_corporate",
    "N_benzin_diesel_overall",
    "N_ev_share_overall", "N_ev_share_private", "N_ev_share_corporate",
]
for _col in _n_interior_cols:
    panel[f"{_col}_imp"] = panel[_col].isna()
    panel[_col] = panel.groupby("AGS8")[_col].transform(
        lambda s: s.interpolate(method="linear", limit=2, limit_area="inside")
    )
    # An _imp cell is only "filled" if it became non-NaN after interpolation.
    panel[f"{_col}_imp"] = panel[f"{_col}_imp"] & panel[_col].notna()

panel["B_elektro_overall_imp"] = panel["B_elektro_overall"].isna()
panel["B_elektro_overall"] = panel.groupby("AGS8")["B_elektro_overall"].transform(
    lambda s: s.interpolate(method="linear", limit=2, limit_direction="both")
)
panel["B_elektro_overall_imp"] = panel["B_elektro_overall_imp"] & panel["B_elektro_overall"].notna()

panel["bev_stock_p100k"]           = panel["B_elektro_overall"]       / panel["xbev"] * 100_000
panel["bev_neuzulassungen_p100k"]  = panel["N_elektro_overall"]       / panel["xbev"] * 100_000
panel["bev_corporate_p100k"]       = panel["N_elektro_corporate"]     / panel["xbev"] * 100_000
panel["bev_private_p100k"]         = panel["N_elektro_private"]       / panel["xbev"] * 100_000
# ICE placebo from counts (plan v2): avoids share-inversion at near-zero shares.
panel["ice_neuzulassungen_p100k"]  = panel["N_benzin_diesel_overall"] / panel["xbev"] * 100_000
panel = panel.drop(columns=["B_elektro_overall"])

_bev_holes = panel["bev_stock_p100k"].isna().sum()
print(f"bev_stock_p100k NaN after interpolation: {_bev_holes} "
      f"({100 * _bev_holes / len(panel):.1f}% of rows — gap > 2 yrs or pre-KBA)")

# ── Derived INKAR variables ───────────────────────────────────────────────────

# Population density (log): population / land area in km²
panel["log_pop_dens"] = np.log(panel["xbev"] / panel["area_qkm"])

# Log Steuerkraft: clip to 1 (handles negatives and exact zeros; no values in (0,1)
# in the estimation sample), then conventional log. Gives clean elasticity reading.
panel["log_steuerkraft"] = np.log(panel["q_gest_bev"].clip(lower=1))

# ── Elections (AGS8 level) ────────────────────────────────────────────────────

panel = panel.merge(elections, on=["AGS8", "year"], how="left")

# ── Personnel (AGS5 level, broadcast to Gemeinden) ───────────────────────────

panel = panel.merge(personal, on=["AGS5", "year"], how="left")

# Fill activity indicators with 0 for never-treated units
for col in ["emk_active", "n_emk_active", "emk_absorbing", "emk_absorbing_n"]:
    panel[col] = panel[col].fillna(0).astype(int)

# ── Lag variables (within AGS8, year-based) ───────────────────────────────────
# Using year-based merge: L1 = value at year t-1, etc. Handles gaps correctly.

LAG_VARS = [
    "q_gest_bev",      # steuerkraft (raw, kept for reference)
    "log_steuerkraft",
    "bev_stock_p100k",
    "ev_chargepoints_p100k",
    "fed_gruene",
    "state_gruene",
    "muni_gruene",
    "n_vze_personal",
    "N_elektro_overall",
    "N_elektro_private",
    "N_elektro_corporate",
    "N_benzin_diesel_overall",
    "N_ev_share_overall",
    "N_ev_share_private",
    "N_ev_share_corporate",
]

panel = panel.sort_values(["AGS8", "year"]).reset_index(drop=True)

for k in [1, 2, 3]:
    lag_src = panel[["AGS8", "year"] + LAG_VARS].copy()
    lag_src["year"] = lag_src["year"] + k
    lag_src = lag_src.rename(columns={c: f"{c}_L{k}" for c in LAG_VARS})
    panel = panel.merge(lag_src, on=["AGS8", "year"], how="left")

# ── Column order ──────────────────────────────────────────────────────────────

lag_cols_all  = [f"{c}_L{k}" for k in [1, 2, 3] for c in LAG_VARS]
derived_cols  = [
    "area_qkm", "log_pop_dens", "log_steuerkraft",
    "bev_stock_p100k", "bev_neuzulassungen_p100k",
    "bev_corporate_p100k", "bev_private_p100k",
    "ice_neuzulassungen_p100k",
    "N_ev_share_overall", "N_ev_share_private", "N_ev_share_corporate",
]
id_cols       = ["AGS8", "AGS5", "AGS2", "year"]
treat_cols    = ["first_treat_direct", "first_treat_broad", "treat_type",
                 "kreis_funded_year", "modellregion_pre2015"]
activity_cols = ["emk_active", "n_emk_active", "emk_absorbing", "emk_absorbing_n"]
ladestation_cols = [
    "ev_stations", "ev_stations_p100k",
    "ev_chargepoints", "ev_chargepoints_p100k",
]
election_cols = [c for c in elections.columns if c not in ["AGS8", "year"]]
inkar_cols    = [c for c in inkar.columns if c not in ("AGS8", "AGS5", "year")]
personal_cols = ["n_vze_personal"]
imp_cols      = [c for c in panel.columns if c.endswith("_imp")]

panel = panel[
    id_cols + treat_cols + activity_cols + ladestation_cols + election_cols +
    inkar_cols + derived_cols + personal_cols + imp_cols + lag_cols_all
]

panel.to_csv(f"{DATA_DIR_FINAL}/emk_inkar_panel_ags8.csv", index=False)

print(f"\nPanel shape:               {panel.shape}")
print(f"AGS8 units:                {panel['AGS8'].nunique()}")
print(f"AGS5 (Kreis) units:        {panel['AGS5'].nunique()}")
print(f"Years:                     {sorted(panel['year'].unique().tolist())}")
print(f"Ever-treated AGS8 (broad): {(panel.groupby('AGS8')['emk_active'].max() == 1).sum()}")
print(f"  direct project match:    {emk_direct['AGS8'].nunique()}")
print(f"  Kreis broadcast only:    {emk_kreis_exp[~emk_kreis_exp['AGS8'].isin(emk_direct['AGS8'])]['AGS8'].nunique()}")
print(f"Never-treated AGS8:        {(panel.groupby('AGS8')['emk_active'].max() == 0).sum()}")
print(f"Max concurrent projects:   {panel['n_emk_active'].max()}")
print(f"\nNew derived columns:       {derived_cols}")
print(f"Imputation flags ({len(imp_cols)}): {imp_cols}")
print(f"Lag columns ({len(lag_cols_all)}): {lag_cols_all[:6]} ...")
