import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

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
        "N_ev_share_overall", "N_ev_share_private", "N_ev_share_corporate",
    ],
)
kba["AGS8"] = kba["AGS8"].str.zfill(8)

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
########   Activity indicators  ########
########################################


# 187 projects matched directly to AGS8; 54 are Kreis-level (AGS8 = NaN).
# Kreis-level projects are broadcast to every Gemeinde in that Kreis so
# that treatment is defined at the Gemeinde level throughout.

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


########################################
########  Time-invariant EMK attrs  ####
########################################


tag_cols   = [c for c in emk.columns if c.startswith("tag_")]
space_cols = [c for c in emk.columns if c.startswith("space_")]
for c in tag_cols + space_cols:
    emk[c] = pd.to_numeric(emk[c], errors="coerce")
emk["gesamtmittel"] = pd.to_numeric(emk["gesamtmittel"], errors="coerce")
emk["bundesmittel"] = pd.to_numeric(emk["bundesmittel"], errors="coerce")

# Aggregate at Kreis level and broadcast to all Gemeinden via AGS5 join
emk_attrs = emk.groupby("AGS5", as_index=False).agg(**{
    "emk_gesamtmittel": ("gesamtmittel", "sum"),
    "emk_bundesmittel": ("bundesmittel", "sum"),
    "n_emk_total":      ("gesamtmittel", "count"),
    **{c: (c, "max") for c in tag_cols + space_cols},
})


########################################
########   Build final panel    ########
########################################


panel = inkar.merge(activity, on=["AGS8", "year"], how="left")
panel = panel.merge(emk_attrs, on="AGS5", how="left")
panel["AGS2"] = panel["AGS8"].str[:2]

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
# (population) stays year-specific. limit=2: fill at most 2 consecutive missing
# years; limit_direction="both" covers interior gaps (linear), leading NAs
# (backward constant from first value), and trailing NAs (forward constant).
# Must happen before eco_index PCA so the fill propagates into the index and
# subsequently into all _L1 lags.
_kba_count_cols = [
    "B_elektro_overall",
    "N_elektro_overall", "N_elektro_private", "N_elektro_corporate",
    "N_ev_share_overall", "N_ev_share_private", "N_ev_share_corporate",
]
panel = panel.sort_values(["AGS8", "year"])
for _col in _kba_count_cols:
    panel[_col] = panel.groupby("AGS8")[_col].transform(
        lambda s: s.interpolate(method="linear", limit=2, limit_direction="both")
    )

panel["bev_stock_p100k"]          = panel["B_elektro_overall"]    / panel["xbev"] * 100_000
panel["bev_neuzulassungen_p100k"] = panel["N_elektro_overall"]   / panel["xbev"] * 100_000
panel["bev_corporate_p100k"]      = panel["N_elektro_corporate"] / panel["xbev"] * 100_000
panel["bev_private_p100k"]        = panel["N_elektro_private"]   / panel["xbev"] * 100_000
panel = panel.drop(columns=["B_elektro_overall"])

_bev_holes = panel["bev_stock_p100k"].isna().sum()
print(f"bev_stock_p100k NaN after interpolation: {_bev_holes} "
      f"({100 * _bev_holes / len(panel):.1f}% of rows — gap > 2 yrs or pre-KBA)")

# ── Derived INKAR variables ───────────────────────────────────────────────────

# Population density (log): population / land area in km²
panel["log_pop_dens"] = np.log(panel["xbev"] / panel["area_qkm"])

# Log Steuerkraft: negatives (rare, fiscal equalization) clipped to 0 before log1p.
panel["log_steuerkraft"] = np.log1p(panel["q_gest_bev"].clip(lower=0))

# ── EV ecosystem index: first PC of (bev_stock_p100k, ev_chargepoints_p100k) ─

eco_vars = ["bev_stock_p100k", "ev_chargepoints_p100k"]
eco_mask = panel[eco_vars].notna().all(axis=1)
eco_data = np.log1p(panel.loc[eco_mask, eco_vars].values)

scaler    = StandardScaler()
pca       = PCA(n_components=1)
eco_scores = pca.fit_transform(scaler.fit_transform(eco_data)).ravel()

panel["eco_index"] = np.nan
panel.loc[eco_mask, "eco_index"] = eco_scores

print(f"\nEco index — PCA explained variance ratio: {pca.explained_variance_ratio_[0]:.1%}")
print(f"  loadings (bev_stock_p100k, ev_chargepoints_p100k): "
      f"{pca.components_[0].round(4).tolist()}")

# ── Elections (AGS8 level) ────────────────────────────────────────────────────

panel = panel.merge(elections, on=["AGS8", "year"], how="left")

# ── Personnel (AGS5 level, broadcast to Gemeinden) ───────────────────────────

panel = panel.merge(personal, on=["AGS5", "year"], how="left")

# Fill activity indicators with 0 for never-treated units
for col in ["emk_active", "n_emk_active", "n_emk_total", "emk_absorbing", "emk_absorbing_n"]:
    panel[col] = panel[col].fillna(0).astype(int)

# ── Lag variables (within AGS8, year-based) ───────────────────────────────────
# Using year-based merge: L1 = value at year t-1, etc. Handles gaps correctly.

LAG_VARS = [
    "q_gest_bev",      # steuerkraft (raw, kept for reference)
    "log_steuerkraft",
    "eco_index",
    "bev_stock_p100k",
    "ev_chargepoints_p100k",
    "fed_gruene",
    "state_gruene",
    "muni_gruene",
    "n_vze_personal",
    "N_elektro_overall",
    "N_elektro_private",
    "N_elektro_corporate",
    "N_ev_share_overall",
    "N_ev_share_private",
    "N_ev_share_corporate"
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
    "N_ev_share_overall", "N_ev_share_private", "N_ev_share_corporate",
    "eco_index",
]
id_cols       = ["AGS8", "AGS5", "AGS2", "year"]
activity_cols = ["emk_active", "n_emk_active", "n_emk_total", "emk_absorbing", "emk_absorbing_n"]
ladestation_cols = [
    "ev_stations", "ev_stations_p100k",
    "ev_chargepoints", "ev_chargepoints_p100k",
]
election_cols = [c for c in elections.columns if c not in ["AGS8", "year"]]
inkar_cols    = [c for c in inkar.columns if c not in ("AGS8", "AGS5", "year")]
personal_cols = ["n_vze_personal"]
emk_attr_cols = [c for c in emk_attrs.columns if c not in ["AGS5"] + activity_cols]

panel = panel[
    id_cols + activity_cols + ladestation_cols + election_cols +
    inkar_cols + derived_cols + personal_cols + emk_attr_cols + lag_cols_all
]

panel.to_csv(f"{DATA_DIR_FINAL}/emk_inkar_panel_ags8.csv", index=False)

print(f"\nPanel shape:               {panel.shape}")
print(f"AGS8 units:                {panel['AGS8'].nunique()}")
print(f"AGS5 (Kreis) units:        {panel['AGS5'].nunique()}")
print(f"Years:                     {sorted(panel['year'].unique().tolist())}")
print(f"Ever-treated AGS8:         {(panel.groupby('AGS8')['emk_active'].max() == 1).sum()}")
print(f"  direct project match:    {emk_direct['AGS8'].nunique()}")
print(f"  Kreis broadcast only:    {emk_kreis_exp[~emk_kreis_exp['AGS8'].isin(emk_direct['AGS8'])]['AGS8'].nunique()}")
print(f"Never-treated AGS8:        {(panel.groupby('AGS8')['emk_active'].max() == 0).sum()}")
print(f"Max concurrent projects:   {panel['n_emk_active'].max()}")
print(f"\nNew derived columns: {derived_cols}")
print(f"Lag columns ({len(lag_cols_all)}): {lag_cols_all[:6]} ...")
