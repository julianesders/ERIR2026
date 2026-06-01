import pandas as pd

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

# Charging stations at Gemeinde (AGS8) level
panel = panel.merge(ladestationen, on=["AGS8", "year"], how="left")
panel["ev_stations"]     = panel["ev_stations"].fillna(0).astype(int)
panel["ev_chargepoints"] = panel["ev_chargepoints"].fillna(0).astype(int)

# Normalised EV rates per 100,000 inhabitants (xbev is Kreis-level, broadcast)
panel["ev_stations_p100k"]     = panel["ev_stations"]     / panel["xbev"] * 100_000
panel["ev_chargepoints_p100k"] = panel["ev_chargepoints"] / panel["xbev"] * 100_000

# Elections at Gemeinde level
panel = panel.merge(elections, on=["AGS8", "year"], how="left")

# Personnel (AGS5 level, broadcast to all Gemeinden in Kreis)
panel = panel.merge(personal, on=["AGS5", "year"], how="left")

# Fill activity indicators with 0 for never-treated units
for col in ["emk_active", "n_emk_active", "n_emk_total", "emk_absorbing", "emk_absorbing_n"]:
    panel[col] = panel[col].fillna(0).astype(int)


# ── Column order ──────────────────────────────────────────────────────────────

id_cols        = ["AGS8", "AGS5", "year"]
activity_cols  = ["emk_active", "n_emk_active", "n_emk_total", "emk_absorbing", "emk_absorbing_n"]
ladestation_cols = [
    "ev_stations", "ev_stations_p100k",
    "ev_chargepoints", "ev_chargepoints_p100k",
]
election_cols  = [c for c in elections.columns if c not in ["AGS8", "year"]]
inkar_cols     = [c for c in inkar.columns if c not in ("AGS8", "AGS5", "year")]
personal_cols  = ["n_vze_personal"]
emk_attr_cols  = [c for c in emk_attrs.columns if c not in ["AGS5"] + activity_cols]

panel = panel[id_cols + activity_cols + ladestation_cols + election_cols + inkar_cols + personal_cols + emk_attr_cols]

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
