import pandas as pd

BASE_PATH             = "/Users/julian/Documents/ERIR2026/wd"
DATA_DIR_RAW          = f"{BASE_PATH}/01_data/01_raw"
DATA_DIR_INTERMEDIATE = f"{BASE_PATH}/01_data/02_intermediate"

CROSSWALK_PATH = f"{DATA_DIR_RAW}/crosswalks/ref-kreise-1990-2024.xlsx"

# ── Load ──────────────────────────────────────────────────────────────────────
# File: Statistik der Beschäftigten der Gemeinden und Gemeindeverbände
# Source: Statistische Ämter des Bundes und der Länder (table 74111-04-04-4-B)
# 11 header rows before data; separator is semicolon; encoding latin-1
#
# Column layout (0-indexed after skiprows):
#   0: date (dd.mm.yyyy, stichtag 30.06.)
#   1: AGS code — mostly 5-digit Kreise, but a small number of entities
#      (Berlin Bezirke, Hannover, Aachen, Saarbrücken) are reported at
#      AGS8 level; those are excluded below.
#   2: Raumeinheit name  (dropped)
#   3: dimension: männlich / weiblich / Insgesamt
#   4: Vollzeitäquivalent der Beschäftigten — Insgesamt
#   5-10: subcategory breakdowns  (dropped)

raw = pd.read_csv(
    f"{DATA_DIR_RAW}/erir26_personal.csv",
    sep=";",
    header=None,
    skiprows=11,
    encoding="latin-1",
    usecols=[0, 1, 3, 4],
    names=["date", "AGS5", "dimension", "n_vze_personal"],
    dtype=str,
    on_bad_lines="skip",
)

# ── Filter ────────────────────────────────────────────────────────────────────

# Keep only "Insgesamt" rows (total across genders)
raw = raw[raw["dimension"] == "Insgesamt"].copy()

# Drop AGS8-level entries (Berlin Bezirke, Hannover, Aachen, Saarbrücken)
raw = raw[raw["AGS5"].str.strip().str.len() == 5].copy()

# Extract year from date string "dd.mm.yyyy"
raw["year"] = pd.to_numeric(raw["date"].str[-4:], errors="coerce")
raw = raw[raw["year"] >= 2005].copy()

# Zero-pad AGS5 to 5 digits
raw["AGS5"] = raw["AGS5"].str.strip().str.zfill(5)

# Parse value: "-" marks statistically suppressed cells
raw["n_vze_personal"] = (
    raw["n_vze_personal"]
    .str.strip()
    .replace("-", pd.NA)
)
raw["n_vze_personal"] = pd.to_numeric(raw["n_vze_personal"], errors="coerce")

# ── Build Kreis reform mapping (2010–2023) ────────────────────────────────────
# The 74111 source reports data as of 30.06. each year, using the AGS5 codes
# active on that date.  Three Kreisreformen fall within the panel period:
#
#   - MV 2011  (September 4):  13001–13062 → 13071–13076  [sheet 2010-2011]
#   - NI 2016  (November 1):   03152+03156 → 03159         [sheet 2015-2016]
#   - Eisenach (July 1, 2021): 16056 → 16063               [sheet 2020-2021]
#
# Since all reforms occurred after June 30 of the transition year, data for
# year Y uses old codes when Y ≤ year_cutoff (= the "to-year" of the sheet).
# Population-proportional weights handle the one split case (13052 → two successors).

def _build_kreis_reform_mapping(path: str) -> pd.DataFrame:
    records = []
    for y in range(2005, 2024):
        sheet = f"{y}-{y+1}"
        df = pd.read_excel(path, sheet_name=sheet, dtype=str)
        w_col    = [c for c in df.columns if "bevölkerungs" in c.lower()][0]
        ags_from = df.iloc[:, 0].str.strip().str.zfill(8).str[:5]
        ags_to   = df.iloc[:, -2].str.strip().str.zfill(8).str[:5]
        weight   = pd.to_numeric(df[w_col], errors="coerce").fillna(1.0)

        changed = ags_from != ags_to
        if changed.any():
            records.append(pd.DataFrame({
                "ags5_from":   ags_from[changed].values,
                "ags5_to":     ags_to[changed].values,
                "weight":      weight[changed].values,
                "year_cutoff": y + 1,
            }))

    if not records:
        return pd.DataFrame(columns=["ags5_from", "ags5_to", "weight", "year_cutoff"])
    return pd.concat(records, ignore_index=True)

reform_map = _build_kreis_reform_mapping(CROSSWALK_PATH)
print(f"Kreis reform mapping: {len(reform_map)} entries")
print(reform_map.to_string(index=False))

# ── Apply crosswalk harmonization ─────────────────────────────────────────────
# Merge reform mapping; for old codes active in a given year, remap to their
# 2023 successor and apply the population weight. Unmatched codes stay as-is.

raw_m = raw.merge(
    reform_map.rename(columns={"ags5_from": "AGS5"}),
    on="AGS5",
    how="left",
)

applies = raw_m["year"] <= raw_m["year_cutoff"]
raw_m["AGS5_2023"]  = raw_m["AGS5"]
raw_m["eff_weight"] = 1.0
raw_m.loc[applies, "AGS5_2023"]  = raw_m.loc[applies, "ags5_to"]
raw_m.loc[applies, "eff_weight"] = raw_m.loc[applies, "weight"]

raw_m["n_vze_personal"] = raw_m["n_vze_personal"] * raw_m["eff_weight"]

panel = (
    raw_m.groupby(["AGS5_2023", "year"], as_index=False)["n_vze_personal"]
    .sum(min_count=1)
    .rename(columns={"AGS5_2023": "AGS5"})
    .assign(year=lambda df: df["year"].astype(int))
    .sort_values(["AGS5", "year"])
    .reset_index(drop=True)
)

# ── Export ────────────────────────────────────────────────────────────────────

panel.to_csv(f"{DATA_DIR_INTERMEDIATE}/inkar/personal_ags5_panel.csv", index=False)

print(f"\nPanel shape:            {panel.shape}")
print(f"AGS5 units:             {panel['AGS5'].nunique()}")
print(f"Years:                  {sorted(panel['year'].unique().tolist())}")
print(f"NaN values:             {panel['n_vze_personal'].isna().sum()}")
print(f"\nn_vze_personal summary:")
print(panel["n_vze_personal"].describe().to_string())

print(f"\nEntities with ANY NaN (Stadtstaaten — data suppressed at source):")
nan_ags5 = panel[panel["n_vze_personal"].isna()]["AGS5"].unique()
for ags in sorted(nan_ags5):
    nan_years = sorted(
        panel[(panel["AGS5"] == ags) & panel["n_vze_personal"].isna()]["year"].tolist()
    )
    print(f"  {ags}  NaN in years: {nan_years}")
