import os
import pandas as pd

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parents[2]))
from config import BASE_PATH
DATA_DIR_RAW          = f"{BASE_PATH}/01_data/01_raw"
DATA_DIR_INTERMEDIATE = f"{BASE_PATH}/01_data/02_intermediate"

# ── Load ──────────────────────────────────────────────────────────────────────
# Source: Statistische Ämter des Bundes und der Länder, "Gebietsfläche in qkm",
# Stichtag 31.12., regionale Tiefe: Gemeinden (Genesis table 11111-01-01-5).
# Encoding: Latin-1. Separator: semicolon.
# Format: date ; ags ; name ; area_qkm
# Header (metadata) and footer (notes) are dropped via the date-pattern filter.

raw = pd.read_csv(
    f"{DATA_DIR_RAW}/11111-01-01-5_area.csv",
    sep=";",
    encoding="latin-1",
    header=None,
    names=["date", "ags", "name", "area_qkm"],
    dtype=str,
    quoting=3,           # QUOTE_NONE: footer contains raw double-quotes
    on_bad_lines="skip",
)

# Keep only data rows: date matches DD.MM.YYYY, ags is purely numeric
raw = raw[raw["date"].str.match(r"^\d{2}\.\d{2}\.\d{4}$", na=False)].copy()
raw["ags"] = raw["ags"].str.strip()
raw = raw[raw["ags"].str.match(r"^\d+$", na=False)].copy()

# Extract year from date string (last 4 characters)
raw["year"] = raw["date"].str[-4:].astype(int)

# Parse area: German decimal comma, potential thousands-separator period
raw["area_qkm"] = (
    raw["area_qkm"]
    .str.strip()
    .str.replace(".", "", regex=False)   # remove thousands separator
    .str.replace(",", ".", regex=False)  # decimal comma → dot
)
raw["area_qkm"] = pd.to_numeric(raw["area_qkm"], errors="coerce")

# ── Filter to Gemeinde-level rows ────────────────────────────────────────────
# Three cases, detected by structural absence of children in the file itself:
#
# 1. 8-digit rows: regular Gemeinden in Landkreisen — include directly, EXCEPT
#    Berlin Bezirke (their 2-digit prefix has no 5-digit children, see below).
#
# 2. 5-digit rows with NO 8-digit children in the same year: kreisfreie Städte
#    (the city IS the only Gemeinde) — assign AGS8 = AGS5 + "000".
#    Bremen (04011, 04012) is handled here.
#
# 3. 2-digit rows with NO 5-digit children: Stadtstaaten where the file skips
#    straight from Bundesland to Bezirke (Hamburg) or Bezirke are 8-digit
#    (Berlin). Assign AGS8 = AGS2 + "000000" and exclude the spurious 8-digit
#    Bezirk rows (they have no matching AGS8 in INKAR).

two_digit  = raw[raw["ags"].str.len() == 2].copy()
five_digit = raw[raw["ags"].str.len() == 5].copy()
eight_digit = raw[raw["ags"].str.len() == 8].copy()

# 2-digit codes that have no 5-digit children = Stadtstaaten (Hamburg, Berlin)
five_prefixes_2 = set(five_digit["ags"].str[:2])
stadtstaaten = two_digit[~two_digit["ags"].isin(five_prefixes_2)].copy()
stadtstaaten["AGS8"] = stadtstaaten["ags"].str.zfill(2) + "000000"

# Exclude 8-digit rows whose 2-digit prefix is a Stadtstaat (Berlin Bezirke)
ss_prefixes = set(stadtstaaten["ags"])
gemeinden = eight_digit[~eight_digit["ags"].str[:2].isin(ss_prefixes)].copy()
gemeinden["AGS8"] = gemeinden["ags"]

# 5-digit rows with no 8-digit children = kreisfreie Städte
has_children_5 = set(zip(gemeinden["ags"].str[:5], gemeinden["year"]))
kreisfrei = five_digit[
    ~five_digit.apply(lambda r: (r["ags"], r["year"]) in has_children_5, axis=1)
].copy()
kreisfrei["AGS8"] = kreisfrei["ags"] + "000"

area = pd.concat(
    [gemeinden[["AGS8", "year", "area_qkm"]],
     kreisfrei[["AGS8", "year", "area_qkm"]],
     stadtstaaten[["AGS8", "year", "area_qkm"]]],
    ignore_index=True,
)
area = area.drop_duplicates(subset=["AGS8", "year"])

print(f"Stadtstaaten identified:       {stadtstaaten['ags'].nunique()} "
      f"({sorted(stadtstaaten['ags'].unique().tolist())})")
print(f"Kreisfreie Städte identified:  {kreisfrei['ags'].nunique()}")

# Zero area is a data artifact; treat as missing
area.loc[area["area_qkm"] == 0, "area_qkm"] = float("nan")

# ── Fill within AGS8 across years ────────────────────────────────────────────
# Area is time-invariant absent boundary reforms. Forward- then backward-fill
# within each AGS8 handles both mid-series and leading/trailing gaps.

area = area.sort_values(["AGS8", "year"]).reset_index(drop=True)
area["area_qkm"] = area.groupby("AGS8")["area_qkm"].transform(
    lambda s: s.ffill().bfill()
)

# ── Export ────────────────────────────────────────────────────────────────────

os.makedirs(f"{DATA_DIR_INTERMEDIATE}/area", exist_ok=True)
area.to_csv(f"{DATA_DIR_INTERMEDIATE}/area/area_ags8_panel.csv", index=False)

n_miss = area["area_qkm"].isna().sum()
print(f"Area panel shape:  {area.shape}")
print(f"AGS8 units:        {area['AGS8'].nunique()}")
print(f"Years:             {sorted(area['year'].unique().tolist())}")
print(f"Missing area:      {n_miss} rows ({100 * n_miss / len(area):.1f}%)")
print(
    area.groupby("year")["area_qkm"]
    .agg(n="count", missing=lambda s: s.isna().sum())
    .to_string()
)
