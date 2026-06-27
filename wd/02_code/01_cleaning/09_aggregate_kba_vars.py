import os
import pandas as pd
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parents[2]))
from config import BASE_PATH

INTERMEDIATE   = f"{BASE_PATH}/01_data/02_intermediate/kba"
CROSSWALK_DIR  = f"{BASE_PATH}/01_data/01_raw/crosswalks"

# Satzart codes — verify against the diagnostic output below if uncertain
SATZART_B = "B"   # Bestand
SATZART_N = "N"   # Neuzulassung (adjust if the actual code differs)

# ── Territorial reform crosswalk (KBA native codes → 2023 AGS8) ───────────────
# KBA delivers AGS8 codes as they existed in the year of the data.  These three
# functions parse the annual BBSR reform sheets, chain them into a forward map,
# and apply population-proportional weights before the variable-building step so
# all downstream aggregations and Abgang calculations operate on 2023 boundaries.

def _parse_cw_sheet(df):
    """
    Extract (ags8_from, ags8_to, bev_weight) from one annual reform sheet.

    Two column layouts exist:
      • Pre-2022 (10 cols): col 0 = from-code, col 3 = bev_prop, col 8 = to-code
      • 2022-2024 (10 cols, extra Regionalschlüssel column at pos 1):
                            col 0 = from-code, col 4 = bev_prop, col 7 = to-code
    Pre-2016 sheets also carry a sub-header row that must be stripped first.
    """
    if str(df.iloc[0, 0]).strip() == "Kennziffer":
        df = df.iloc[1:].reset_index(drop=True)

    col1 = str(df.columns[1]).lower()
    if "schlüssel" in col1 or "regional" in col1:
        col_from, col_bev, col_to = 0, 4, 7
    else:
        col_from, col_bev, col_to = 0, 3, 8

    ags_from = pd.to_numeric(df.iloc[:, col_from], errors="coerce")
    ags_to   = pd.to_numeric(df.iloc[:, col_to],   errors="coerce")
    bev_w    = pd.to_numeric(df.iloc[:, col_bev],  errors="coerce").fillna(1.0)

    mask = ags_from.notna() & ags_to.notna()
    return pd.DataFrame({
        "ags8_from":  ags_from[mask].astype(int).astype(str).str.zfill(8),
        "ags8_to":    ags_to[mask].astype(int).astype(str).str.zfill(8),
        "bev_weight": bev_w[mask].clip(0.0, 1.0).values,
    })


def _build_transitions(cw_dir, year_min, year_max):
    """
    Parse both crosswalk Excel files and return annual transition DataFrames
    for year_min ≤ year_from ≤ year_max.
    Sheet naming: "2013-2014" → year_from=2013; "2020" → year_from=2020.
    """
    files = [
        os.path.join(cw_dir, "ref-gemeinden-2010-2020.xlsx"),
        os.path.join(cw_dir, "ref-gemeinden-2020-2024.xlsx"),
    ]
    transitions = {}
    for fpath in files:
        wb = pd.read_excel(fpath, sheet_name=None)
        for sheet, df in wb.items():
            year_from = int(sheet.split("-")[0]) if "-" in sheet else int(sheet)
            if year_from < year_min or year_from > year_max:
                continue
            transitions[year_from] = _parse_cw_sheet(df)
    return transitions


def _build_forward_map(transitions, year_min, year_max):
    """
    Chain annual transitions into a forward map:
        {year: {ags8_native: {ags8_2023: weight}}}
    Weights are composed multiplicatively across splits; mergers keep weight=1.
    Codes absent from the crosswalk for a given year are handled by an identity
    fallback in apply_kba_crosswalk (left-join NaN → filled with source code).
    """
    # Convert DataFrames to nested dicts for fast lookup
    trans = {}
    for yr, df in transitions.items():
        d = {}
        for _, row in df.iterrows():
            af, at, w = row["ags8_from"], row["ags8_to"], row["bev_weight"]
            if af not in d:
                d[af] = {}
            d[af][at] = d[af].get(at, 0.0) + w
        trans[yr] = d

    fwd = {}
    # Seed from the last transition (year_max → year_max+1 = 2023)
    fwd[year_max] = {af: dict(succ) for af, succ in trans.get(year_max, {}).items()}

    for yr in range(year_max - 1, year_min - 1, -1):
        fwd[yr] = {}
        for af, succ_mid in trans.get(yr, {}).items():
            final = {}
            for ags_mid, w1 in succ_mid.items():
                if ags_mid in fwd.get(yr + 1, {}):
                    for ags_final, w2 in fwd[yr + 1][ags_mid].items():
                        final[ags_final] = final.get(ags_final, 0.0) + w1 * w2
                else:
                    final[ags_mid] = final.get(ags_mid, 0.0) + w1
            fwd[yr][af] = final

    return fwd


def apply_kba_crosswalk(panel, cw_dir):
    """
    Remap native-year AGS8 codes in kba_panel to 2023 boundaries.
    - Mergers:  sum counts of all old units that converge on the same 2023 code.
    - Splits:   multiply anzahl by bevölkerungsprop. weight before summing.
    - Identity: codes not in the crosswalk are passed through unchanged.
    Returns a panel with the same columns as input but on 2023-boundary AGS8.
    """
    year_min = int(panel["year"].min())
    year_max = int(panel["year"].max()) - 1   # last required transition: max-1 → max

    transitions = _build_transitions(cw_dir, year_min, year_max)
    fwd = _build_forward_map(transitions, year_min, year_max)

    # Flatten forward map to a DataFrame for a vectorised merge
    xw_rows = [
        {"year": yr, "AGS8": ags_from, "AGS8_2023": ags_to, "weight": w}
        for yr, ags_map in fwd.items()
        for ags_from, targets in ags_map.items()
        for ags_to, w in targets.items()
    ]
    xw_df = (
        pd.DataFrame(xw_rows)
        if xw_rows
        else pd.DataFrame(columns=["year", "AGS8", "AGS8_2023", "weight"])
    )

    merged = panel.merge(xw_df, on=["year", "AGS8"], how="left")
    merged["AGS8_2023"] = merged["AGS8_2023"].fillna(merged["AGS8"])   # identity fallback
    merged["weight"]    = merged["weight"].fillna(1.0)
    merged["anzahl"]    = (merged["anzahl"] * merged["weight"])

    grp_cols = ["AGS8_2023", "year", "satzart", "energiequelle", "zulassungsart"]
    result = (
        merged
        .groupby(grp_cols, as_index=False)["anzahl"]
        .sum()
        .rename(columns={"AGS8_2023": "AGS8"})
    )
    result["anzahl"] = result["anzahl"].round().astype("int32")

    # Diagnostics
    n_remapped = (merged["AGS8"] != merged["AGS8_2023"]).sum()
    n_splits   = int((merged["weight"] < 1.0).sum())
    print(f"  AGS8 units before crosswalk: {panel['AGS8'].nunique()}")
    print(f"  AGS8 units after  crosswalk: {result['AGS8'].nunique()}")
    print(f"  Rows remapped (code changed): {n_remapped}")
    print(f"  Rows from splits (weight<1):  {n_splits}")

    return result


# ── Load ───────────────────────────────────────────────────────────────────────

panel = pd.read_csv(
    f"{INTERMEDIATE}/kba_panel.csv",
    dtype={"AGS8": str, "satzart": str, "energiequelle": str, "zulassungsart": str},
)

print(f"Loaded kba_panel: {panel.shape}")
print(f"Satzarten found:  {sorted(panel['satzart'].unique().tolist())}")
print(f"Energiequellen:   {sorted(panel['energiequelle'].unique().tolist())}")
print(f"Zulassungsarten:  {sorted(panel['zulassungsart'].unique().tolist())}")

# ── Apply crosswalk ────────────────────────────────────────────────────────────

print("\nApplying territorial reform crosswalk...")
panel = apply_kba_crosswalk(panel, CROSSWALK_DIR)

# ── Grouping definitions ───────────────────────────────────────────────────────

FUEL_GROUPS = {
    "benzin_diesel": ["01", "02"],
    "elektro":       ["04"],
    "hybrid":        ["05", "06"],
    "gas_other":     ["03", "07"],
}

OWNERSHIP = {
    "overall":   None,   # no zulassungsart filter
    "private":   "01",
    "corporate": "02",
}

IDX = ["AGS8", "year"]

# ── Variable builder ───────────────────────────────────────────────────────────

def build_vars(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    spine = df[IDX].drop_duplicates().reset_index(drop=True)

    def _agg(fuel_codes=None, own_code=None) -> pd.DataFrame:
        sub = df
        if fuel_codes is not None:
            sub = sub[sub["energiequelle"].isin(fuel_codes)]
        if own_code is not None:
            sub = sub[sub["zulassungsart"] == own_code]
        return sub.groupby(IDX, as_index=False)["anzahl"].sum()

    specs: dict[str, dict] = {}
    for own_name, own_code in OWNERSHIP.items():
        specs[f"{prefix}_total_{own_name}"] = dict(fuel_codes=None, own_code=own_code)
    for fuel_name, fuel_codes in FUEL_GROUPS.items():
        for own_name, own_code in OWNERSHIP.items():
            specs[f"{prefix}_{fuel_name}_{own_name}"] = dict(fuel_codes=fuel_codes, own_code=own_code)

    wide = spine.copy()
    for col_name, kwargs in specs.items():
        agg = _agg(**kwargs).rename(columns={"anzahl": col_name})
        wide = wide.merge(agg, on=IDX, how="left")
        wide[col_name] = wide[col_name].fillna(0).astype("int32")

    for own_name in OWNERSHIP:
        wide[f"{prefix}_ev_share_{own_name}"] = (
            wide[f"{prefix}_elektro_{own_name}"] / wide[f"{prefix}_total_{own_name}"]
        )

    return wide

# ── Build per satzart and join ─────────────────────────────────────────────────

wide_b = build_vars(panel[panel["satzart"] == SATZART_B], prefix="B")
wide_n = build_vars(panel[panel["satzart"] == SATZART_N], prefix="N")

output = wide_b.merge(wide_n, on=IDX, how="outer")
output = output.sort_values(IDX).reset_index(drop=True)

# ── Abgänge: A(t) = B(t-1) + N(t) – B(t) ────────────────────────────────────
# Computable for years where all three exist: 2015–2022
# (B available 2014–2023, N available 2013–2022 → overlap with B(t-1) starts 2015)

b_count_cols = [c for c in wide_b.columns if c not in IDX and "ev_share" not in c]
dim_names    = [c[2:] for c in b_count_cols]  # strip "B_" → e.g. "total_overall"

b_lagged = wide_b[IDX + b_count_cols].copy()
b_lagged["year"] = b_lagged["year"] + 1  # B(t-1) appears at year t
b_lagged = b_lagged.rename(columns={c: f"Blag_{c[2:]}" for c in b_count_cols})

abgang = wide_n.merge(wide_b[IDX + b_count_cols], on=IDX, how="inner")
abgang = abgang.merge(b_lagged, on=IDX, how="inner")

for dim in dim_names:
    abgang[f"A_{dim}"] = (
        abgang[f"Blag_{dim}"] + abgang[f"N_{dim}"] - abgang[f"B_{dim}"]
    ).astype("int32")

for own_name in OWNERSHIP:
    abgang[f"A_ev_share_{own_name}"] = (
        abgang[f"A_elektro_{own_name}"] / abgang[f"A_total_{own_name}"]
    )

a_cols = [c for c in abgang.columns if c.startswith("A_")]
output = output.merge(abgang[IDX + a_cols], on=IDX, how="left")
output = output.sort_values(IDX).reset_index(drop=True)

# ── Export ─────────────────────────────────────────────────────────────────────

output.to_csv(f"{INTERMEDIATE}/kba_ags8_panel.csv", index=False)

print(f"\nOutput shape: {output.shape}")
print(f"AGS8 units:   {output['AGS8'].nunique()}")
print(f"Years (B):    {sorted(wide_b['year'].unique().tolist())}")
print(f"Years (N):    {sorted(wide_n['year'].unique().tolist())}")
print(f"Years (A):    {sorted(abgang['year'].unique().tolist())}")
print(f"\nAbgänge sanity (A_total_overall, negative rows): "
      f"{(abgang['A_total_overall'] < 0).sum()} / {len(abgang)}")
print(f"\nColumns:\n  " + "\n  ".join(output.columns.tolist()))
