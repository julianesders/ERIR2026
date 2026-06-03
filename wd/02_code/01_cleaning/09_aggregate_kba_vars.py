import pandas as pd
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parents[2]))
from config import BASE_PATH

INTERMEDIATE = f"{BASE_PATH}/01_data/02_intermediate/kba"

# Satzart codes — verify against the diagnostic output below if uncertain
SATZART_B = "B"   # Bestand
SATZART_N = "N"   # Neuzulassung (adjust if the actual code differs)

# ── Load ───────────────────────────────────────────────────────────────────────

panel = pd.read_csv(
    f"{INTERMEDIATE}/kba_panel.csv",
    dtype={"AGS8": str, "satzart": str, "energiequelle": str, "zulassungsart": str},
)

print(f"Loaded kba_panel: {panel.shape}")
print(f"Satzarten found:  {sorted(panel['satzart'].unique().tolist())}")
print(f"Energiequellen:   {sorted(panel['energiequelle'].unique().tolist())}")
print(f"Zulassungsarten:  {sorted(panel['zulassungsart'].unique().tolist())}")

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
