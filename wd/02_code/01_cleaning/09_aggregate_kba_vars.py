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

# ── Export ─────────────────────────────────────────────────────────────────────

output.to_csv(f"{INTERMEDIATE}/kba_ags8_panel.csv", index=False)

print(f"\nOutput shape: {output.shape}")
print(f"AGS8 units:   {output['AGS8'].nunique()}")
print(f"Years (B):    {sorted(wide_b['year'].unique().tolist())}")
print(f"Years (N):    {sorted(wide_n['year'].unique().tolist())}")
print(f"\nColumns:\n  " + "\n  ".join(output.columns.tolist()))
