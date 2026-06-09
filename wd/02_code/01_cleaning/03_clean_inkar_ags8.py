import pandas as pd

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parents[2]))
from config import BASE_PATH
INTERMEDIATE = f"{BASE_PATH}/01_data/02_intermediate"

# ── Load ───────────────────────────────────────────────────────────────────────

joint = pd.read_csv(
    f"{INTERMEDIATE}/inkar/inkar_joint_panel.csv",
    dtype={"id": str},
)
joint["id"] = joint["id"].str.strip()

ags_xw = pd.read_csv(
    f"{INTERMEDIATE}/anschriftenverzeichnis/anschriftenverzeichnis_ags8.csv",
    dtype=str,
    usecols=["AGS8", "AGS5"],
)
ags_xw["AGS8"] = ags_xw["AGS8"].str.zfill(8)
ags_xw["AGS5"] = ags_xw["AGS5"].str.zfill(5)

# ── Split by level ─────────────────────────────────────────────────────────────

kreise    = joint[joint["raumbezug"] == "Kreise"].drop(columns="raumbezug").copy()
gemeinden = joint[joint["raumbezug"] == "Gemeinden"].drop(columns="raumbezug").copy()

kreise["id"]    = kreise["id"].str.zfill(5)
gemeinden["id"] = gemeinden["id"].str.zfill(8)

kreise    = kreise.rename(columns={"id": "AGS5"})
gemeinden = gemeinden.rename(columns={"id": "AGS8"})

# Variable columns per level — ignore columns that are all-NaN for that level
# (they belong to the other level in the joint panel)
kreise_vars    = [c for c in kreise.columns    if c not in ("AGS5", "year") and kreise[c].notna().any()]
gemeinden_vars = [c for c in gemeinden.columns if c not in ("AGS8", "year") and gemeinden[c].notna().any()]

# ── Spine: all AGS8 × union of all years present in either level ───────────────

all_years = sorted(
    set(kreise["year"].dropna().astype(int)) |
    set(gemeinden["year"].dropna().astype(int))
)

spine = (
    ags_xw[["AGS8", "AGS5"]]
    .assign(_k=1)
    .merge(pd.DataFrame({"year": all_years, "_k": 1}), on="_k")
    .drop(columns="_k")
)

# ── Broadcast Kreise → every AGS8 in that Kreis ───────────────────────────────

panel = spine.merge(
    kreise[["AGS5", "year"] + kreise_vars],
    on=["AGS5", "year"],
    how="left",
)

# ── Merge Gemeinden variables ──────────────────────────────────────────────────

panel = panel.merge(
    gemeinden[["AGS8", "year"] + gemeinden_vars],
    on=["AGS8", "year"],
    how="left",
)

panel = panel.sort_values(["AGS8", "year"]).reset_index(drop=True)

# Area (TN23-kataster_qkm) is superseded by the dedicated area panel built in
# 11_clean_area_ags8.py and merged in 02_merge_emk_panel_ags8.py.
panel = panel.drop(columns=["TN23-kataster_qkm"], errors="ignore")

# ── Export ─────────────────────────────────────────────────────────────────────

panel.to_csv(f"{INTERMEDIATE}/inkar/inkar_ags8_panel.csv", index=False)

print(f"Panel shape:   {panel.shape}")
print(f"AGS8 units:    {panel['AGS8'].nunique()}")
print(f"AGS5 units:    {panel['AGS5'].nunique()}")
print(f"Years:         {sorted(panel['year'].unique().tolist())}")
print(f"Kreise vars:   {kreise_vars}")
print(f"Gemeinden vars:{gemeinden_vars}")
