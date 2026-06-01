import pandas as pd

BASE_PATH    = "/Users/julian/Documents/ERIR2026/wd"
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

# ── Backward-fill area variable ────────────────────────────────────────────────
# Bodenfläche (TN23-kataster_qkm) is only available from 2016 in INKAR.
# Land area is constant absent boundary reforms, so the earliest observed value
# is propagated back to fill pre-2016 years.
BFILL_VARS = ["TN23-kataster_qkm"]
for v in BFILL_VARS:
    if v in panel.columns:
        panel[v] = panel.groupby("AGS8")[v].transform("bfill")

# ── Export ─────────────────────────────────────────────────────────────────────

panel.to_csv(f"{INTERMEDIATE}/inkar/inkar_ags8_panel.csv", index=False)

print(f"Panel shape:   {panel.shape}")
print(f"AGS8 units:    {panel['AGS8'].nunique()}")
print(f"AGS5 units:    {panel['AGS5'].nunique()}")
print(f"Years:         {sorted(panel['year'].unique().tolist())}")
print(f"Kreise vars:   {kreise_vars}")
print(f"Gemeinden vars:{gemeinden_vars}")
