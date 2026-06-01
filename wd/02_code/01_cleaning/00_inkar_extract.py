import pandas as pd

from config import BASE_PATH
INKAR_FILE  = f"{BASE_PATH}/01_data/05_inkar_raw/inkar_2025.csv"
OUT_FILE = f"{BASE_PATH}/01_data/02_intermediate/inkar/inkar_joint_panel.csv"

# ── Variable selection ─────────────────────────────────────────────────────────
# Use INKAR Kuerzel codes (short identifier, e.g. "xbev", "q_alo").
# AGS5_VARS → filtered from Raumbezug == "Kreise",   id truncated to 5 digits
# AGS8_VARS → filtered from Raumbezug == "Gemeinden", id kept as 8 digits

AGS5_VARS: list[str] = [
    # insert Kuerzel codes here
]

AGS8_VARS: list[str] = [
    "TN23-kataster_qkm", "xbev", "q_kaufkraft", "a_landwirtschaft", "q_gest_bev", "q_sach", "q_investZ", "q_pendlersaldo"
]

# ── Load ───────────────────────────────────────────────────────────────────────

ALL_VARS      = set(AGS5_VARS) | set(AGS8_VARS)
RAUMBEZUG_MAP = {"Kreise": AGS5_VARS, "Gemeinden": AGS8_VARS}

chunks = []
for chunk in pd.read_csv(
    INKAR_FILE,
    sep=";",
    dtype=str,
    usecols=["Kuerzel", "Raumbezug", "Kennziffer", "Zeitbezug", "Wert"],
    chunksize=2_000_000,
):
    mask = (
        chunk["Raumbezug"].isin(RAUMBEZUG_MAP) &
        chunk["Kuerzel"].isin(ALL_VARS)
    )
    if mask.any():
        chunks.append(chunk.loc[mask])

raw = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame(
    columns=["Kuerzel", "Raumbezug", "Kennziffer", "Zeitbezug", "Wert"]
)

# ── Parse ──────────────────────────────────────────────────────────────────────

raw["Wert"] = (
    raw["Wert"]
    .str.replace(".", "", regex=False)
    .str.replace(",", ".", regex=False)
)
raw["Wert"]      = pd.to_numeric(raw["Wert"],      errors="coerce")
raw["Zeitbezug"] = pd.to_numeric(raw["Zeitbezug"], errors="coerce")
raw              = raw[raw["Zeitbezug"].notna()]

# ── Build panels per level and pivot wide ─────────────────────────────────────

def _pivot(df: pd.DataFrame, id_len: int) -> pd.DataFrame:
    df = df.copy()
    df["id"] = df["Kennziffer"].str.strip().str.zfill(8).str[:id_len]
    return (
        df.pivot_table(index=["id", "Zeitbezug"], columns="Kuerzel",
                       values="Wert", aggfunc="first")
        .reset_index()
        .rename(columns={"Zeitbezug": "year"})
    )

parts = []

kreise_raw = raw[(raw["Raumbezug"] == "Kreise") & raw["Kuerzel"].isin(AGS5_VARS)]
if not kreise_raw.empty:
    p5 = _pivot(kreise_raw, id_len=5)
    p5.insert(1, "raumbezug", "Kreise")
    parts.append(p5)

gemeinden_raw = raw[(raw["Raumbezug"] == "Gemeinden") & raw["Kuerzel"].isin(AGS8_VARS)]
if not gemeinden_raw.empty:
    p8 = _pivot(gemeinden_raw, id_len=8)
    p8.insert(1, "raumbezug", "Gemeinden")
    parts.append(p8)

panel = (
    pd.concat(parts, ignore_index=True, sort=False)
    .sort_values(["raumbezug", "id", "year"])
    .reset_index(drop=True)
)
panel.columns.name = None

# ── Export ─────────────────────────────────────────────────────────────────────

panel.to_csv(OUT_FILE, index=False)

print(f"Panel shape:  {panel.shape}")
for rbz, grp in panel.groupby("raumbezug"):
    id_col = "id"
    print(f"  {rbz:12s}: {grp[id_col].nunique():5d} units × {grp['year'].nunique():2d} years"
          f"  cols: {[c for c in grp.columns if c not in ('id','raumbezug','year')]}")
