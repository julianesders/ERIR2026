import os
import zipfile
import pandas as pd

from config import BASE_PATH
DELIVERY     = f"{BASE_PATH}/01_data/00_delivery/kba"
RAW          = f"{BASE_PATH}/01_data/01_raw/kba"
INTERMEDIATE = f"{BASE_PATH}/01_data/02_intermediate/kba"

os.makedirs(INTERMEDIATE, exist_ok=True)

# ── Unpack delivery ────────────────────────────────────────────────────────────

def _unzip_delivery(folder: str, destination: str) -> None:
    for fname in sorted(os.listdir(folder)):
        if not fname.endswith(".zip"):
            continue
        if "Bestand" in fname:
            subfolder = "bestaende"
        elif "Neuzulassung" in fname:
            subfolder = "neuzulassungen"
        else:
            print(f"  WARNING: unrecognised delivery file: {fname}")
            continue
        out_dir = os.path.join(destination, subfolder)
        os.makedirs(out_dir, exist_ok=True)
        with zipfile.ZipFile(os.path.join(folder, fname), "r") as z:
            z.extractall(out_dir)
        print(f"  unpacked: {fname} → {subfolder}/")

_unzip_delivery(DELIVERY, RAW)

# ── FWF spec ───────────────────────────────────────────────────────────────────
# Source: KBA data delivery (Bestandsdaten / Neuzulassungen)
# Fixed-width layout per column:
#   satzart (1), jahr (4), AGS8 (8), fab_code (3), fab_text (15),
#   mod_code (4), mod_text (34), energiequelle (2), zulassungsart (2), anzahl (7)

COLSPECS = [
    (0, 1), (1, 5), (5, 13), (13, 16), (16, 31),
    (31, 35), (35, 69), (69, 71), (71, 73), (73, 80),
]
COLNAMES = [
    "satzart", "jahr", "AGS8", "fab_code", "fab_text",
    "mod_code", "mod_text", "energiequelle", "zulassungsart", "anzahl",
]
DTYPES = {
    "satzart":       "category",
    "jahr":          "int16",
    "AGS8":          "string",
    "fab_code":      "string",
    "fab_text":      "string",
    "mod_code":      "string",
    "mod_text":      "string",
    "energiequelle": "category",
    "zulassungsart": "category",
    "anzahl":        "int16",
}

GROUP_COLS = ["satzart", "jahr", "AGS8", "energiequelle", "zulassungsart"]

def _read_fwf(filepath: str) -> pd.DataFrame:
    return pd.read_fwf(
        filepath,
        colspecs=COLSPECS,
        names=COLNAMES,
        dtype=DTYPES,
        encoding="latin-1",
        header=None,
    )

def _process_folder(folder: str) -> pd.DataFrame:
    dfs = []
    for fname in sorted(os.listdir(folder)):
        if not fname.endswith(".txt"):
            continue
        print(f"  reading {fname}...")
        df = _read_fwf(os.path.join(folder, fname))
        df = df.drop(columns=["fab_code", "fab_text", "mod_code", "mod_text"])
        df = df.groupby(GROUP_COLS, observed=True)["anzahl"].sum().reset_index()
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)

# ── Read and aggregate ─────────────────────────────────────────────────────────

print("Processing Bestandsdaten...")
bestaende = _process_folder(f"{RAW}/bestaende")

print("\nProcessing Neuzulassungen...")
neuzulassungen = _process_folder(f"{RAW}/neuzulassungen")

panel = pd.concat([bestaende, neuzulassungen], ignore_index=True)
panel = panel.groupby(GROUP_COLS, observed=True)["anzahl"].sum().reset_index()

# ── Clean AGS8 ─────────────────────────────────────────────────────────────────
# Drop "Sonstige" entries: AGS8 ending in 999 (unassigned within Kreis)
# and 99999999 (nationally unassigned).

panel["AGS8"] = panel["AGS8"].str.strip().str.zfill(8)
panel = panel[
    ~panel["AGS8"].str.endswith("999") &
    (panel["AGS8"] != "99999999")
].copy()

panel = (
    panel
    .rename(columns={"jahr": "year"})
    .astype({"anzahl": "int32"})
    .sort_values(["AGS8", "year", "satzart", "energiequelle", "zulassungsart"])
    .reset_index(drop=True)
)

# ── Export ─────────────────────────────────────────────────────────────────────

panel.to_csv(f"{INTERMEDIATE}/kba_panel.csv", index=False, encoding="utf-8-sig")

print(f"\nPanel shape:    {panel.shape}")
print(f"AGS8 units:     {panel['AGS8'].nunique()}")
print(f"Years:          {sorted(panel['year'].unique().tolist())}")
print(f"Satzarten:      {sorted(panel['satzart'].unique().tolist())}")
print(f"Energiequellen: {sorted(panel['energiequelle'].unique().tolist())}")
print(f"Zulassungsarten:{sorted(panel['zulassungsart'].unique().tolist())}")
