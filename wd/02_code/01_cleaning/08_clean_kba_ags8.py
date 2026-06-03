import os
import zipfile
import pandas as pd

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parents[2]))
from config import BASE_PATH
DELIVERY     = f"{BASE_PATH}/01_data/00_delivery/kba"
RAW          = f"{BASE_PATH}/01_data/01_raw/kba"
INTERMEDIATE = f"{BASE_PATH}/01_data/02_intermediate/kba"

os.makedirs(INTERMEDIATE, exist_ok=True)

# KBA Sonstige codes: fab_code="000" = unknown brand, mod_code="0000" = unknown model
SONSTIGE_FAB_CODE = "000"
SONSTIGE_MOD_CODE = "0000"

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
MODEL_COLS = ["jahr", "AGS8", "fab_code", "fab_text", "mod_code", "mod_text",
              "energiequelle", "zulassungsart"]


def _read_fwf(filepath: str) -> pd.DataFrame:
    return pd.read_fwf(
        filepath,
        colspecs=COLSPECS,
        names=COLNAMES,
        dtype=DTYPES,
        encoding="latin-1",
        header=None,
    )


def _clean_ags8(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["AGS8"] = df["AGS8"].str.strip().str.zfill(8)
    return df[
        ~df["AGS8"].str.endswith("999") &
        (df["AGS8"] != "99999999")
    ].copy()


def _process_folder(folder: str, keep_model: bool = False) -> pd.DataFrame:
    dfs = []
    for fname in sorted(os.listdir(folder)):
        if not fname.endswith(".txt"):
            continue
        print(f"  reading {fname}...")
        df = _read_fwf(os.path.join(folder, fname))
        if keep_model:
            df = df.groupby(MODEL_COLS, observed=True)["anzahl"].sum().reset_index()
        else:
            df = df.drop(columns=["fab_code", "fab_text", "mod_code", "mod_text"])
            df = df.groupby(GROUP_COLS, observed=True)["anzahl"].sum().reset_index()
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)


# ── Read and aggregate ─────────────────────────────────────────────────────────

print("Processing Bestandsdaten...")
bestaende = _process_folder(f"{RAW}/bestaende")

print("\nProcessing Neuzulassungen (aggregated)...")
neuzulassungen = _process_folder(f"{RAW}/neuzulassungen")

print("\nProcessing Neuzulassungen (model level)...")
neuz_model = _process_folder(f"{RAW}/neuzulassungen", keep_model=True)


# ── Main panel ─────────────────────────────────────────────────────────────────

panel = pd.concat([bestaende, neuzulassungen], ignore_index=True)
panel = panel.groupby(GROUP_COLS, observed=True)["anzahl"].sum().reset_index()
panel = _clean_ags8(panel)
panel = (
    panel
    .rename(columns={"jahr": "year"})
    .astype({"anzahl": "int32"})
    .sort_values(["AGS8", "year", "satzart", "energiequelle", "zulassungsart"])
    .reset_index(drop=True)
)


# ── Model-level Neuzulassungen ─────────────────────────────────────────────────

neuz_model = _clean_ags8(neuz_model)
neuz_model["fab_code"] = neuz_model["fab_code"].str.strip().str.zfill(3)
neuz_model["mod_code"] = neuz_model["mod_code"].str.strip().str.zfill(4)
neuz_model["fab_text"] = neuz_model["fab_text"].str.strip()
neuz_model["mod_text"] = neuz_model["mod_text"].str.strip()
neuz_model = (
    neuz_model
    .rename(columns={"jahr": "year"})
    .astype({"anzahl": "int32"})
    .sort_values(["AGS8", "year", "fab_code", "mod_code", "energiequelle", "zulassungsart"])
    .reset_index(drop=True)
)


# ── Sonstige diagnostic ────────────────────────────────────────────────────────

total = neuz_model["anzahl"].sum()
n_sonstige_fab = neuz_model[
    neuz_model["fab_code"] == SONSTIGE_FAB_CODE
]["anzahl"].sum()
n_sonstige_mod = neuz_model[
    (neuz_model["fab_code"] != SONSTIGE_FAB_CODE) &
    (neuz_model["mod_code"] == SONSTIGE_MOD_CODE)
]["anzahl"].sum()
n_matchable = total - n_sonstige_fab - n_sonstige_mod

print(f"\nSonstige breakdown (Neuzulassungen, post AGS8 cleaning):")
print(f"  Sonstige brand  (fab_code={SONSTIGE_FAB_CODE}):         "
      f"{n_sonstige_fab:>10,}  ({n_sonstige_fab / total * 100:.2f}%)")
print(f"  Sonstige model in known brand (mod_code={SONSTIGE_MOD_CODE}):  "
      f"{n_sonstige_mod:>10,}  ({n_sonstige_mod / total * 100:.2f}%)")
print(f"  Directly matchable (neither Sonstige):        "
      f"{n_matchable:>10,}  ({n_matchable / total * 100:.2f}%)")


# ── Export ─────────────────────────────────────────────────────────────────────

panel.to_csv(f"{INTERMEDIATE}/kba_panel.csv", index=False, encoding="utf-8-sig")
neuz_model.to_csv(f"{INTERMEDIATE}/kba_neuz_model_panel.csv", index=False, encoding="utf-8-sig")

print(f"\nkba_panel:                  {panel.shape}")
print(f"  AGS8 units:               {panel['AGS8'].nunique()}")
print(f"  Years:                    {sorted(panel['year'].unique().tolist())}")
print(f"  Satzarten:                {sorted(panel['satzart'].unique().tolist())}")
print(f"  Energiequellen:           {sorted(panel['energiequelle'].unique().tolist())}")
print(f"  Zulassungsarten:          {sorted(panel['zulassungsart'].unique().tolist())}")
print(f"\nkba_neuz_model_panel:       {neuz_model.shape}")
print(f"  AGS8 units:               {neuz_model['AGS8'].nunique()}")
print(f"  Years:                    {sorted(neuz_model['year'].unique().tolist())}")
print(f"  Unique brand codes:       {neuz_model['fab_code'].nunique()}")
print(f"  Unique model codes:       {neuz_model['mod_code'].nunique()}")
print(f"  Unique brand × model:     {neuz_model[['fab_code', 'mod_code']].drop_duplicates().shape[0]}")
