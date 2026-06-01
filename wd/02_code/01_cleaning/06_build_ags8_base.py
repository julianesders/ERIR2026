import pandas as pd

from config import BASE_PATH
DATA_DIR_INTERMEDIATE = f"{BASE_PATH}/01_data/02_intermediate"

ags8 = pd.read_csv(
    f"{DATA_DIR_INTERMEDIATE}/anschriftenverzeichnis/anschriftenverzeichnis_ags8.csv",
    dtype=str,
    usecols=["AGS8", "AGS5"],
)
ags8["AGS8"] = ags8["AGS8"].str.zfill(8)
ags8["AGS5"] = ags8["AGS5"].str.zfill(5)

inkar = pd.read_csv(
    f"{DATA_DIR_INTERMEDIATE}/inkar/inkar_ags8_panel.csv",
    dtype={"AGS8": str, "AGS5": str},
)

inkar_units  = set(inkar["AGS8"].unique())
registry     = set(ags8["AGS8"].unique())
missing      = registry - inkar_units
extra        = inkar_units - registry
inkar_vars   = [c for c in inkar.columns if c not in ("AGS8", "AGS5", "year")]

print(f"Anschriftenverzeichnis units : {len(registry)}")
print(f"INKAR AGS8 units             : {len(inkar_units)}")
print(f"In registry, not in INKAR   : {len(missing)}")
if missing:
    print(f"  {sorted(missing)[:10]}{'...' if len(missing) > 10 else ''}")
print(f"In INKAR, not in registry   : {len(extra)}")
print(f"Years in INKAR panel         : {sorted(inkar['year'].unique().tolist())}")
print(f"INKAR variables ({len(inkar_vars)}): {inkar_vars}")
