import pandas as pd

BASE_PATH = "/Users/julian/Documents/ERIR2026/wd"

DATA_DIR_RAW          = f"{BASE_PATH}/01_data/01_raw"
DATA_DIR_INTERMEDIATE = f"{BASE_PATH}/01_data/02_intermediate"

# Load only relevant columns + rows from excel
anschriften = pd.read_excel(
    f"{DATA_DIR_RAW}/anschriftenverzeichnis-5119101237005.xlsx",
    sheet_name="Anschriften_31_01_2023",
    skiprows=8,
    skipfooter=8,
    dtype=str,
    usecols="A:H",
    names=[
        "AGS2", "Bundesland", "Satzart",
        "TKZ", "TKZ_name", "ARS",
        "AGS8", "Gemeinde_Stadt"
    ]
)

# Keep only the first instance of each AGS number (all duplicates only differ for the non-relevant columns)
anschriften = anschriften[anschriften["AGS8"].isna() | ~anschriften["AGS8"].duplicated()]
anschriften["AGS2"] = anschriften["AGS2"].str.zfill(2)


########################################
########        GEMEINDEN       ########
########################################


# Filter by Gemeinden
gemeinden = anschriften.dropna(axis=0, subset="AGS8")

# Extract AGS3 and AGS5
gemeinden["AGS3"] = gemeinden["AGS8"].str[:3].str.zfill(3)
gemeinden["AGS5"] = gemeinden["AGS8"].str[:5].str.zfill(5)


########################################
########        Kreise          ########
########################################


# Filter by Kreise
kreise = anschriften[anschriften["Satzart"] == "40"]      # Satzart 40 is equal to Kreise

# Give all Duplicate entries (where both Kreis and kreisfreie Stadt exists) a suffix on the Stadt entry
dupes = kreise["Gemeinde_Stadt"].duplicated(keep=False)
mask = dupes & (kreise["TKZ_name"] == "Kreisfreie Stadt")
kreise.loc[mask, "Gemeinde_Stadt"] = kreise.loc[mask, "Gemeinde_Stadt"] + ", Stadt"

# Create AGS8: for kreisfreie Städte, AGS8 = ARS + "000"
mask_kfs = kreise["TKZ_name"] == "Kreisfreie Stadt"
kreise.loc[mask_kfs, "AGS8"] = kreise.loc[mask_kfs, "ARS"] + "000"

# Extract AGS3 and AGS5
kreise["AGS3"] = kreise["ARS"].str[:3].str.zfill(3)
kreise["AGS5"] = kreise["ARS"].str[:5].str.zfill(5)


########################################
########    Append and Export   ########
########################################


# Make it prettier
cols = ["Bundesland", "Gemeinde_Stadt", "TKZ_name", "AGS8", "AGS5", "AGS3", "AGS2"]
gemeinden = gemeinden.drop(columns=["Satzart", "TKZ", "ARS"])[cols]
kreise = kreise[cols]

# Append kreise and gemeinden; kreise takes priority for duplicates on AGS8
combined = (
    pd.concat([kreise, gemeinden], ignore_index=True)
    .assign(AGS8=lambda df: df["AGS8"].astype(str).replace("nan", pd.NA))
    .sort_values("AGS8", na_position="last")
    .pipe(lambda df: df[df["AGS8"].isna() | ~df["AGS8"].duplicated(keep="first")])
    .reset_index(drop=True)
)

OUT = f"{DATA_DIR_INTERMEDIATE}/anschriftenverzeichnis"
combined.to_csv(f"{OUT}/anschriftenverzeichnis_combined.csv", sep=",", index=False, encoding="utf-8")
gemeinden.to_csv(f"{OUT}/anschriftenverzeichnis_ags8.csv",   sep=",", index=False, encoding="utf-8")
