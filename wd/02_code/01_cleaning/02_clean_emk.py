import pandas as pd

from config import BASE_PATH
DATA_DIR_RAW          = f"{BASE_PATH}/01_data/01_raw"
DATA_DIR_INTERMEDIATE = f"{BASE_PATH}/01_data/02_intermediate"

# Load the scraped data
emk_raw = pd.read_csv(f"{DATA_DIR_RAW}/scraping_raw_60/emkonzepte_database.csv")

# Drop unnecessary columns and convert all columns to string
emk_raw = emk_raw.drop(columns=["title", "url", "pdf_url", "source_page", "raw_text"]).astype("string")

# Extract funding duration columns (and convert to datetime)
emk_raw[["laufzeit_start", "laufzeit_end"]] = emk_raw["laufzeit"].str.split(" - ", expand=True)
emk_raw["laufzeit_start"] = pd.to_datetime(emk_raw["laufzeit_start"], dayfirst=True)
emk_raw["laufzeit_end"] = pd.to_datetime(emk_raw["laufzeit_end"], dayfirst=True)

# Extract funding duration (in days)
emk_raw["laufzeit_days"] = (emk_raw["laufzeit_end"] - emk_raw["laufzeit_start"]).dt.days

emk_raw = emk_raw.drop(columns=["laufzeit"])

# Convert funding amount columns into floats
emk_raw[["gesamtmittel", "bundesmittel"]] = emk_raw[["gesamtmittel", "bundesmittel"]].apply(
    lambda col: col
    .str.replace("€", "", regex=False)
    .str.replace(".", "", regex=False)
    .str.replace(",", ".", regex=False)
    .str.strip()
    .astype(float)
)

# Extract dummies from schlagwoerter variable
tags = emk_raw["schlagwoerter"].str.get_dummies(sep="; ")

tags.columns = (tags.columns
    .str.strip()
    .str.lower()
    .str.replace("ä", "ae", regex=False)
    .str.replace("ö", "oe", regex=False)
    .str.replace("ü", "ue", regex=False)
    .str.replace("ß", "ss", regex=False)
    .str.replace("(", "", regex=False)
    .str.replace(")", "", regex=False)
    .str.replace("/", "_", regex=False)
    .str.replace("-", "_", regex=False)
    .str.replace(" ", "_", regex=False)
)

tags = tags.add_prefix("tag_")
emk_raw = pd.concat([emk_raw, tags], axis=1).drop(columns=["schlagwoerter"])

# Extract raumtyp dummies
raumtyp_dummies = pd.get_dummies(emk_raw["raumtyp"], prefix="space").rename(columns=lambda c: c
    .replace("ä", "ae")
    .replace("ö", "oe")
    .replace("ü", "ue")
    .replace("ß", "ss")
    .replace("(", "")
    .replace(")", "")
    .replace("/", "_")
    .replace("-", "_")
    .replace(" ", "_")
    .replace(",", "")
    .lower()
).astype(int)

emk_raw = pd.concat([emk_raw, raumtyp_dummies], axis=1)

# Manual empfaenger replacement for Nordseeheilbad Borkum GmbH, as it is a 100% subsidiary of the city of Borkum
emk_raw["empfaenger"] = emk_raw["empfaenger"].replace({"GebietskörperschaftUnternehmen": "Gebietskörperschaft"})

# Rename antragsteller of format "Land/Stadt [NAME], vertreten durch __"
emk_raw["antragsteller"] = emk_raw["antragsteller"].where(
    emk_raw["empfaenger"] != "Gebietskörperschaft",
    emk_raw["antragsteller"].str.split(",").str[0]
)

# Create antragsteller_id variable
emk_raw["antragsteller_id"] = emk_raw["antragsteller"].astype("category").cat.codes

# Rearrange order of columns for a clean export
first_cols = ["foerderkennzeichen", "antragsteller_id", "antragsteller", "empfaenger", "bundesland",
              "raumtyp", "laufzeit_start", "laufzeit_end", "laufzeit_days",
              "gesamtmittel", "bundesmittel"]

remaining_cols = [c for c in emk_raw.columns if c not in first_cols]

emk_raw = emk_raw[first_cols + remaining_cols]

# Export to csv
emk_raw.to_csv(f"{DATA_DIR_INTERMEDIATE}/emk/emkonzepte_cleaned.csv", sep=",", index=False)
