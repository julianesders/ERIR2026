import pandas as pd
import geopandas as gpd

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parents[2]))
from config import BASE_PATH
DATA_DIR_RAW          = f"{BASE_PATH}/01_data/01_raw"
DATA_DIR_INTERMEDIATE = f"{BASE_PATH}/01_data/02_intermediate"
GPKG_PATH = (
    f"{BASE_PATH}/01_data/04_shapefiles/"
    "vg250_01-01.utm32s.gpkg.ebenen/vg250_ebenen_0101/DE_VG250.gpkg"
)

# ── Load BNetzA Ladesäulenregister ────────────────────────────────────────────

raw = pd.read_excel(
    f"{DATA_DIR_RAW}/Ladesaeulenregister_BNetzA_2026-03-25.xlsx",
    header=10,
    usecols=[
        "Ladeeinrichtungs-ID",
        "Inbetriebnahmedatum",
        "Anzahl Ladepunkte",
        "Breitengrad",
        "Längengrad",
    ],
)

raw = raw.rename(columns={
    "Ladeeinrichtungs-ID": "station_id",
    "Anzahl Ladepunkte":   "n_ladepunkte_station",
})

raw["Inbetriebnahmedatum"] = pd.to_datetime(raw["Inbetriebnahmedatum"], errors="coerce")
raw["year_commissioned"]   = raw["Inbetriebnahmedatum"].dt.year.astype("Int64")

# German decimal format uses comma; convert to float
raw["lat"] = raw["Breitengrad"].astype(str).str.replace(",", ".", regex=False)
raw["lon"] = raw["Längengrad"].astype(str).str.replace(",", ".", regex=False)
raw["lat"] = pd.to_numeric(raw["lat"], errors="coerce")
raw["lon"] = pd.to_numeric(raw["lon"], errors="coerce")

raw = raw.dropna(subset=["lat", "lon", "year_commissioned"]).copy()
raw["year_commissioned"] = raw["year_commissioned"].astype(int)

# ── Load VG250 Gemeinde boundary layer ────────────────────────────────────────

gem = gpd.read_file(GPKG_PATH, layer="vg250_gem")
gem = gem[gem["GF"] == 4][["AGS", "geometry"]].rename(columns={"AGS": "AGS8"})

# ── Spatial join: WGS84 points → UTM32N to match VG250 ───────────────────────

gdf = gpd.GeoDataFrame(
    raw,
    geometry=gpd.points_from_xy(raw["lon"], raw["lat"]),
    crs="EPSG:4326",
).to_crs("EPSG:25832")

# AGS8: join to Gemeinden
gdf = gpd.sjoin(gdf, gem, how="left", predicate="within")
gdf = gdf[~gdf.index.duplicated(keep="first")].drop(columns="index_right")

# ── Assemble station-level dataset ────────────────────────────────────────────

stations = (
    pd.DataFrame(gdf.drop(columns=["geometry", "lat", "lon", "Breitengrad", "Längengrad"],
                           errors="ignore"))
    .assign(AGS8=lambda df: df["AGS8"].astype(str).str.zfill(8))
)

print(f"Total stations:       {len(stations)}")
print(f"Matched to AGS8:      {stations['AGS8'].notna().sum()}")
print(f"Unmatched AGS8:       {stations['AGS8'].isna().sum()}")
print(f"Year range:           {stations['year_commissioned'].min()} – {stations['year_commissioned'].max()}")

# ── Build AGS8 × year cumulative stock panel ──────────────────────────────────
#
# For each (AGS8, year): count all stations commissioned up to and including
# that year, and sum their charging-point capacity.  This is a stock measure
# (how many stations/points existed at year-end), not a flow.

s8 = stations.dropna(subset=["AGS8"]).copy()
year_min = max(int(s8["year_commissioned"].min()), 2005)
year_max = int(s8["year_commissioned"].max())

ags8_rows = []
for y in range(year_min, year_max + 1):
    stock = (
        s8[s8["year_commissioned"] <= y]
        .groupby("AGS8", as_index=False)
        .agg(
            n_ladestationen=("station_id",           "count"),
            n_ladepunkte   =("n_ladepunkte_station", "sum"),
        )
        .assign(year=y)
    )
    ags8_rows.append(stock)

ags8_panel = pd.concat(ags8_rows, ignore_index=True)
ags8_panel["AGS8"] = ags8_panel["AGS8"].str.zfill(8)
ags8_panel = ags8_panel.rename(columns={
    "n_ladestationen": "ev_stations",
    "n_ladepunkte":    "ev_chargepoints",
})

# ── Export ────────────────────────────────────────────────────────────────────

ags8_panel.to_csv(f"{DATA_DIR_INTERMEDIATE}/ladestationen/ladestationen_ags8_panel.csv", index=False)

print()
print(f"AGS8 panel shape:     {ags8_panel.shape}")
print()
print("ev_stations by year (cumulative stock, all Gemeinden):")
print(
    ags8_panel.groupby("year")["ev_stations"].sum()
    .rename("total_stations_stock")
    .to_string()
)
