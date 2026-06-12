"""
Build two parallel spatial-neighbour structures and per-year neighbour-
treatment indicators:

  GRANULAR  — Queen contiguity at the AGS8 (Gemeinde) level.
              Treatment definition: DIRECT only (own AGS8 has an EMK
              project; broadcast units do NOT count).
              Output columns: direct_treated_(any_)nbrs_gem_{1,2,3}

  AGGREGATED — Queen contiguity at the AGS5 (Kreis) level: two Kreise are
              neighbours if any of their constituent Gemeinden share a
              boundary. Treatment definition: BROAD (any AGS8 in the Kreis
              under direct ∪ Kreis-broadcast coverage).
              Output columns: broad_treated_(any_)nbrs_kreis

Outputs (01_data/03_final/spatial_neighbors_ags8.csv): one row per
(AGS8, year) with both the granular and the aggregated indicators.
"""
import json

import geopandas as gpd
import libpysal
import numpy as np
import pandas as pd
import scipy.sparse as sp

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parents[2]))
from config import BASE_PATH

DATA_DIR_INTERMEDIATE = f"{BASE_PATH}/01_data/02_intermediate"
DATA_DIR_FINAL        = f"{BASE_PATH}/01_data/03_final"
GPKG_PATH = (
    f"{BASE_PATH}/01_data/04_shapefiles/"
    "vg250_01-01.utm32s.gpkg.ebenen/vg250_ebenen_0101/DE_VG250.gpkg"
)

# ── Load panel + geometries ───────────────────────────────────────────────────

panel = pd.read_csv(
    f"{DATA_DIR_FINAL}/emk_inkar_panel_ags8.csv",
    dtype={"AGS8": str, "AGS5": str},
)
panel["AGS8"] = panel["AGS8"].str.zfill(8)
panel["AGS5"] = panel["AGS5"].str.zfill(5)

gem = gpd.read_file(GPKG_PATH, layer="vg250_gem")
gem = gem[gem["GF"] == 4][["AGS", "geometry"]].rename(columns={"AGS": "AGS8"})
gem["AGS8"] = gem["AGS8"].str.zfill(8)

panel_units = sorted(panel["AGS8"].unique())
gdf = gem[gem["AGS8"].isin(panel_units)].copy().reset_index(drop=True)

# AGS5 is the first 5 chars of AGS8
gdf["AGS5"] = gdf["AGS8"].str[:5]
print(f"Geometry units loaded: {len(gdf)} AGS8 in {gdf['AGS5'].nunique()} AGS5")

# ── Direct-treatment indicator at the AGS8 level (granular) ───────────────────
# `first_treat_direct` is constant within AGS8 (NA for never-direct).
# direct_active(t) = 1{first_treat_direct <= t} (absorbing).

dir_treat = (
    panel[["AGS8", "year", "first_treat_direct"]]
    .assign(
        direct_active=lambda d: (
            d["first_treat_direct"].notna()
            & (d["first_treat_direct"] <= d["year"])
        ).astype(int)
    )
    [["AGS8", "year", "direct_active"]]
)

# ── Broad-treatment indicator at the AGS5 level (aggregated) ──────────────────
# Kreis is "treated" in year t if ANY AGS8 in that Kreis is under broad
# coverage (direct ∪ broadcast) by year t.

broad_treat_ags8 = (
    panel[["AGS8", "AGS5", "year", "first_treat_broad"]]
    .assign(
        broad_active=lambda d: (
            d["first_treat_broad"].notna()
            & (d["first_treat_broad"] <= d["year"])
        ).astype(int)
    )
)
broad_treat_ags5 = (
    broad_treat_ags8
    .groupby(["AGS5", "year"], as_index=False)["broad_active"].max()
    .rename(columns={"broad_active": "broad_active_kreis"})
)


# ──────────────────────────────────────────────────────────────────────────────
# 1) GRANULAR: AGS8 queen contiguity + direct-treatment neighbours
# ──────────────────────────────────────────────────────────────────────────────

print("\nBuilding Gemeinde-level queen contiguity ...")
w_gem = libpysal.weights.Queen.from_dataframe(gdf, ids="AGS8")
w_gem.transform = "b"
print(f"  Gemeinde queen weights. Islands: {len(w_gem.islands)}")

units    = [str(u) for u in gdf["AGS8"].tolist()]
unit_idx = {u: i for i, u in enumerate(units)}

nbrs_1 = {u: set(str(n) for n in w_gem.neighbors.get(u, [])) for u in units}

nbrs_2, nbrs_3 = {}, {}
for u in units:
    n2 = set()
    for n1 in nbrs_1[u]:
        n2.update(nbrs_1.get(n1, set()))
    n2 -= nbrs_1[u]
    n2.discard(u)
    nbrs_2[u] = n2
for u in units:
    n3 = set()
    for nn in nbrs_2[u]:
        n3.update(nbrs_1.get(nn, set()))
    n3 -= nbrs_1[u]
    n3 -= nbrs_2[u]
    n3.discard(u)
    nbrs_3[u] = n3

for deg, nd in [(1, nbrs_1), (2, nbrs_2), (3, nbrs_3)]:
    counts = [len(v) for v in nd.values()]
    print(f"  Gemeinde ring {deg}: min={min(counts)}  mean={np.mean(counts):.1f}  "
          f"median={np.median(counts):.0f}  max={max(counts)}")

with open(f"{DATA_DIR_INTERMEDIATE}/spatial/queen_neighbors_ags8.json", "w") as f:
    json.dump({u: sorted(v) for u, v in nbrs_1.items()}, f)


def _to_sparse(nbrs: dict, idx: dict, n: int) -> sp.csr_matrix:
    rows, cols = [], []
    for u, nbr_set in nbrs.items():
        i = idx[u]
        for nb in nbr_set:
            j = idx.get(nb)
            if j is not None:
                rows.append(i)
                cols.append(j)
    return sp.csr_matrix(
        (np.ones(len(rows), dtype=np.float32), (rows, cols)),
        shape=(n, n),
    )


n_gem = len(units)
A_gem = {d: _to_sparse(nbrs_d, unit_idx, n_gem)
         for d, nbrs_d in [(1, nbrs_1), (2, nbrs_2), (3, nbrs_3)]}
n_nbrs_gem = {d: np.array(A_gem[d].sum(axis=1)).flatten().astype(int)
              for d in (1, 2, 3)}


# ──────────────────────────────────────────────────────────────────────────────
# 2) AGGREGATED: AGS5 queen contiguity + broad-treatment neighbours
# ──────────────────────────────────────────────────────────────────────────────

print("\nDissolving Gemeinden -> Kreise and building Kreis queen contiguity ...")
gdf5 = gdf.dissolve(by="AGS5", as_index=False)[["AGS5", "geometry"]]
print(f"  AGS5 polygons: {len(gdf5)}")

w_kr = libpysal.weights.Queen.from_dataframe(gdf5, ids="AGS5")
w_kr.transform = "b"
print(f"  Kreis queen weights. Islands: {len(w_kr.islands)}")

kr_units    = [str(u) for u in gdf5["AGS5"].tolist()]
kr_unit_idx = {u: i for i, u in enumerate(kr_units)}
kr_nbrs_1   = {u: set(str(n) for n in w_kr.neighbors.get(u, [])) for u in kr_units}

counts = [len(v) for v in kr_nbrs_1.values()]
print(f"  Kreis ring 1: min={min(counts)}  mean={np.mean(counts):.1f}  "
      f"median={np.median(counts):.0f}  max={max(counts)}")

n_kr = len(kr_units)
A_kr = _to_sparse(kr_nbrs_1, kr_unit_idx, n_kr)
n_nbrs_kreis_by_ags5 = pd.Series(
    np.array(A_kr.sum(axis=1)).flatten().astype(int),
    index=kr_units,
    name="n_nbrs_kreis",
)


# ──────────────────────────────────────────────────────────────────────────────
# 3) Per-year neighbour indicators
# ──────────────────────────────────────────────────────────────────────────────

panel_years = sorted(panel["year"].unique())
records = []

# Wide pivots for fast year-by-year lookup
dir_wide   = dir_treat.pivot(index="AGS8", columns="year",
                              values="direct_active").fillna(0)
broad_wide = (
    broad_treat_ags5.pivot(index="AGS5", columns="year",
                            values="broad_active_kreis").fillna(0)
)

for year in panel_years:
    # Granular: direct treatment of Gemeinde neighbours
    if year in dir_wide.columns:
        t_gem = dir_wide.reindex(units, fill_value=0)[year].values.astype(np.float32)
    else:
        t_gem = np.zeros(n_gem, dtype=np.float32)

    # Aggregated: broad treatment of Kreis neighbours
    if year in broad_wide.columns:
        t_kr = broad_wide.reindex(kr_units, fill_value=0)[year].values.astype(np.float32)
    else:
        t_kr = np.zeros(n_kr, dtype=np.float32)

    row = {
        "AGS8":         units,
        "year":         year,
        "n_nbrs_gem_1": n_nbrs_gem[1],
        "n_nbrs_gem_2": n_nbrs_gem[2],
        "n_nbrs_gem_3": n_nbrs_gem[3],
    }

    for deg in (1, 2, 3):
        sums = np.array(A_gem[deg] @ t_gem).flatten()
        row[f"direct_treated_nbrs_gem_{deg}"]     = sums.astype(int)
        row[f"direct_treated_any_nbrs_gem_{deg}"] = (sums > 0).astype(int)

    # Map Kreis-level indicator back to AGS8 (broadcast via prefix)
    sums_kr_per_ags5 = np.array(A_kr @ t_kr).flatten()
    kr_count = pd.Series(sums_kr_per_ags5, index=kr_units).astype(int)
    ags5_of_units = pd.Series([u[:5] for u in units])
    row["broad_treated_nbrs_kreis"]     = ags5_of_units.map(kr_count).fillna(0).astype(int).values
    row["broad_treated_any_nbrs_kreis"] = (row["broad_treated_nbrs_kreis"] > 0).astype(int)
    row["n_nbrs_kreis"]                 = ags5_of_units.map(n_nbrs_kreis_by_ags5).fillna(0).astype(int).values

    records.append(pd.DataFrame(row))

spatial = pd.concat(records, ignore_index=True)
spatial["year"] = spatial["year"].astype(int)

# ── Export ────────────────────────────────────────────────────────────────────

spatial.to_csv(f"{DATA_DIR_FINAL}/spatial_neighbors_ags8.csv", index=False)

print(f"\nSpatial neighbour panel shape: {spatial.shape}")
print(f"AGS8 units: {spatial['AGS8'].nunique()}")
print(f"Years:      {sorted(spatial['year'].unique().tolist())}")

print("\nAverage neighbour-set sizes (time-invariant):")
print(f"  Gemeinde ring 1: {spatial['n_nbrs_gem_1'].mean():.1f}")
print(f"  Gemeinde ring 2: {spatial['n_nbrs_gem_2'].mean():.1f}")
print(f"  Gemeinde ring 3: {spatial['n_nbrs_gem_3'].mean():.1f}")
print(f"  Kreis           : {spatial['n_nbrs_kreis'].mean():.1f}")

print("\nNeighbour-treatment coverage (2021):")
y = spatial[spatial["year"] == 2021]
print(f"  direct Gemeinde nbr (ring 1) treated: "
      f"{100 * (y['direct_treated_any_nbrs_gem_1'] == 1).mean():.1f}% of AGS8")
print(f"  any broad Kreis nbr treated         : "
      f"{100 * (y['broad_treated_any_nbrs_kreis'] == 1).mean():.1f}% of AGS8")
