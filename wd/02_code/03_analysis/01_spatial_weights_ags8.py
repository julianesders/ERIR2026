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

TREAT_COLS = ["emk_active", "emk_absorbing"]


# ── Geometry ──────────────────────────────────────────────────────────────────

gem = gpd.read_file(GPKG_PATH, layer="vg250_gem")
gem = gem[gem["GF"] == 4][["AGS", "geometry"]].rename(columns={"AGS": "AGS8"})
gem["AGS8"] = gem["AGS8"].str.zfill(8)

# Restrict to units present in the panel (all 10,779 match)
panel = pd.read_csv(
    f"{DATA_DIR_FINAL}/emk_inkar_panel_ags8.csv",
    dtype={"AGS8": str},
)
panel_units = sorted(panel["AGS8"].unique())
gdf = gem[gem["AGS8"].isin(panel_units)].copy().reset_index(drop=True)
print(f"Geometry units loaded: {len(gdf)}")


# ── Queen contiguity weights ──────────────────────────────────────────────────

print("Building queen contiguity weights ...")
w = libpysal.weights.Queen.from_dataframe(gdf, ids="AGS8")
w.transform = "b"  # binary (0/1)
print(f"Queen weights built. Islands (no neighbors): {len(w.islands)}")


# ── Exclusive neighbor sets: 1st, 2nd, 3rd degree ────────────────────────────

units    = [str(u) for u in gdf["AGS8"].tolist()]
unit_idx = {u: i for i, u in enumerate(units)}

nbrs_1 = {u: set(str(n) for n in w.neighbors.get(u, [])) for u in units}

nbrs_2 = {}
for u in units:
    n2 = set()
    for n1 in nbrs_1[u]:
        n2.update(nbrs_1.get(n1, set()))
    n2 -= nbrs_1[u]
    n2.discard(u)
    nbrs_2[u] = n2

nbrs_3 = {}
for u in units:
    n3 = set()
    for n2 in nbrs_2[u]:
        n3.update(nbrs_1.get(n2, set()))
    n3 -= nbrs_1[u]
    n3 -= nbrs_2[u]
    n3.discard(u)
    nbrs_3[u] = n3

# Summary
for deg, nd in [(1, nbrs_1), (2, nbrs_2), (3, nbrs_3)]:
    counts = [len(v) for v in nd.values()]
    print(f"Degree {deg}: min={min(counts)}  mean={np.mean(counts):.1f}  "
          f"median={np.median(counts):.0f}  max={max(counts)}")


# ── Save neighbor structure ───────────────────────────────────────────────────

with open(f"{DATA_DIR_INTERMEDIATE}/spatial/queen_neighbors_ags8.json", "w") as f:
    json.dump({u: sorted(v) for u, v in nbrs_1.items()}, f)
print("Saved spatial/queen_neighbors_ags8.json (1st-degree neighbor lists)")


# ── Sparse adjacency matrices ─────────────────────────────────────────────────

def _to_sparse(nbrs: dict, unit_idx: dict, n: int) -> sp.csr_matrix:
    rows, cols = [], []
    for u, nbr_set in nbrs.items():
        i = unit_idx[u]
        for nb in nbr_set:
            j = unit_idx.get(nb)
            if j is not None:
                rows.append(i)
                cols.append(j)
    return sp.csr_matrix(
        (np.ones(len(rows), dtype=np.float32), (rows, cols)),
        shape=(n, n),
    )

n = len(units)
A = {1: _to_sparse(nbrs_1, unit_idx, n),
     2: _to_sparse(nbrs_2, unit_idx, n),
     3: _to_sparse(nbrs_3, unit_idx, n)}

# Time-invariant neighbor counts
n_nbrs = {deg: np.array(A[deg].sum(axis=1)).flatten().astype(int) for deg in [1, 2, 3]}


# ── Compute neighbor treatment indicators for each year ───────────────────────

panel_years = sorted(panel["year"].unique())
records = []

for year in panel_years:
    sub = panel[panel["year"] == year].set_index("AGS8")

    row = {"AGS8": units, "year": year}
    row["n_nbrs_1"] = n_nbrs[1]
    row["n_nbrs_2"] = n_nbrs[2]
    row["n_nbrs_3"] = n_nbrs[3]

    for tc in TREAT_COLS:
        t_vec = sub[tc].reindex(units).fillna(0).values.astype(np.float32)
        for deg in [1, 2, 3]:
            sums = np.array(A[deg] @ t_vec).flatten()
            row[f"{tc}_nbrs_{deg}"]     = sums.astype(int)
            row[f"{tc}_any_nbrs_{deg}"] = (sums > 0).astype(int)

    records.append(pd.DataFrame(row))

spatial = pd.concat(records, ignore_index=True)
spatial["year"] = spatial["year"].astype(int)

# ── Export ────────────────────────────────────────────────────────────────────

spatial.to_csv(f"{DATA_DIR_FINAL}/spatial_neighbors_ags8.csv", index=False)

print(f"\nSpatial neighbor panel shape: {spatial.shape}")
print(f"AGS8 units: {spatial['AGS8'].nunique()}")
print(f"Years:      {sorted(spatial['year'].unique().tolist())}")
print(f"\nNeighbor counts (avg across all units):")
print(f"  1st degree: {spatial['n_nbrs_1'].mean():.1f}")
print(f"  2nd degree: {spatial['n_nbrs_2'].mean():.1f}")
print(f"  3rd degree: {spatial['n_nbrs_3'].mean():.1f}")
print(f"\nTreatment spill-over (2021, 1st degree):")
y2021 = spatial[spatial["year"] == 2021]
for tc in TREAT_COLS:
    pct = (y2021[f"{tc}_any_nbrs_1"] == 1).mean() * 100
    print(f"  any {tc} neighbor: {pct:.1f}% of Gemeinden")
