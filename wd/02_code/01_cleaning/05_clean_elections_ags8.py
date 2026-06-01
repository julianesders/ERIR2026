import pandas as pd

BASE_PATH             = "/Users/julian/Documents/ERIR2026/wd"
DATA_DIR_RAW          = f"{BASE_PATH}/01_data/01_raw"
DATA_DIR_INTERMEDIATE = f"{BASE_PATH}/01_data/02_intermediate"

MAIN_PARTIES  = ["cdu_csu", "spd", "linke_pds", "gruene", "afd", "fdp"]
PARTY_COLS    = MAIN_PARTIES + ["other"]
CROSSWALK_PATH = f"{DATA_DIR_RAW}/crosswalks/ref-gemeinden-2020-2024.xlsx"

# ── AGS8 unit list ────────────────────────────────────────────────────────────

ags8_units = pd.read_csv(
    f"{DATA_DIR_INTERMEDIATE}/anschriftenverzeichnis/anschriftenverzeichnis_ags8.csv",
    dtype=str,
    usecols=["AGS8"],
).drop_duplicates()


# ── Crosswalk: 2021 → 2023 Gemeinde boundaries ───────────────────────────────
# Chain annual sheets 2021→2022 and 2022→2023.
# Population-proportional weight handles both mergers (weight=1 for all predecessors)
# and the rare splits (weight<1, proportional to population share).

def _build_crosswalk(path: str) -> pd.DataFrame:
    s21 = pd.read_excel(path, sheet_name="2021", header=0, dtype=str)
    cw_a = pd.DataFrame({
        "ags_from": s21.iloc[:, 0].str.strip().str.zfill(8),
        "ags_mid":  s21.iloc[:, 8].str.strip().str.zfill(8),
        "w_a":      pd.to_numeric(s21.iloc[:, 3], errors="coerce"),
    }).dropna(subset=["ags_from", "ags_mid"])

    s22 = pd.read_excel(path, sheet_name="2022", header=0, dtype=str)
    cw_b = pd.DataFrame({
        "ags_mid":  s22.iloc[:, 0].str.strip().str.zfill(8),
        "ags_to":   s22.iloc[:, 7].str.strip().str.zfill(8),
        "w_b":      pd.to_numeric(s22.iloc[:, 4], errors="coerce"),
    }).dropna(subset=["ags_mid", "ags_to"])

    cw = cw_a.merge(cw_b, on="ags_mid", how="left")
    cw["ags_to"] = cw["ags_to"].fillna(cw["ags_mid"])
    cw["w_b"]    = cw["w_b"].fillna(1.0)
    cw["weight"] = cw["w_a"] * cw["w_b"]
    return cw[["ags_from", "ags_to", "weight"]]

crosswalk = _build_crosswalk(CROSSWALK_PATH)
n_changes = (crosswalk["ags_from"] != crosswalk["ags_to"]).sum()
print(f"Crosswalk 2021→2023: {len(crosswalk)} entries, {n_changes} boundary changes")


# ── Helpers ───────────────────────────────────────────────────────────────────

def _compute_other(df: pd.DataFrame) -> pd.Series:
    """Other = 1 − sum(six main parties); NaN if all main parties are missing."""
    has_data = df[MAIN_PARTIES].notna().any(axis=1)
    known    = df[MAIN_PARTIES].fillna(0).sum(axis=1)
    return (1 - known).clip(lower=0).where(has_data)


def _harmonize(df: pd.DataFrame, cw: pd.DataFrame) -> pd.DataFrame:
    """
    Convert election DataFrame from 2021 to 2023 Gemeinde boundaries.

    df must have: AGS8, year, eligible_voters, number_voters, valid_votes,
                  + MAIN_PARTIES as vote shares.
    Returns same structure with AGS8 on 2023 boundaries and shares recomputed
    from aggregated raw counts. Count columns are dropped on return.
    """
    df = df.copy()

    # Back-calculate raw party vote counts (share × valid_votes)
    for col in MAIN_PARTIES:
        df[f"_n_{col}"] = df[col] * df["valid_votes"]

    # Attach crosswalk; unmatched units map to themselves with weight 1
    df = df.merge(cw.rename(columns={"ags_from": "AGS8"}), on="AGS8", how="left")
    df["ags_to"] = df["ags_to"].fillna(df["AGS8"])
    df["weight"] = df["weight"].fillna(1.0)

    # Apply weight to all count columns
    cnt_cols = ["eligible_voters", "number_voters", "valid_votes"] + [f"_n_{c}" for c in MAIN_PARTIES]
    for col in cnt_cols:
        df[col] = df[col] * df["weight"]

    # Aggregate to 2023 successor Gemeinden
    df = (
        df.groupby(["ags_to", "year"], as_index=False)[cnt_cols]
        .sum()
        .rename(columns={"ags_to": "AGS8"})
    )

    # Recompute shares from aggregated counts
    df["turnout"] = df["number_voters"] / df["eligible_voters"]
    for col in MAIN_PARTIES:
        df[col] = df[f"_n_{col}"] / df["valid_votes"]
        df.drop(columns=[f"_n_{col}"], inplace=True)

    return df.drop(columns=["eligible_voters", "number_voters", "valid_votes"])


def _forward_fill(elections: pd.DataFrame, prefix: str) -> pd.DataFrame:
    """
    Merge election results onto a contiguous AGS8 × year spine and LOCF.
    elections must have: AGS8, year, turnout, + PARTY_COLS.
    Spine runs from the earliest election year through 5 years after the last,
    so the final LOCF value is available well beyond the most recent election.
    """
    elections = elections.rename(
        columns={c: f"{prefix}_{c}" for c in ["turnout"] + PARTY_COLS}
    )
    p_cols = [f"{prefix}_{c}" for c in ["turnout"] + PARTY_COLS]

    min_yr = max(int(elections["year"].dropna().min()), 2005)
    max_yr = int(elections["year"].dropna().max()) + 5
    spine = ags8_units.assign(_k=1).merge(
        pd.DataFrame({"year": range(min_yr, max_yr + 1), "_k": 1}), on="_k"
    ).drop(columns="_k")

    merged = spine.merge(elections, on=["AGS8", "year"], how="left")
    merged = merged.sort_values(["AGS8", "year"])
    merged[p_cols] = merged.groupby("AGS8")[p_cols].transform(lambda s: s.ffill())

    return merged.reset_index(drop=True)


# ── Federal elections ─────────────────────────────────────────────────────────

fed = pd.read_csv(
    f"{DATA_DIR_RAW}/election_GERDA/federal_muni_harm_21.csv",
    usecols=["ags", "election_year", "eligible_voters", "number_voters", "valid_votes",
             "cdu", "csu", "spd", "linke_pds", "gruene", "afd", "fdp"],
    dtype=str,
    low_memory=False,
)
fed = fed.rename(columns={"ags": "AGS8", "election_year": "year"})
fed["AGS8"] = fed["AGS8"].str.zfill(8)
fed["year"] = pd.to_numeric(fed["year"], errors="coerce")

for col in ["eligible_voters", "number_voters", "valid_votes",
            "cdu", "csu", "spd", "linke_pds", "gruene", "afd", "fdp"]:
    fed[col] = pd.to_numeric(fed[col], errors="coerce")

# CDU and CSU run in mutually exclusive states
fed["cdu_csu"] = fed["cdu"].fillna(0) + fed["csu"].fillna(0)
fed.loc[fed["cdu"].isna() & fed["csu"].isna(), "cdu_csu"] = float("nan")

fed = fed[["AGS8", "year", "eligible_voters", "number_voters", "valid_votes"] + MAIN_PARTIES]
fed = _harmonize(fed, crosswalk)
fed["other"] = _compute_other(fed)
fed = fed[["AGS8", "year", "turnout"] + PARTY_COLS]
fed_panel = _forward_fill(fed, "fed")

print(f"Federal: {fed['year'].nunique()} election years, {fed['AGS8'].nunique()} units")
print(f"Federal panel NaN turnout: {fed_panel['fed_turnout'].isna().sum()} of {len(fed_panel)}")


# ── State elections ───────────────────────────────────────────────────────────

state = pd.read_csv(
    f"{DATA_DIR_RAW}/election_GERDA/state_harm_21.csv",
    usecols=["ags", "election_year", "eligible_voters", "number_voters", "valid_votes",
             "cdu_csu", "spd", "linke_pds", "gruene", "afd", "fdp",
             "flag_unsuccessful_naive_merge", "flag_harm_turnout_above_1"],
    dtype=str,
    low_memory=False,
)
state = state[
    (state["flag_unsuccessful_naive_merge"] != "1") &
    (state["flag_harm_turnout_above_1"]      != "1")
].copy()
state = state.rename(columns={"ags": "AGS8", "election_year": "year"})
state["AGS8"] = state["AGS8"].str.zfill(8)
state["year"] = pd.to_numeric(state["year"], errors="coerce")

for col in ["eligible_voters", "number_voters", "valid_votes",
            "cdu_csu", "spd", "linke_pds", "gruene", "afd", "fdp"]:
    state[col] = pd.to_numeric(state[col], errors="coerce")

state = state[["AGS8", "year", "eligible_voters", "number_voters", "valid_votes"] + MAIN_PARTIES]
state = _harmonize(state, crosswalk)
state["other"] = _compute_other(state)
state = state[["AGS8", "year", "turnout"] + PARTY_COLS]
state_panel = _forward_fill(state, "state")

print(f"State: {state['year'].nunique()} election years, {state['AGS8'].nunique()} units")
print(f"State panel NaN turnout: {state_panel['state_turnout'].isna().sum()} of {len(state_panel)}")


# ── Municipal elections ───────────────────────────────────────────────────────

muni = pd.read_csv(
    f"{DATA_DIR_RAW}/election_GERDA/municipal_harm.csv",
    usecols=["ags", "election_year", "eligible_voters", "number_voters", "valid_votes",
             "cdu_csu", "spd", "linke_pds", "gruene", "afd", "fdp",
             "flag_unsuccessful_naive_merge", "election_type"],
    dtype=str,
    low_memory=False,
)
# Exclude failed geo-harmonizations and Berlin's Abgeordnetenhauswahl
muni = muni[
    (muni["flag_unsuccessful_naive_merge"] != "1") &
    (muni["election_type"] != "Abgeordnetenhauswahl (Zweitstimmen)")
].copy()
muni = muni.rename(columns={"ags": "AGS8", "election_year": "year"})
muni["AGS8"] = muni["AGS8"].str.zfill(8)
muni["year"] = pd.to_numeric(muni["year"], errors="coerce")

for col in ["eligible_voters", "number_voters", "valid_votes",
            "cdu_csu", "spd", "linke_pds", "gruene", "afd", "fdp"]:
    muni[col] = pd.to_numeric(muni[col], errors="coerce")

muni = muni[["AGS8", "year", "eligible_voters", "number_voters", "valid_votes"] + MAIN_PARTIES]
muni = _harmonize(muni, crosswalk)
muni["other"] = _compute_other(muni)
muni = muni[["AGS8", "year", "turnout"] + PARTY_COLS]
muni_panel = _forward_fill(muni, "muni")

print(f"Municipal: {muni['year'].nunique()} election years, {muni['AGS8'].nunique()} units")
print(f"Municipal panel NaN turnout: {muni_panel['muni_turnout'].isna().sum()} of {len(muni_panel)}")


# ── Combine and export ────────────────────────────────────────────────────────

panel = (
    fed_panel
    .merge(state_panel, on=["AGS8", "year"], how="left")
    .merge(muni_panel,  on=["AGS8", "year"], how="left")
    .sort_values(["AGS8", "year"])
    .reset_index(drop=True)
)

panel.to_csv(f"{DATA_DIR_INTERMEDIATE}/elections/elections_ags8_panel.csv", index=False)

print(f"\nFinal panel shape:  {panel.shape}")
print(f"AGS8 units:         {panel['AGS8'].nunique()}")
print(f"Years:              {sorted(panel['year'].unique().tolist())}")
print(f"\nNaN share by column:")
for col in ["fed_turnout", "state_turnout", "muni_turnout"]:
    pct = panel[col].isna().mean() * 100
    print(f"  {col:25s}: {pct:.1f}%")
