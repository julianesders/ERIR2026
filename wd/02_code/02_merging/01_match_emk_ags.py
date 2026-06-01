import pandas as pd
import linktransformer as lt

BASE_PATH             = "/Users/julian/Documents/ERIR2026/wd"
DATA_DIR_INTERMEDIATE = f"{BASE_PATH}/01_data/02_intermediate"

# ── Load data ─────────────────────────────────────────────────────────────────
emk = pd.read_csv(f"{DATA_DIR_INTERMEDIATE}/emk/emkonzepte_cleaned.csv", delimiter=",")
anschriften = pd.read_csv(
    f"{DATA_DIR_INTERMEDIATE}/anschriftenverzeichnis/anschriftenverzeichnis_combined.csv",
    delimiter=",", dtype=str
)
emk["Bundesland"] = emk["bundesland"]

emk_gk = emk[emk["empfaenger"] == "Gebietskörperschaft"]

# ── Deduplicate before merge ──────────────────────────────────────────────────
key_matches = ["bundesland", "antragsteller"]
emk_gk_unique = emk_gk.drop_duplicates(subset=key_matches)
emk_gk_dupes  = emk_gk[emk_gk.duplicated(subset=key_matches, keep="first")]
print(f"Unique: {len(emk_gk_unique)}  |  Duplicates: {len(emk_gk_dupes)}")

# ── LinkTransformer merge ─────────────────────────────────────────────────────
matches = lt.merge_blocking(
    emk_gk_unique, anschriften,
    merge_type="m:1",
    left_on=["antragsteller"],
    right_on=["Gemeinde_Stadt"],
    model="T-Systems-onsite/cross-en-de-roberta-sentence-transformer",
    blocking_vars=["Bundesland"],
)

# ── Rejoin duplicates with matched columns from representative ────────────────
new_cols = [c for c in matches.columns if c not in emk_gk.columns]
lookup = matches.set_index(key_matches)[new_cols]
dupes_filled = emk_gk_dupes.join(lookup, on=key_matches)
matches_full = pd.concat([matches, dupes_filled], ignore_index=True)
print(f"Total rows after rejoin: {len(matches_full)}")

# ── Flag ambiguous matches (duplicate keys in right df) ───────────────────────
dupes_right = (
    anschriften[anschriften.duplicated(subset=["Bundesland", "Gemeinde_Stadt"], keep=False)]
    [["Bundesland", "Gemeinde_Stadt"]].drop_duplicates()
)
ambiguous = matches_full[
    matches_full.set_index(["bundesland", "Gemeinde_Stadt"]).index.isin(
        dupes_right.set_index(["Bundesland", "Gemeinde_Stadt"]).index
    )
]
print(f"Ambiguous matches: {len(ambiguous)}")

# ── Manual overrides ──────────────────────────────────────────────────────────
OVERRIDES: dict[str, dict] = {
    "03EMK255":  {"AGS5": "09184", "AGS8": None,         "Gemeinde_Stadt": "München",                                       "TKZ_name": "Landkreis"},
    "03EMK3058": {"AGS5": "01054", "AGS8": None,         "Gemeinde_Stadt": "Nordfriesland",                                 "TKZ_name": "Landkreis"},
    "03EMK280":  {"AGS5": "06634", "AGS8": "06634009",   "Gemeinde_Stadt": "Homberg (Efze), Reformationsstadt, Kreisstadt", "TKZ_name": "Stadt"},
    "03EMK3043": {"AGS5": "08116", "AGS8": None,         "Gemeinde_Stadt": "Esslingen",                                     "TKZ_name": "Landkreis"},
    "03EMK269":  {"AGS5": "08118", "AGS8": "08118011",   "Gemeinde_Stadt": "Ditzingen, Stadt",                              "TKZ_name": "Große Kreisstadt"},
    "03EMK206":  {"AGS5": "09178", "AGS8": "09178120",   "Gemeinde_Stadt": "Eching",                                        "TKZ_name": "Kreisangehörige Gemeinde"},
    "03EMK4100": {"AGS5": "08111", "AGS8": "08111000",   "Gemeinde_Stadt": "Stuttgart, Landeshauptstadt",                   "TKZ_name": "Stadtkreis"},
    "03EMK273":  {"AGS5": "09773", "AGS8": "09773182",   "Gemeinde_Stadt": "Wertingen, St",                                 "TKZ_name": "Stadt"},
    "03EMK216":  {"AGS5": "01053", "AGS8": None,         "Gemeinde_Stadt": "Herzogtum Lauenburg",                           "TKZ_name": "Kreis"},
    "03EMK104":  {"AGS5": "09475", "AGS8": None,         "Gemeinde_Stadt": "Hof",                                           "TKZ_name": "Landkreis"},
    "03EMK102":  {"AGS5": "05974", "AGS8": None,         "Gemeinde_Stadt": "Soest",                                         "TKZ_name": "Kreis"},
    "03EMK3092": {"AGS5": "08316", "AGS8": "08316043",   "Gemeinde_Stadt": "Teningen",                                      "TKZ_name": "Kreisangehörige Gemeinde"},
    "03EMK3017": {"AGS5": "08315", "AGS8": "08315047",   "Gemeinde_Stadt": "Gundelfingen",                                  "TKZ_name": "Kreisangehörige Gemeinde"},
    "03EMK5052": {"AGS5": "03452", "AGS8": None,         "Gemeinde_Stadt": "Aurich",                                        "TKZ_name": "Landkreis"},
    "03EMK5017": {"AGS5": "03457", "AGS8": "03457002",   "Gemeinde_Stadt": "Borkum, Stadt",                                 "TKZ_name": "Stadt"},
    "03EMK4213": {"AGS5": "09162", "AGS8": "09162000",   "Gemeinde_Stadt": "München, Landeshauptstadt",                     "TKZ_name": "Kreisfreie Stadt"},
    "03EMK5073": {"AGS5": "05515", "AGS8": "05515000",   "Gemeinde_Stadt": "Münster, Stadt",                                "TKZ_name": "Kreisfreie Stadt"},
    "03EMK265":  {"AGS5": "09776", "AGS8": "09776125",   "Gemeinde_Stadt": "Scheidegg, M.",                                 "TKZ_name": "Markt"},
    "03EMK4097": {"AGS5": "07313", "AGS8": "07313000",   "Gemeinde_Stadt": "Landau in der Pfalz, kreisfreie Stadt",         "TKZ_name": "Kreisfreie Stadt"},
    "03EMK3073": {"AGS5": "06531", "AGS8": None,         "Gemeinde_Stadt": "Gießen",                                        "TKZ_name": "Landkreis"},
    "03EMK034":  {"AGS5": "12060", "AGS8": "12060005",   "Gemeinde_Stadt": "Ahrensfelde",                                   "TKZ_name": "Kreisangehörige Gemeinde"},
    "03EMK3055": {"AGS5": "14522", "AGS8": "14522230",   "Gemeinde_Stadt": "Hainichen, Stadt",                              "TKZ_name": "Stadt"},
    "03EMK4089": {"AGS5": "03361", "AGS8": None,         "Gemeinde_Stadt": "Verden",                                        "TKZ_name": "Landkreis"},
}

known_fkz = set(matches_full["foerderkennzeichen"].unique())
for fkz in OVERRIDES:
    if fkz not in known_fkz:
        print(f"WARNING: {fkz} not found in matches_full")

for fkz, vals in OVERRIDES.items():
    mask = matches_full["foerderkennzeichen"] == fkz
    for col, val in vals.items():
        matches_full.loc[mask, col] = val

print("Overrides applied.")

# ── Export ────────────────────────────────────────────────────────────────────
id_cols    = ["foerderkennzeichen", "antragsteller_id", "antragsteller", "empfaenger", "bundesland"]
ags_cols   = ["Gemeinde_Stadt", "TKZ_name", "AGS8", "AGS5", "AGS3", "AGS2", "score"]
meta_cols  = ["raumtyp", "laufzeit_start", "laufzeit_end", "laufzeit_days", "gesamtmittel", "bundesmittel"]
tag_cols   = [c for c in emk_gk.columns if c.startswith("tag_")]
space_cols = [c for c in emk_gk.columns if c.startswith("space_")]
export_cols = id_cols + ags_cols + meta_cols + tag_cols + space_cols

(
    matches_full[export_cols]
    .sort_values("score", ascending=True)
    .to_csv(f"{DATA_DIR_INTERMEDIATE}/emk/emk_ags_matched.csv", index=False, encoding="utf-8")
)
print("Exported emk/emk_ags_matched.csv")

# ── Review export (non-overridden automatic matches, worst-first) ─────────────
review_cols = ["foerderkennzeichen", "antragsteller", "Gemeinde_Stadt", "TKZ_name",
               "AGS8", "AGS5", "bundesland", "score"]
(
    matches_full[review_cols]
    .sort_values("score", ascending=True)
    .to_excel(f"{DATA_DIR_INTERMEDIATE}/emk/emk_ags_review.xlsx", index=False)
)
print("Exported emk/emk_ags_review.xlsx")
