#!/usr/bin/env python3
"""
Scrape ADAC car catalog for list prices, production year windows, and fuel types.

Input:  kba_neuz_model_panel.csv  (produced by 08_clean_kba_ags8.py)
Output: adac_prices_raw.csv

For each unique (fab_text, mod_text) pair in the KBA Neuzulassungen data
(Sonstige excluded at both brand and model level), the script fetches the
corresponding ADAC model-family page:
  /marken-modelle/{brand-slug}/{model-slug}/
and parses all variant rows: production year window, fuel type, list price.

Multiple KBA pairs that normalise to the same ADAC slug are fetched only once
to avoid redundant requests. Failed slugs are printed at the end with suggested
overrides so BRAND_SLUG_OVERRIDES / MODEL_SLUG_OVERRIDES can be filled in.

Runtime: ~1–2 s per unique ADAC slug (REQUEST_DELAY + network). Expect several
hours for a full KBA model list; run overnight or reduce REQUEST_DELAY at your
own risk.
"""

import re
import time
import requests
import pandas as pd
from collections import defaultdict
from lxml import html as lxml_html

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parents[2]))
from config import BASE_PATH

INTERMEDIATE  = f"{BASE_PATH}/01_data/02_intermediate/kba"
ADAC_BASE     = "https://www.adac.de/rund-ums-fahrzeug/autokatalog/marken-modelle"
REQUEST_DELAY = 1.0   # seconds between requests

SONSTIGE_FAB_CODE = "000"
SONSTIGE_MOD_CODE = "0000"

# ── Fuel type mapping: ADAC Kraftstoff label → KBA energiequelle code ─────────

FUEL_MAP: dict[str, str] = {
    "Super":                    "01",
    "Super Plus":               "01",
    "Benzin":                   "01",
    "Diesel":                   "02",
    "Strom":                    "04",
    "Autogas (LPG)":            "03",
    "Erdgas (CNG)":             "03",
    "Mildhybrid (Benzin)":      "05",
    "Mildhybrid (Diesel)":      "05",
    "Plug-in-Hybrid":           "06",
    "Plug-in-Hybrid (Benzin)":  "06",
    "Plug-in-Hybrid (Diesel)":  "06",
    "Wasserstoff":              "07",
}

# ── Slug overrides ─────────────────────────────────────────────────────────────
# KBA fab_text is capped at 15 chars, so truncated names don't normalise to the
# right ADAC slug automatically. Add model overrides as (fab_text, mod_text) → slug.

BRAND_SLUG_OVERRIDES: dict[str, str] = {
    "MERCEDES-B":    "mercedes-benz",
    "ROLLS-ROYCE":   "rolls-royce",
    "ALFA ROMEO":    "alfa-romeo",
    "ASTON MARTI":   "aston-martin",
    "LAND ROVER":    "land-rover",
    "GREAT WALL":    "great-wall",
    "DS AUTOMOBIL":  "ds-automobiles",
}

MODEL_SLUG_OVERRIDES: dict[tuple[str, str], str] = {
    # Example: ("MERCEDES-B", "C-KLASSE"): "c-klasse",
}


# ── Slug normalisation ─────────────────────────────────────────────────────────

def _to_slug(text: str) -> str:
    s = text.lower().strip()
    s = s.replace("ä", "ae").replace("ö", "oe").replace("ü", "ue").replace("ß", "ss")
    s = re.sub(r"[^a-z0-9]+", "-", s).strip("-")
    return s


def _brand_slug(fab_text: str) -> str:
    return BRAND_SLUG_OVERRIDES.get(fab_text, _to_slug(fab_text))


def _model_slug(fab_text: str, mod_text: str) -> str:
    return MODEL_SLUG_OVERRIDES.get((fab_text, mod_text), _to_slug(mod_text))


# ── HTTP / HTML helpers ────────────────────────────────────────────────────────

_SESSION = requests.Session()
_SESSION.headers.update({"User-Agent": "Mozilla/5.0 (compatible; research-scraper)"})


def _fetch(url: str) -> lxml_html.HtmlElement | None:
    try:
        resp = _SESSION.get(url, timeout=30)
        if resp.status_code == 200:
            return lxml_html.fromstring(resp.content)
    except requests.RequestException as exc:
        print(f"    request error: {exc}")
    return None


def _max_page(tree: lxml_html.HtmlElement) -> int:
    hrefs = tree.xpath('//a[contains(@href, "pageNumber=")]/@href')
    nums = [
        int(m.group(1))
        for href in hrefs
        if (m := re.search(r"pageNumber=(\d+)", href))
    ]
    return max(nums, default=1)


def _parse_year(token: str) -> int | None:
    m = re.fullmatch(r"\d{2}/(\d{2})", token.strip())
    return (2000 + int(m.group(1))) if m else None


def _parse_price(s: str) -> float | None:
    s = s.replace("€", "").replace("\xa0", "").strip()
    if not s or s == "-":
        return None
    try:
        return float(s.replace(".", "").replace(",", "."))
    except ValueError:
        return None


def _parse_rows(tree: lxml_html.HtmlElement) -> list[dict]:
    rows_out = []
    for row in tree.xpath('//tr[@data-testid="carpages:generation:model:row"]'):
        cells = {td.get("data-th"): td for td in row.xpath(".//td[@data-th]")}

        fahrzeug = cells.get("Fahrzeug")
        if fahrzeug is None:
            continue

        model_raw  = " ".join(fahrzeug.text_content().split())
        kraftstoff = cells["Kraftstoff"].text_content().strip() if "Kraftstoff" in cells else ""
        preis_raw  = cells["Listenpreis"].text_content().strip() if "Listenpreis" in cells else ""

        year_from = year_to = None
        yr = re.search(
            r"\((\d{2}/\d{2})\s*-\s*(\d{2}/\d{2}|[Hh]eute|[Pp]resent|k\.A\.)\)",
            model_raw,
        )
        if yr:
            year_from = _parse_year(yr.group(1))
            year_to   = _parse_year(yr.group(2))   # None when "heute"

        model_clean = re.sub(r"\s*\(\d{2}/\d{2}.*?\)\s*$", "", model_raw).strip()

        energiequelle = FUEL_MAP.get(kraftstoff)

        rows_out.append({
            "model_raw":      model_raw,
            "model_clean":    model_clean,
            "fuel_type_adac": kraftstoff,
            "energiequelle":  energiequelle,
            "year_from":      year_from,
            "year_to":        year_to,
            "list_price_eur": _parse_price(preis_raw),
        })
    return rows_out


def _scrape_slug(brand_slug: str, model_slug: str) -> list[dict]:
    """Fetch all paginated rows for one (brand, model) ADAC page."""
    base = f"{ADAC_BASE}/{brand_slug}/{model_slug}/?sort=ALPHABETIC_ASC&type=klassik-autosuche"

    tree = _fetch(f"{base}&pageNumber=1")
    if tree is None:
        return []

    rows = _parse_rows(tree)
    if not rows:
        return []

    for page in range(2, _max_page(tree) + 1):
        time.sleep(REQUEST_DELAY)
        t = _fetch(f"{base}&pageNumber={page}")
        if t is not None:
            rows.extend(_parse_rows(t))

    return rows


# ── Load unique (brand, model) pairs ──────────────────────────────────────────

neuz_model = pd.read_csv(
    f"{INTERMEDIATE}/kba_neuz_model_panel.csv",
    dtype={"AGS8": str, "fab_code": str, "mod_code": str},
)

pairs = (
    neuz_model[
        (neuz_model["fab_code"] != SONSTIGE_FAB_CODE) &
        (neuz_model["mod_code"] != SONSTIGE_MOD_CODE)
    ][["fab_text", "mod_text"]]
    .dropna()
    .drop_duplicates()
    .sort_values(["fab_text", "mod_text"])
    .reset_index(drop=True)
)
print(f"Unique (brand, model) pairs from KBA data: {len(pairs):,}")

# Group KBA pairs by their resolved ADAC slug — avoids redundant requests when
# e.g. "GOLF" and "GOLF GTI" both normalise to slug "golf".
slug_to_kba: dict[tuple[str, str], list[tuple[str, str]]] = defaultdict(list)
for _, row in pairs.iterrows():
    fab, mod = row["fab_text"], row["mod_text"]
    slug_to_kba[(_brand_slug(fab), _model_slug(fab, mod))].append((fab, mod))

print(f"Unique ADAC slugs to fetch:                {len(slug_to_kba):,}\n")


# ── Scrape ─────────────────────────────────────────────────────────────────────

all_dfs: list[pd.DataFrame] = []
failed:  list[tuple[str, str]] = []   # (brand_slug, model_slug) pairs

for (bslug, mslug), kba_pairs in slug_to_kba.items():
    rows = _scrape_slug(bslug, mslug)
    time.sleep(REQUEST_DELAY)

    if not rows:
        failed.append((bslug, mslug))
        kba_labels = ", ".join(f"{f}/{m}" for f, m in kba_pairs)
        print(f"  [MISS] {bslug}/{mslug}  ← KBA: {kba_labels}")
        continue

    df = pd.DataFrame(rows)
    # Tag every row with all KBA (fab_text, mod_text) pairs that map here
    for fab, mod in kba_pairs:
        tagged = df.copy()
        tagged.insert(0, "fab_text_kba", fab)
        tagged.insert(1, "mod_text_kba", mod)
        tagged.insert(2, "adac_brand_slug", bslug)
        tagged.insert(3, "adac_model_slug", mslug)
        all_dfs.append(tagged)

    print(f"  [OK]   {bslug}/{mslug}  {len(rows)} entries  "
          f"(KBA pairs: {len(kba_pairs)})", flush=True)


# ── Export ─────────────────────────────────────────────────────────────────────

adac_raw = pd.concat(all_dfs, ignore_index=True) if all_dfs else pd.DataFrame()
adac_raw.to_csv(f"{INTERMEDIATE}/adac_prices_raw.csv", index=False, encoding="utf-8-sig")

n_scraped = adac_raw["adac_brand_slug"].apply(lambda x: f"{x}").nunique() if not adac_raw.empty else 0
print(f"\nScraped {len(adac_raw):,} rows for "
      f"{adac_raw[['fab_text_kba','mod_text_kba']].drop_duplicates().shape[0]:,} KBA pairs")
print(f"Price coverage: "
      f"{adac_raw['list_price_eur'].notna().mean()*100:.1f}% of entries have a price")

unmapped = (
    adac_raw[adac_raw["energiequelle"].isna() & adac_raw["fuel_type_adac"].notna()]
    ["fuel_type_adac"].value_counts()
)
if not unmapped.empty:
    print(f"\nUnmapped Kraftstoff values — add to FUEL_MAP:")
    print(unmapped.to_string())

if failed:
    print(f"\nFailed slugs ({len(failed)}) — add to BRAND_SLUG_OVERRIDES / MODEL_SLUG_OVERRIDES:")
    for bslug, mslug in sorted(failed):
        kba_labels = ", ".join(
            f"'{f}'/'{m}'" for f, m in slug_to_kba[(bslug, mslug)]
        )
        print(f"  {bslug}/{mslug}  ← {kba_labels}")
