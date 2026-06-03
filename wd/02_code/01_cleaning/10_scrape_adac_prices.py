"""
Scrape the full ADAC car catalogue for list prices, fuel types, and production
year windows.  Results are written page-by-page to adac_prices_raw.csv so the
run can be interrupted and resumed without losing progress.

Inputs (optional — used only for coverage reporting at the end):
  kba_neuz_model_panel.csv   (produced by 08_clean_kba_ags8.py)

Output:
  adac_prices_raw.csv

URL strategy
------------
ADAC uses opaque brand IDs in its URLs (not human-readable slugs), and the
brand-filter GET parameter does not filter server-side.  The only reliable
endpoint is the alphabetical full-catalogue listing:
  /autosuche/?sort=ALPHABETIC_ASC&type=klassik-autosuche&pageNumber=N
Pagination continues until a page returns no car rows (~8 000 pages total).
At REQUEST_DELAY=1.0 s this takes roughly 2–2.5 h; run overnight on the server.
The checkpoint file lets you resume an interrupted run.
"""

import re
import time
import requests
import pandas as pd
from pathlib import Path
from lxml import html as lxml_html

import sys
sys.path.insert(0, str(Path(__file__).parents[2]))
from config import BASE_PATH

INTERMEDIATE  = f"{BASE_PATH}/01_data/02_intermediate/kba"
OUT_CSV       = f"{INTERMEDIATE}/adac_prices_raw.csv"
CHECKPOINT    = f"{INTERMEDIATE}/adac_scrape_checkpoint.txt"

ADAC_SEARCH   = (
    "https://www.adac.de/rund-ums-fahrzeug/autokatalog/marken-modelle"
    "/autosuche/?sort=ALPHABETIC_ASC&type=klassik-autosuche"
)
REQUEST_DELAY = 1.0   # seconds between requests — be polite

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


# ── Helpers ────────────────────────────────────────────────────────────────────

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
            year_to   = _parse_year(yr.group(2))

        model_clean = re.sub(r"\s*\(\d{2}/\d{2}.*?\)\s*$", "", model_raw).strip()

        rows_out.append({
            "model_raw":      model_raw,
            "model_clean":    model_clean,
            "fuel_type_adac": kraftstoff,
            "energiequelle":  FUEL_MAP.get(kraftstoff),
            "year_from":      year_from,
            "year_to":        year_to,
            "list_price_eur": _parse_price(preis_raw),
        })
    return rows_out


# ── Resumable scrape ───────────────────────────────────────────────────────────

def scrape_catalogue() -> None:
    start_page = 1
    if Path(CHECKPOINT).exists():
        start_page = int(Path(CHECKPOINT).read_text().strip()) + 1
        print(f"Resuming from page {start_page}  "
              f"(delete {CHECKPOINT} to restart from scratch)")

    mode   = "a" if start_page > 1 else "w"
    header = (start_page == 1)
    t0     = time.time()
    n_rows = 0

    print(f"Scraping ADAC catalogue starting at page {start_page} ...")
    print(f"Output: {OUT_CSV}\n")

    page = start_page
    while True:
        url  = f"{ADAC_SEARCH}&pageNumber={page}"
        tree = _fetch(url)

        if tree is None:
            print(f"Page {page}: fetch failed — stopping.")
            break

        rows = _parse_rows(tree)
        if not rows:
            print(f"Page {page}: no rows — catalogue exhausted after page {page - 1}.")
            break

        df = pd.DataFrame(rows)
        df.to_csv(OUT_CSV, mode=mode, header=header, index=False, encoding="utf-8-sig")
        mode   = "a"
        header = False
        n_rows += len(rows)

        Path(CHECKPOINT).write_text(str(page))

        if page % 100 == 0:
            elapsed = time.time() - t0
            rate    = (page - start_page + 1) / elapsed
            print(f"  page {page:>5}  |  {n_rows:>8,} rows  |  "
                  f"{rate:.1f} pages/s  |  {elapsed/60:.0f} min elapsed",
                  flush=True)

        page += 1
        time.sleep(REQUEST_DELAY)

    if Path(CHECKPOINT).exists():
        Path(CHECKPOINT).unlink()

    print(f"\nScrape complete: {n_rows:,} rows across {page - start_page} pages.")


# ── Run ────────────────────────────────────────────────────────────────────────

scrape_catalogue()


# ── Coverage report against KBA pairs ─────────────────────────────────────────

kba_path = Path(f"{INTERMEDIATE}/kba_neuz_model_panel.csv")
if not kba_path.exists():
    print("\nkba_neuz_model_panel.csv not found — skipping coverage report.")
else:
    adac = pd.read_csv(OUT_CSV)
    kba  = pd.read_csv(kba_path, dtype={"AGS8": str, "fab_code": str, "mod_code": str})

    kba_pairs = (
        kba[
            (kba["fab_code"] != SONSTIGE_FAB_CODE) &
            (kba["mod_code"] != SONSTIGE_MOD_CODE)
        ][["fab_text", "mod_text"]]
        .dropna()
        .drop_duplicates()
    )
    print(f"\nKBA (brand, model) pairs (non-Sonstige): {len(kba_pairs):,}")
    print(f"ADAC rows scraped:                        {len(adac):,}")
    print(f"ADAC unique model_clean values:           {adac['model_clean'].nunique():,}")

    unmapped = (
        adac[adac["energiequelle"].isna() & adac["fuel_type_adac"].notna()]
        ["fuel_type_adac"].value_counts()
    )
    if not unmapped.empty:
        print(f"\nUnmapped Kraftstoff labels — add to FUEL_MAP:")
        print(unmapped.to_string())
