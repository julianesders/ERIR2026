"""
Scrape the ADAC car catalogue for list prices and production year windows,
iterating over motorart (fuel/drive type) values and paginating within each
until the catalogue is exhausted.

Each row is tagged directly with the motorart and the corresponding KBA
energiequelle code — no post-hoc Kraftstoff-label mapping needed.

Checkpoint format: "{motorart}:{last_completed_page}" — updated after every
page so a run can be resumed after interruption without losing progress.

Inputs (optional — used only for coverage reporting at the end):
  kba_neuz_model_panel.csv   (produced by 08_clean_kba_ags8.py)

Output:
  adac_prices_raw.csv
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

ADAC_URL      = (
    "https://www.adac.de/rund-ums-fahrzeug/autokatalog/marken-modelle/autosuche/"
)
ADAC_PARAMS   = {"sort": "ALPHABETIC_ASC", "type": "klassik-autosuche"}
REQUEST_DELAY = 1.0   # seconds between requests — be polite

SONSTIGE_FAB_CODE = "000"
SONSTIGE_MOD_CODE = "0000"

# motorart URL values in order → KBA energiequelle code
MOTORART_MAP: dict[str, str] = {
    "Otto":                  "01",
    "Diesel":                "02",
    "Gas":                   "03",
    "Elektro":               "04",
    "Hybrid":                "05",
    "PlugIn-Hybrid":         "06",
    "Wankel":                "07",
    "Wasserstoff (E-Motor)": "07",
}


# ── Helpers ────────────────────────────────────────────────────────────────────

_SESSION = requests.Session()
_SESSION.headers.update({"User-Agent": "Mozilla/5.0 (compatible; research-scraper)"})


def _fetch(motorart: str, page: int) -> lxml_html.HtmlElement | None:
    params = {**ADAC_PARAMS, "motorart": motorart, "pageNumber": page}
    try:
        resp = _SESSION.get(ADAC_URL, params=params, timeout=30)
        if resp.status_code == 200:
            return lxml_html.fromstring(resp.content)
        print(f"    HTTP {resp.status_code}")
    except requests.RequestException as exc:
        print(f"    request error: {exc}")
    return None


def _parse_year(token: str) -> int | None:
    m = re.fullmatch(r"\d{2}/(\d{2})", token.strip())
    return (2000 + int(m.group(1))) if m else None


def _parse_price(s: str) -> float | None:
    s = s.replace("€", "").replace("\xa0", "").strip()
    if not s or s in ("-", "k.A."):
        return None
    try:
        return float(s.replace(".", "").replace(",", "."))
    except ValueError:
        return None


def _parse_rows(
    tree: lxml_html.HtmlElement,
    motorart: str,
    energiequelle: str,
) -> list[dict]:
    rows_out = []
    for row in tree.xpath('//tr[@data-testid="carpages:generation:model:row"]'):
        cells = {td.get("data-th"): td for td in row.xpath(".//td[@data-th]")}

        fahrzeug = cells.get("Fahrzeug")
        if fahrzeug is None:
            continue

        model_raw = " ".join(fahrzeug.text_content().split())
        preis_raw = cells["Listenpreis"].text_content().strip() if "Listenpreis" in cells else ""

        year_from = year_to = None
        yr = re.search(
            r"\((\d{2}/\d{2})\s*-\s*(\d{2}/\d{2}|[Hh]eute|[Pp]resent|k\.A\.)\)",
            model_raw,
        )
        if yr:
            year_from = _parse_year(yr.group(1))
            year_to   = _parse_year(yr.group(2))
            if year_to is None:  # "heute" / "present" / "k.A." → still on sale
                year_to = 9999

        model_clean = re.sub(r"\s*\(\d{2}/\d{2}.*?\)\s*$", "", model_raw).strip()

        rows_out.append({
            "model_raw":      model_raw,
            "model_clean":    model_clean,
            "motorart":       motorart,
            "energiequelle":  energiequelle,
            "year_from":      year_from,
            "year_to":        year_to,
            "list_price_eur": _parse_price(preis_raw),
        })
    return rows_out


# ── Resumable scrape ───────────────────────────────────────────────────────────

def scrape_catalogue() -> None:
    resume_motorart = None
    resume_page     = 0

    if Path(CHECKPOINT).exists():
        ckpt = Path(CHECKPOINT).read_text().strip()
        resume_motorart, _p = ckpt.rsplit(":", 1)
        resume_page = int(_p)
        print(f"Resuming from motorart={resume_motorart!r}, next page={resume_page + 1}")
        print(f"(delete {CHECKPOINT} to restart from scratch)\n")

    first_write = not Path(OUT_CSV).exists()
    skip        = resume_motorart is not None
    t0          = time.time()
    total_rows  = 0

    for motorart, energiequelle in MOTORART_MAP.items():

        if skip:
            if motorart == resume_motorart:
                skip       = False
                start_page = resume_page + 1
            else:
                print(f"Skipping motorart={motorart!r} (already completed)")
                continue
        else:
            start_page = 1

        print(f"\n── motorart={motorart!r}  energiequelle={energiequelle}  "
              f"starting at page {start_page} ──")

        page   = start_page
        n_rows = 0

        while True:
            tree = _fetch(motorart, page)

            if tree is None:
                print(f"  page {page}: fetch failed — stopping this motorart.")
                break

            rows = _parse_rows(tree, motorart, energiequelle)
            if not rows:
                print(f"  page {page}: empty — motorart exhausted after page {page - 1}.")
                break

            df = pd.DataFrame(rows)
            df.to_csv(OUT_CSV, mode="a", header=first_write, index=False, encoding="utf-8-sig")
            first_write  = False
            n_rows      += len(rows)
            total_rows  += len(rows)

            Path(CHECKPOINT).write_text(f"{motorart}:{page}")

            if page % 50 == 0:
                elapsed = time.time() - t0
                print(f"  page {page:>4}  |  {n_rows:>6,} rows this motorart  |  "
                      f"{total_rows:>8,} total  |  {elapsed/60:.0f} min elapsed",
                      flush=True)

            page += 1
            time.sleep(REQUEST_DELAY)

        print(f"  motorart={motorart!r} done: {n_rows:,} rows")

    if Path(CHECKPOINT).exists():
        Path(CHECKPOINT).unlink()

    elapsed = time.time() - t0
    print(f"\nScrape complete: {total_rows:,} total rows in {elapsed/60:.0f} min.")


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
    print(f"\nRows by motorart:")
    print(adac.groupby("motorart")[["model_clean", "list_price_eur"]]
          .agg(rows=("model_clean", "count"), with_price=("list_price_eur", "count"))
          .to_string())
