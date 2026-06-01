#!/usr/bin/env python3
"""
Scraper for the Elektromobilitätskonzepte database at
https://elektromobilitaet-now.de/emkonzepte/#datenbank-emkonzepte

Extracts all ~328 entries via the WordPress admin-ajax.php endpoint,
parses the HTML response cards, saves structured data to CSV,
and optionally downloads all linked PDF reports.

Requirements:
    pip install requests beautifulsoup4

Usage:
    python scrape_emkonzepte.py                  # scrape data + download PDFs
    python scrape_emkonzepte.py --no-pdfs         # scrape data only, skip PDFs
    python scrape_emkonzepte.py --output-dir ./my_data  # custom output directory
"""

import argparse
import csv
import json
import os
import random
import re
import sys
import time
from pathlib import Path
from urllib.parse import urljoin

import requests
from bs4 import BeautifulSoup

BASE_URL = "https://elektromobilitaet-now.de"
AJAX_URL = f"{BASE_URL}/wp-admin/admin-ajax.php"
DATABASE_PAGE_URL = f"{BASE_URL}/emkonzepte/"

# Rotate through realistic browser User-Agent strings to avoid fingerprinting
_USER_AGENTS = [
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36",
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/125.0.0.0 Safari/537.36",
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36",
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 14_4_1) AppleWebKit/605.1.15 "
    "(KHTML, like Gecko) Version/17.4.1 Safari/605.1.15",
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:125.0) Gecko/20100101 Firefox/125.0",
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36",
]

# Base browser headers that accompany every AJAX request
_BASE_HEADERS = {
    "Accept": "application/json, text/javascript, */*; q=0.01",
    "Accept-Language": "de-DE,de;q=0.9,en-US;q=0.8,en;q=0.7",
    "Accept-Encoding": "gzip, deflate, br",
    "Referer": DATABASE_PAGE_URL,
    "X-Requested-With": "XMLHttpRequest",
    "Connection": "keep-alive",
    "DNT": "1",
}

# Delay between individual page requests (jittered around this value)
REQUEST_DELAY = 1.0  # seconds


def _pick_headers() -> dict:
    """Return headers with a freshly chosen User-Agent."""
    return {**_BASE_HEADERS, "User-Agent": random.choice(_USER_AGENTS)}


def _jitter(base: float, spread: float = 0.4) -> float:
    """Return base ± spread seconds."""
    return base + random.uniform(-spread, spread)


def warmup_session(session: requests.Session) -> None:
    """
    Visit the database page before issuing AJAX requests.

    This lets the server set any cookies or nonces a real browser would
    receive, and avoids the pattern of hitting the AJAX endpoint cold.
    """
    try:
        headers = {**_pick_headers(), "Accept": "text/html,application/xhtml+xml,*/*;q=0.8"}
        resp = session.get(DATABASE_PAGE_URL, headers=headers, timeout=30)
        resp.raise_for_status()
        print(f"  Session warmup: {resp.status_code} ({len(resp.content)} bytes)")
        time.sleep(_jitter(1.5, 0.5))  # brief pause after page load, like a real browser
    except Exception as e:
        print(f"  Session warmup failed (continuing anyway): {e}")


# Sort/order combinations to rotate across passes.
# Different orderings shift page boundaries, so gaps in one pass are
# covered by another rather than waiting for random server-side drift.
_SORT_COMBOS = [
    ("title", "asc"),
    ("title", "desc"),
    ("date", "asc"),
    ("date", "desc"),
    ("modified", "asc"),
    ("modified", "desc"),
]


def fetch_page(page: int, session: requests.Session,
               sort: str = "title", order: str = "asc") -> dict:
    """Fetch a single page of results from the AJAX endpoint."""
    params = {
        "action": "get_data",
        "postType": "emkonzepte",
        "category": "",
        "searchQuery": "",
        "page": page,
        "sort": sort,
        "order": order,
        "archive": 0,
    }
    resp = session.get(AJAX_URL, params=params, headers=_pick_headers(), timeout=30)
    resp.raise_for_status()
    return resp.json()


def parse_cards(html: str, page: int = 0) -> list[dict]:
    """
    Parse the HTML card fragments returned by the AJAX endpoint.

    Since we can't inspect the exact card HTML ahead of time, this parser
    is written to be flexible. It tries multiple strategies:
    1. Look for structured card elements (article, div with class containing 'card'/'item'/'post')
    2. Extract all text, links, and taxonomy data from each card
    3. Try to identify fields by common patterns (Bundesland, Förderkennzeichen, etc.)
    """
    soup = BeautifulSoup(html, "html.parser")
    entries = []

    # Prefer the known structure of this site to avoid duplicates
    cards = soup.select("article.emkonzepte-item")

    # Fallback for potential template changes
    if not cards:
        cards = (
            soup.find_all("article") or
            soup.find_all("div", class_=re.compile(r"card|item|post|konzept|entry|result", re.I))
        )

    if not cards:
        # Fallback: if no obvious card structure, try top-level children
        cards = [child for child in soup.children if hasattr(child, 'get_text') and child.get_text(strip=True)]

    dropped = []
    for card in cards:
        entry = extract_entry_from_card(card)
        if entry and entry.get("title"):
            entry["source_page"] = page
            entries.append(entry)
        else:
            dropped.append(card)

    if dropped:
        print(f"  [debug] page {page}: {len(dropped)} card(s) dropped (no title extracted)")
        debug_path = Path(f"debug_dropped_page{page}.html")
        if not debug_path.exists():
            debug_path.write_text(
                "\n\n<!-- NEXT CARD -->\n\n".join(str(c) for c in dropped),
                encoding="utf-8"
            )
            print(f"  [debug] saved dropped cards to {debug_path}")

    return entries


def extract_entry_from_card(card) -> dict:
    """Extract structured data from a single card element."""
    entry = {
        "source_page": 0,
        "title": "",
        "antragsteller": "",
        "url": "",
        "empfaenger": "",
        "bundesland": "",
        "raumtyp": "",
        "schlagwoerter": "",
        "foerderkennzeichen": "",
        "gesamtmittel": "",
        "bundesmittel": "",
        "laufzeit": "",
        "pdf_url": "",
        "raw_text": "",
    }

    # --- Title / Antragsteller ---
    title_el = card.find("h4") or card.find(re.compile(r"^h[1-6]$"))
    if title_el:
        entry["antragsteller"] = title_el.get_text(" ", strip=True)
        entry["title"] = entry["antragsteller"]

    # --- Link / URL ---
    # Note: on this site, wrapping <a> tags always point to PDFs, not detail pages.
    # There are no individual detail page URLs — leave entry["url"] empty.
    if card.name == "a" and card.get("href") and ".pdf" not in card["href"].lower():
        entry["url"] = urljoin(BASE_URL, card["href"])

    # --- PDF link ---
    pdf_link = card.find("a", href=re.compile(r"\.pdf", re.I))
    if pdf_link:
        entry["pdf_url"] = urljoin(BASE_URL, pdf_link["href"])
    elif card.parent and card.parent.name == "a" and card.parent.get("href"):
        parent_href = card.parent["href"]
        if ".pdf" in parent_href.lower():
            entry["pdf_url"] = urljoin(BASE_URL, parent_href)
    elif entry["url"] and ".pdf" in entry["url"].lower():
        entry["pdf_url"] = entry["url"]
    else:
        # Check all links for PDF references
        for a in card.find_all("a", href=True):
            if ".pdf" in a["href"].lower():
                entry["pdf_url"] = urljoin(BASE_URL, a["href"])
                break

    # --- Schlagwörter / Tags ---
    tags = [t.get_text(" ", strip=True) for t in card.select(".item-tags .tag") if t.get_text(strip=True)]
    if tags:
        entry["schlagwoerter"] = "; ".join(tags)

    # --- Structured key/value metadata ---
    label_map = {
        "bundesland": "bundesland",
        "empfänger": "empfaenger",
        "empfaenger": "empfaenger",
        "raumtyp": "raumtyp",
        "fkz": "foerderkennzeichen",
        "gesamtmittel": "gesamtmittel",
        "davon bundesmittel": "bundesmittel",
        "laufzeit": "laufzeit",
    }

    for p in card.select(".emkonzept-content p"):
        spans = p.find_all("span")
        if len(spans) < 2:
            continue

        label = spans[0].get_text(" ", strip=True).rstrip(":").strip().lower()
        value = spans[1].get_text(" ", strip=True)
        field = label_map.get(label)
        if field and value:
            entry[field] = value

    # --- Förderkennzeichen (funding reference) ---
    full_text = card.get_text(" ", strip=True)
    entry["raw_text"] = full_text
    fkz_match = re.search(r"03EMK\d+", full_text)
    if fkz_match:
        entry["foerderkennzeichen"] = fkz_match.group()

    # --- Taxonomy fallback: class/data matching (only if still missing) ---
    taxonomy_mappings = {
        "empfaenger": ["empfaenger", "empfänger", "recipient"],
        "bundesland": ["bundesland", "state", "land"],
        "raumtyp": ["raumtyp", "raum", "settlement"],
        "schlagwoerter": ["schlagwort", "schlagwörter", "tag", "keyword"],
    }

    for field, patterns in taxonomy_mappings.items():
        if entry.get(field):
            continue
        for pattern in patterns:
            els = card.find_all(attrs={"data-taxonomy": re.compile(pattern, re.I)})
            if els:
                entry[field] = "; ".join(el.get_text(strip=True) for el in els)
                break
            els = card.find_all(class_=re.compile(pattern, re.I))
            if els:
                values = [el.get_text(" ", strip=True) for el in els if el.get_text(strip=True)]
                if values:
                    entry[field] = "; ".join(values)
                    break

    # --- Fallback: try to detect taxonomy values from known lists in the text ---
    known_bundeslaender = [
        "Baden-Württemberg", "Bayern", "Berlin", "Brandenburg", "Bremen",
        "Hamburg", "Hessen", "Mecklenburg-Vorpommern", "Niedersachsen",
        "Nordrhein-Westfalen", "Rheinland-Pfalz", "Saarland", "Sachsen",
        "Sachsen-Anhalt", "Schleswig-Holstein", "Thüringen"
    ]
    known_empfaenger = ["Gebietskörperschaft", "Unternehmen"]

    if not entry["bundesland"]:
        for bl in known_bundeslaender:
            if bl in full_text:
                if entry["bundesland"]:
                    entry["bundesland"] += f"; {bl}"
                else:
                    entry["bundesland"] = bl

    if not entry["empfaenger"]:
        for emp in known_empfaenger:
            if emp in full_text:
                entry["empfaenger"] = emp
                break

    return entry


def _entry_key(entry: dict) -> str:
    """Stable deduplication key: prefer FKZ, fall back to title."""
    return entry.get("foerderkennzeichen") or entry.get("title", "")


def _fetch_one_pass(session: requests.Session, pass_num: int,
                    sort: str = "title", order: str = "asc") -> tuple[list[dict], int]:
    """Fetch all pages in a single pass. Returns (entries, total_data)."""
    warmup_session(session)
    first_response = fetch_page(1, session, sort=sort, order=order)

    total_pages = first_response.get("totalPages", 1)
    total_data = first_response.get("totalData", 0)

    if not first_response.get("content"):
        print(f"  Pass {pass_num}, page 1: No content returned!")
        print(f"  Raw response keys: {list(first_response.keys())}")
        print(f"  Response preview: {json.dumps(first_response, ensure_ascii=False)[:500]}")
        with open("debug_response_page1.json", "w", encoding="utf-8") as f:
            json.dump(first_response, f, ensure_ascii=False, indent=2)
        return [], total_data

    all_entries = parse_cards(first_response["content"], page=1)
    print(f"  Pass {pass_num}, page 1/{total_pages}: {len(all_entries)} entries")

    for page in range(2, total_pages + 1):
        time.sleep(_jitter(REQUEST_DELAY))
        try:
            response = fetch_page(page, session, sort=sort, order=order)
            if response.get("content"):
                entries = parse_cards(response["content"], page=page)
                all_entries.extend(entries)
                print(f"  Pass {pass_num}, page {page}/{total_pages}: {len(entries)} entries")
            else:
                print(f"  Pass {pass_num}, page {page}/{total_pages}: No content")
        except Exception as e:
            print(f"  Pass {pass_num}, page {page}/{total_pages}: ERROR - {e}")

    return all_entries, total_data


def scrape_all_entries(session: requests.Session, passes: int = 20,
                       duration_minutes: float = 60) -> list[dict]:
    """
    Fetch all pages across multiple passes and union results by FKZ/title.

    The server uses OFFSET-based pagination without a stable sort order, so
    entries near page boundaries shift between requests — some appear on two
    pages (duplicates) while others fall into gaps (missing). Running multiple
    passes spread over `duration_minutes` and taking the union recovers the
    full dataset.
    """
    seen: dict[str, dict] = {}  # key -> entry
    total_data = 0
    interval = duration_minutes * 60 / passes  # seconds between pass starts

    for pass_num in range(1, passes + 1):
        sort, order = _SORT_COMBOS[(pass_num - 1) % len(_SORT_COMBOS)]
        pass_start = time.time()
        print(f"\n--- Pass {pass_num}/{passes} (sort={sort} {order}) ---")
        entries, total_data = _fetch_one_pass(session, pass_num, sort=sort, order=order)

        before = len(seen)
        for entry in entries:
            key = _entry_key(entry)
            if key and key not in seen:
                seen[key] = entry
        new_this_pass = len(seen) - before
        print(f"  Pass {pass_num} result: {len(entries)} raw, +{new_this_pass} new unique "
              f"(total unique so far: {len(seen)})")

        if pass_num < passes:
            wait = max(0.0, interval - (time.time() - pass_start))
            if wait > 0:
                next_pass_time = time.strftime("%H:%M:%S", time.localtime(time.time() + wait))
                print(f"  Next pass at {next_pass_time} (waiting {wait:.0f}s)...")
                time.sleep(wait)

    all_entries = list(seen.values())
    print(f"\nExpected entries (from server): {total_data}")
    print(f"Unique entries collected:        {len(all_entries)}")
    if len(all_entries) < total_data:
        print(f"[!] Still missing {total_data - len(all_entries)} entries after {passes} passes. "
              f"Try --passes {passes + 2}.")

    # Save first page raw HTML for debugging/inspection
    debug_path = "debug_first_page_content.html"
    if not Path(debug_path).exists():
        pass  # already written by last _fetch_one_pass if needed

    return all_entries


def download_pdfs(entries: list[dict], output_dir: Path, session: requests.Session):
    """Download all PDF reports linked in the entries."""
    pdf_entries = [e for e in entries if e.get("pdf_url")]
    if not pdf_entries:
        print("No PDF links found in entries.")
        return

    pdf_dir = output_dir / "pdfs"
    pdf_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nDownloading {len(pdf_entries)} PDFs to '{pdf_dir}/'...")

    for i, entry in enumerate(pdf_entries, 1):
        url = entry["pdf_url"]
        filename = url.split("/")[-1]
        if not filename.endswith(".pdf"):
            filename += ".pdf"

        filepath = pdf_dir / filename

        if filepath.exists():
            print(f"  [{i}/{len(pdf_entries)}] Already exists: {filename}")
            continue

        try:
            time.sleep(REQUEST_DELAY)
            resp = session.get(url, headers=_pick_headers(), timeout=60)
            resp.raise_for_status()
            filepath.write_bytes(resp.content)
            size_mb = len(resp.content) / (1024 * 1024)
            print(f"  [{i}/{len(pdf_entries)}] Downloaded: {filename} ({size_mb:.1f} MB)")
        except Exception as e:
            print(f"  [{i}/{len(pdf_entries)}] FAILED: {filename} - {e}")


def save_csv(entries: list[dict], output_path: Path):
    """Save entries to CSV."""
    if not entries:
        print("No entries to save.")
        return

    fieldnames = [
        "source_page", "title", "antragsteller", "url", "empfaenger", "bundesland", "raumtyp",
        "laufzeit", "gesamtmittel", "bundesmittel", "schlagwoerter",
        "foerderkennzeichen", "pdf_url", "raw_text"
    ]

    with open(output_path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(entries)

    print(f"\nSaved {len(entries)} entries to '{output_path}'")


def try_detail_pages(entries: list[dict], session: requests.Session) -> list[dict]:
    """
    If entries have individual URLs, fetch them to extract richer data.
    Only runs if the initial card parsing left fields mostly empty.
    """
    fields_to_check = ["empfaenger", "bundesland", "raumtyp", "schlagwoerter"]
    filled = sum(1 for e in entries if any(e.get(f) for f in fields_to_check))

    if filled > len(entries) * 0.3:
        # More than 30% of entries have taxonomy data from cards — skip detail pages
        return entries

    urls_available = [e for e in entries if e.get("url") and "/emkonzepte" in e["url"]]
    if not urls_available:
        return entries

    print(f"\nCard data is sparse — fetching {len(urls_available)} detail pages for richer data...")

    for i, entry in enumerate(urls_available, 1):
        try:
            time.sleep(REQUEST_DELAY)
            resp = session.get(entry["url"], headers=_pick_headers(), timeout=30)
            resp.raise_for_status()
            soup = BeautifulSoup(resp.text, "html.parser")

            full_text = soup.get_text(" ", strip=True)

            # Try to find PDF link on detail page
            if not entry.get("pdf_url"):
                pdf_link = soup.find("a", href=re.compile(r"\.pdf", re.I))
                if pdf_link:
                    entry["pdf_url"] = urljoin(BASE_URL, pdf_link["href"])

            # Try to find Förderkennzeichen
            if not entry.get("foerderkennzeichen"):
                fkz = re.search(r"03EMK\d{4}", full_text)
                if fkz:
                    entry["foerderkennzeichen"] = fkz.group()

            # Try to find taxonomy terms from meta or tag elements
            for tax_class, field in [
                ("empfaenger", "empfaenger"), ("bundesland", "bundesland"),
                ("raumtyp", "raumtyp"), ("schlagwort", "schlagwoerter")
            ]:
                if not entry.get(field):
                    els = soup.find_all(class_=re.compile(tax_class, re.I))
                    if els:
                        entry[field] = "; ".join(el.get_text(strip=True) for el in els)

            if i % 10 == 0 or i == len(urls_available):
                print(f"  Fetched {i}/{len(urls_available)} detail pages")

        except Exception as e:
            print(f"  Detail page {i} FAILED: {e}")

    return entries


def main():
    parser = argparse.ArgumentParser(description="Scrape Elektromobilitätskonzepte database")
    parser.add_argument("--no-pdfs", action="store_true", help="Skip downloading PDF reports")
    parser.add_argument("--output-dir", type=str, default=".",
                        help="Output directory (default: current working directory)")
    parser.add_argument("--delay", type=float, default=1.0,
                        help="Delay between requests in seconds (default: 1.0)")
    parser.add_argument("--passes", type=int, default=20,
                        help="Number of full scrape passes to union results (default: 20)")
    parser.add_argument("--duration", type=float, default=60,
                        help="Total duration in minutes to spread passes over (default: 60)")
    args = parser.parse_args()

    global REQUEST_DELAY
    REQUEST_DELAY = args.delay

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    session = requests.Session()

    print("=" * 60)
    print("  Elektromobilitätskonzepte Database Scraper")
    print("=" * 60)
    print(f"  Target:  {AJAX_URL}")
    print(f"  Output:  {output_dir.resolve()}")
    print(f"  Delay:   {REQUEST_DELAY}s between requests")
    print(f"  Passes:  {args.passes} over {args.duration:.0f} min "
          f"(every {args.duration * 60 / args.passes:.0f}s)")
    print(f"  PDFs:    {'skip' if args.no_pdfs else 'download'}")
    print("=" * 60 + "\n")

    # Step 1: Scrape all entries (multi-pass to handle unstable server-side ordering)
    entries = scrape_all_entries(session, passes=args.passes, duration_minutes=args.duration)

    if not entries:
        print("\n[!] No entries were parsed. Check debug files for the raw response.")
        print("    You may need to adjust the parse_cards() function based on the actual HTML.")
        sys.exit(1)

    # Step 2: Optionally enrich via detail pages
    entries = try_detail_pages(entries, session)

    # Step 3: Save all entries to CSV
    csv_path = output_dir / "emkonzepte_database.csv"
    save_csv(entries, csv_path)

    from collections import Counter
    fkz_counts = Counter(e.get("foerderkennzeichen") or "" for e in entries)
    no_fkz_entries = [e for e in entries if not e.get("foerderkennzeichen")]
    duplicate_fkz_entries = [e for e in entries if fkz_counts[e.get("foerderkennzeichen") or ""] > 1]

    # Step 4: Optionally download PDFs
    if not args.no_pdfs:
        download_pdfs(entries, output_dir, session)

    # Summary
    print("\n" + "=" * 60)
    print("  DONE")
    print("=" * 60)
    print(f"  Total entries:       {len(entries)}")
    print(f"  Duplicate FKZ:       {len(duplicate_fkz_entries)}")
    print(f"  No FKZ:              {len(no_fkz_entries)}")
    print(f"  With PDF links:      {sum(1 for e in entries if e.get('pdf_url'))}")
    print(f"  With Bundesland:     {sum(1 for e in entries if e.get('bundesland'))}")
    print(f"  CSV:                 {csv_path.resolve()}")
    if not args.no_pdfs:
        print(f"  PDFs saved to:    {(output_dir / 'pdfs').resolve()}")


if __name__ == "__main__":
    main()
