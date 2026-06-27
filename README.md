# ERIR2026 — Replication Files

Replication package for "Communal Projects and EV Charging Inequality" (Julian Esders, Seminar ERIR 2026). The paper estimates the effect of German federal mobility-concept (EMK) funding on Battery-Electric-Vehicle (BEV) adoption at the German municipality (Gemeinde / AGS8) level using a staggered Callaway-Sant'Anna difference-in-differences design, complemented by discrete-time hazard models of funding onset.

---

## Software requirements

| Tool | Version tested | Purpose |
|------|---------------|---------|
| Python | ≥ 3.10 | Data cleaning and merging |
| R | ≥ 4.3 | Analysis and table/figure output |
| LaTeX | any modern TeX Live / MiKTeX | Compiling the paper |

**Python packages:** `pandas`, `numpy`, `requests`, `beautifulsoup4`, `openpyxl`, `geopandas`, `libpysal`, `scipy`

**R packages:** `data.table`, `did`, `fixest`, `ggplot2`, `marginaleffects`, `patchwork`, `ragg`, `scales`, `sf`

---

## Data requirements

All data is included in the repository **except** for two sources that are either large or subject to licensing restrictions:

### INKAR panel (BBSR) — download required

Population, density, fiscal capacity, purchasing power, and demographic structure at AGS8/AGS5 level, 1995–2023. Download the raw export from the [INKAR portal](https://www.inkar.de) (free registration required). The expected file is a single large CSV; place it at:

```
wd/01_data/01_raw/INKAR_2025.csv
```

Update the filename in `wd/02_code/01_cleaning/00_inkar_extract.py` if the download date differs.

### KBA Fahrzeugregister — place delivery files

BEV and ICE new registrations and vehicle stock at AGS8 level by fuel type and holder (corporate/private), 2012–2023. Place the KBA delivery zip archives in:

```
wd/01_data/01_raw/kba/
```

The cleaning script `wd/02_code/01_cleaning/08_clean_kba_ags8.py` expects fixed-width text files unpacked from those zips.

---

## Setup

1. **Clone the repository.**

2. **Configure the base path:**
   ```
   cp wd/config.example.py wd/config.py
   # Edit BASE_PATH in config.py to point to your local wd/ directory.
   ```

3. **Place the two external data sources** as described above.

---

## Run order

All steps can be triggered via the top-level orchestrator:

```bash
python wd/main.py
```

Or run each stage individually:

```bash
# Stage 1 — Python: clean each source (run in numbered order)
python wd/02_code/01_cleaning/00_inkar_extract.py
python wd/02_code/01_cleaning/01_clean_anschriftenverzeichnis.py
# ... through 10_clean_area_ags8.py

# Stage 2 — Python: match EMK to AGS and build the panel
python wd/02_code/02_merging/01_match_emk_ags.py
python wd/02_code/02_merging/02_merge_emk_panel_ags8.py

# Stage 3 — R: build estimation frames
Rscript wd/02_code/03_analysis/00_prep_analysis.R

# Stage 4 — R: run analyses (any order within this stage)
Rscript wd/02_code/03_analysis/01_descriptives.R
Rscript wd/02_code/03_analysis/02_hazard.R
Rscript wd/02_code/03_analysis/03_logit_uncensored.R
Rscript wd/02_code/03_analysis/04_did_main.R
Rscript wd/02_code/03_analysis/05_did_anticipation.R
Rscript wd/02_code/03_analysis/06_heterogeneity.R
Rscript wd/02_code/03_analysis/07_spillovers.R
```

All R outputs (tables as CSV + TeX, figures as PNG) are written to `wd/03_output/<script_stem>/`.

---

## Repository layout

```
ERIR2026/
├── wd/
│   ├── 01_data/
│   │   ├── 01_raw/          # raw source files (KBA and INKAR not tracked)
│   │   ├── 02_clean/        # intermediate cleaned files per source
│   │   └── 03_final/        # merged panel + estimation frames
│   ├── 02_code/
│   │   ├── 00_scraping/     # EMK web scraper
│   │   ├── 01_cleaning/     # per-source Python cleaners
│   │   ├── 02_merging/      # EMK–AGS matching + panel build
│   │   └── 03_analysis/     # R analysis pipeline
│   ├── 03_output/           # generated tables, figures, CSVs (not tracked)
│   ├── config.example.py    # path config template
│   └── main.py              # pipeline orchestrator
├── documents/               # paper drafts, slides, methodology notes
└── materials/               # source PDFs and reference material
```

For a detailed walkthrough of each cleaning and merging step see `documents/pipeline_summary.md`. For the full econometric specification see `documents/econometric_specification.md`.
