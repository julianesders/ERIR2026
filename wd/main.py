"""
ERIR2026 pipeline runner.

Usage:
    python main.py            # run full pipeline (INKAR, KBA, ADAC scraping skipped)
    python main.py --inkar    # re-extract INKAR from source (slow, large file)
    python main.py --kba      # include KBA steps (requires delivery files on K: drive)
    python main.py --adac     # include ADAC price scraping (~1–2 h, requires network)
    python main.py --scrape   # include the EMK scraping step (slow, network)
    python main.py --step 4   # resume from step 4 onwards
"""
import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

CODE          = Path(__file__).parent / "02_code"
DATA          = Path(__file__).parent / "01_data" / "02_intermediate"
INKAR_EXTRACT = DATA / "inkar" / "inkar_joint_panel.csv"
INKAR_PANEL   = DATA / "inkar" / "inkar_ags8_panel.csv"

STEPS = [
    # (step_number, label, script_path, skip_note_or_None)
    (0,  "Scrape EMK",                  CODE / "00_scraping/scrape_emkonzepte.py",                       "scrape-only"),
    (1,  "Clean Anschriftenverzeichnis", CODE / "01_cleaning/01_clean_anschriftenverzeichnis.py",         None),
    (2,  "Clean EMK",                   CODE / "01_cleaning/02_clean_emk.py",                            None),
    (3,  "Extract INKAR",              CODE / "01_cleaning/00_inkar_extract.py",                         "inkar-only"),
    (4,  "Build INKAR panel (AGS8)",   CODE / "01_cleaning/03_clean_inkar_ags8.py",                      "inkar-only"),
    (5,  "Clean Ladestationen (AGS8)", CODE / "01_cleaning/04_clean_ladestationen_ags8.py",              None),
    (6,  "Clean Elections (AGS8)",     CODE / "01_cleaning/05_clean_elections_ags8.py",                  None),
    (7,  "Validate AGS8 base panel",   CODE / "01_cleaning/06_build_ags8_base.py",                       None),
    (8,  "Clean Personnel (AGS5)",     CODE / "01_cleaning/07_clean_personal_ags5.py",                   None),
    (9,  "Unpack & clean KBA (AGS8)", CODE / "01_cleaning/08_clean_kba_ags8.py",                        "kba-only"),
    (10, "Aggregate KBA variables",   CODE / "01_cleaning/09_aggregate_kba_vars.py",                    "kba-only"),
    (11, "Scrape ADAC prices",        CODE / "01_cleaning/10_scrape_adac_prices.py",                    "adac-only"),
    (12, "Match EMK → AGS",           CODE / "02_merging/01_match_emk_ags.py",                          None),
    (13, "Merge EMK panel (AGS8)",    CODE / "02_merging/02_merge_emk_panel_ags8.py",                   None),
    (14, "Spatial weights (AGS8)",    CODE / "03_analysis/01_spatial_weights_ags8.py",                  None),
]


def run(label: str, script: Path) -> None:
    print(f"\n{'─'*60}")
    print(f"  Step: {label}")
    print(f"  Script: {script.relative_to(Path(__file__).parent)}")
    print(f"{'─'*60}")
    t0 = time.time()
    env = os.environ.copy()
    env["PYTHONPATH"] = str(Path(__file__).parent) + os.pathsep + env.get("PYTHONPATH", "")
    subprocess.run([sys.executable, str(script)], check=True, env=env)
    print(f"\n  Done in {time.time() - t0:.1f}s")


def main() -> None:
    parser = argparse.ArgumentParser(description="ERIR2026 pipeline runner")
    parser.add_argument(
        "--scrape", action="store_true",
        help="Include the EMK scraping step (requires network access)",
    )
    parser.add_argument(
        "--inkar", action="store_true",
        help="Re-extract INKAR variables from source (slow — reads 63M-row file)",
    )
    parser.add_argument(
        "--kba", action="store_true",
        help="Include KBA steps (requires delivery files on K: drive)",
    )
    parser.add_argument(
        "--adac", action="store_true",
        help="Include ADAC price scraping step (~1–2 h, requires network)",
    )
    parser.add_argument(
        "--step", type=int, default=0,
        help="Resume from this step number (0 = start from beginning)",
    )
    args = parser.parse_args()

    if not args.inkar:
        missing = [p for p in (INKAR_EXTRACT, INKAR_PANEL) if not p.exists()]
        if missing:
            for p in missing:
                print(f"  ERROR: missing INKAR file: {p.relative_to(Path(__file__).parent)}")
            print("  Run with --inkar to extract from source.")
            sys.exit(1)

    for step_no, label, script, skip_flag in STEPS:
        if step_no < args.step and skip_flag not in ("scrape-only", "inkar-only", "kba-only"):
            print(f"  [skip] Step {step_no}: {label}")
            continue

        if skip_flag == "scrape-only":
            if args.scrape and step_no >= args.step:
                run(label, script)
            else:
                print(f"  [skip] Step {step_no}: {label}  (pass --scrape to enable)")
            continue

        if skip_flag == "inkar-only":
            if args.inkar and step_no >= args.step:
                run(label, script)
            else:
                print(f"  [skip] Step {step_no}: {label}  (pass --inkar to re-extract)")
            continue

        if skip_flag == "kba-only":
            if args.kba and step_no >= args.step:
                run(label, script)
            else:
                print(f"  [skip] Step {step_no}: {label}  (pass --kba to enable)")
            continue

        if skip_flag == "adac-only":
            if args.adac and step_no >= args.step:
                run(label, script)
            else:
                print(f"  [skip] Step {step_no}: {label}  (pass --adac to enable)")
            continue

        run(label, script)

    print(f"\n{'='*60}")
    print("  Pipeline complete.")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
