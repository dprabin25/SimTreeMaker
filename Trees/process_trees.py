"""
process_trees.py — batch-converts every .trees file sitting in this folder
into Newick + a per-cell SNV CSV + tree PNGs, using the same
slim_utilities.py functions TumorGrowth/simtreemaker.py and
CaseStudy/casestudy.py use after their own SLiM runs. No SLiM run happens
here at all — this is purely post-processing for .trees files you already
have (from any source: TG.slim, TG_Migration, the CaseStudy script,
elsewhere), so there's no slim_config.txt/config.txt for this one.

Usage — from inside this Trees folder:

    python process_trees.py

Scans this folder (not subfolders) for every *.trees file and, for each
one not already processed (i.e. no matching newick/<stem>.nwk yet), writes:

    newick/<stem>.nwk
    snv/<stem>_snv.csv
    pngTree/<stem>_horizontal_labels.png
    pngTree/<stem>_horizontal_no_labels.png

Already-processed files are skipped automatically — safe to drop new
.trees files in here and re-run any time without redoing old ones. Pass
--force to reprocess everything anyway:

    python process_trees.py --force

Or name specific files to process (skips the auto-scan, ignores whether
they're already processed):

    python process_trees.py someFile.trees anotherFile.trees

sim_end_tick is intentionally left for slim_utilities to auto-detect from
each .trees file's own SLiM metadata (ts.metadata["SLiM"]["tick"]), since
unlike TG.slim/casestudy.py this script has no config value telling it
what tick any given file was written at — files here can come from
different scripts/runs with different endpoints. mutation_position is
likewise mostly moot: mutation_type_label() in slim_utilities.py prefers
each mutation's own SLiM metadata over any position-based guess, and real
SLiM .trees files always carry that metadata.

Every run appends to log.txt in this folder (one line per file processed,
or the failure reason).
"""

import os
import sys
import glob
import datetime

import slim_utilities as su

MUTATION_POSITION_FALLBACK = 10000  # only used if a .trees file somehow lacks SLiM mutation metadata


def already_processed(script_dir, stem):
    return os.path.exists(os.path.join(script_dir, "newick", stem + ".nwk"))


def process_one(script_dir, trees_path, tree_opts, log_path):
    stem = os.path.splitext(os.path.basename(trees_path))[0]

    dir_newick = os.path.join(script_dir, "newick")
    dir_png    = os.path.join(script_dir, "pngTree")
    dir_snv    = os.path.join(script_dir, "snv")
    os.makedirs(dir_newick, exist_ok=True)
    os.makedirs(dir_png, exist_ok=True)
    os.makedirs(dir_snv, exist_ok=True)
    nwk_path = os.path.join(dir_newick, stem + ".nwk")
    snv_path = os.path.join(dir_snv, stem + "_snv.csv")

    horizontal_only = ["horizontal_labels", "horizontal_no_labels"]

    def _log(status, reason=None, outputs=None):
        with open(log_path, "a", encoding="utf-8") as lf:
            lf.write(f"\n=== process_trees.py: {stem}  {datetime.datetime.now()} ===\n")
            lf.write(f"Source .trees: {trees_path}\n")
            lf.write(f"Status       : {status}\n")
            if reason:
                lf.write(f"Reason       : {reason}\n")
            if outputs:
                lf.write("Outputs:\n")
                for line in outputs:
                    lf.write(line + "\n")

    print(f"\n--- Processing: {os.path.basename(trees_path)} ---")

    if not su.convert_to_newick(trees_path, nwk_path, mutation_position=MUTATION_POSITION_FALLBACK,
                                 tree_opts=tree_opts, snv_csv_path=snv_path, sim_end_tick=None):
        print(f"[ERROR] Newick conversion failed for {trees_path}")
        _log("FAILED", reason="Newick conversion raised an error:\n" +
             (su._last_newick_error or "(no traceback captured)"))
        return False

    r_ok = su.plot_png_r(nwk_path, dir_png, stem, "process_trees", "unknown", script_dir,
                          mode="Snapshot", variants=horizontal_only)
    if not r_ok:
        print("[INFO] Rscript unavailable — using matplotlib.")
        su.plot_png(nwk_path, dir_png, stem, "process_trees", "unknown",
                    mode="Snapshot", variants=horizontal_only)

    outputs = [
        f"  {nwk_path}",
        f"  {snv_path}",
        f"  {dir_png}/{stem}_horizontal_labels.png",
        f"  {dir_png}/{stem}_horizontal_no_labels.png",
    ]
    _log("OK", outputs=outputs)
    print(f"[OK] {stem}: newick, SNV CSV, and PNGs written.")
    return True


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    log_path   = os.path.join(script_dir, "log.txt")

    args  = sys.argv[1:]
    force = "--force" in args
    named = [a for a in args if a != "--force"]

    if named:
        trees_files = [
            a if os.path.isabs(a) else os.path.join(script_dir, a)
            for a in named
        ]
        missing = [p for p in trees_files if not os.path.exists(p)]
        if missing:
            for p in missing:
                print(f"[ERROR] Not found: {p}")
            sys.exit(1)
    else:
        trees_files = sorted(glob.glob(os.path.join(script_dir, "*.trees")))

    if not trees_files:
        print(f"[INFO] No .trees files found in {script_dir}.")
        return

    tree_opts = su.load_tree_options(script_dir)

    processed = skipped = failed = 0
    for trees_path in trees_files:
        stem = os.path.splitext(os.path.basename(trees_path))[0]
        if not force and not named and already_processed(script_dir, stem):
            print(f"[SKIP] {stem}: already processed (newick/{stem}.nwk exists). Use --force to redo.")
            skipped += 1
            continue
        ok = process_one(script_dir, trees_path, tree_opts, log_path)
        processed += 1 if ok else 0
        failed += 0 if ok else 1

    print(f"\n{'='*60}")
    print(f"  Done. {processed} processed, {skipped} skipped (already done), {failed} failed.")
    print(f"  Log: {log_path}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
