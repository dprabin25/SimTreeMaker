"""
casestudy.py — runs Simulation4_ClonalHematopoiesis.slim (a patient-calibrated
clonal hematopoiesis case study — see that file's own header for the model
and its source: Mitchell et al. 2022, Nature), then converts the resulting
.trees file to Newick + a per-cell SNV CSV + tree PNGs using the same
slim_utilities.py functions TumorGrowth/simtreemaker.py uses for TG.slim.

This is a separate, self-contained pipeline — nothing here reads or writes
anything under TumorGrowth/, so it can't affect the migration/no-migration
setup over there.

No config.txt: this case study is just one .slim file with three -d
constants, so AGE/DRIVER_TICK/REP are edited straight from the terminal —
no config layer needed for something this small.

    AGE          generation the simulation stops at (the patient's age at
                 sampling; must be <= 500, the script's own hardcoded
                 tick-range ceiling)
    DRIVER_TICK  generation a driver (m2) mutation is introduced into a
                 random newborn cell
    REP          replicate label, used only to keep output filenames from
                 different runs of the same AGE apart (e.g. running the
                 same scenario multiple times to see stochastic
                 variability) — note this does NOT seed SLiM's RNG, so two
                 runs with the same AGE/DRIVER_TICK/REP are not guaranteed
                 to reproduce identical results; SLiM draws its own random
                 seed each run unless you pass its separate global -seed
                 flag (not one of these three -d constants).

  1. Directly from the terminal, supplying every -d constant yourself —
     run from inside this CaseStudy folder, since the .slim file writes
     its own output using a relative path ("CaseStudyOutputs/..."), not a
     treesPath -d flag like TG.slim has:

       cd CaseStudy
       slim -d AGE=81 -d DRIVER_TICK=38 -d REP=1 Simulation4_ClonalHematopoiesis.slim

  2. Through this script, which reads slim_config.txt for the SLiM
     executable path (that's the only file it reads — see above), then
     also does the Newick/SNV/PNG conversion afterward:

       python casestudy.py                                   — uses the
           defaults below (AGE=81, DRIVER_TICK=38, REP=1)

       python casestudy.py AGE=85 DRIVER_TICK=40 REP=2        — override
           any of the three straight on the command line, in any order,
           any subset

     This script prints the equivalent terminal command it built before
     running it, and changes its own working directory to this folder
     first (os.chdir below) so the .slim file's relative output path
     always resolves here, regardless of where `python casestudy.py` was
     actually invoked from.

Every run appends an entry to CaseStudyOutputs/log.txt (the SLiM command
used, SLiM's own console output, and the resulting output paths, or the
failure reason).
"""

import os
import sys
import datetime

import slim_utilities as su

SCRIPT_NAME   = "Simulation4_ClonalHematopoiesis"
SLIM_FILE     = f"{SCRIPT_NAME}.slim"
OUTPUT_SUBDIR = "CaseStudyOutputs"

# Defaults, straight from the .slim file's own example usage comment.
# Override any of these from the terminal, e.g.: python casestudy.py AGE=85
DEFAULTS = {
    "AGE": "81",
    "DRIVER_TICK": "38",
    "REP": "1",
}

# Matches genome.addNewDrawnMutation(m2, 10000) inside the .slim file —
# used only as the fallback label position if a .trees file is ever read
# back without SLiM's own mutation-type metadata (real runs always carry
# that metadata, via slim_utilities.mutation_type_label(), so this has no
# effect in practice).
MUTATION_POSITION = 10000


def build_d_args(params):
    d_args = []

    def _d(name, value):
        d_args.extend(["-d", f"{name}={value}"])

    _d("AGE", params["AGE"])
    _d("DRIVER_TICK", params["DRIVER_TICK"])
    _d("REP", params["REP"])
    return d_args


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Simulation4_ClonalHematopoiesis.slim writes its own output path as a
    # literal relative string ("CaseStudyOutputs/...") — it has no
    # treesPath -d flag the way TG.slim does. Changing this process's cwd
    # to this folder is the simplest way to make that relative path always
    # resolve here, regardless of where `python casestudy.py` was invoked
    # from.
    os.chdir(script_dir)

    config   = su.load_config(os.path.join(script_dir, "slim_config.txt"))
    slim_exe = config.get("SLIM_EXE", "").strip()
    if not slim_exe:
        print("[ERROR] slim_config.txt must define SLIM_EXE.")
        sys.exit(1)
    if sys.platform == "win32" and not slim_exe.lower().endswith(".exe"):
        slim_exe += ".exe"

    # AGE/DRIVER_TICK/REP: start from DEFAULTS, then apply any KEY=VALUE
    # given on the command line, e.g. `python casestudy.py AGE=85 REP=2`.
    params = dict(DEFAULTS)
    cli_overrides = {}
    for a in sys.argv[1:]:
        if "=" in a:
            k, _, v = a.partition("=")
            cli_overrides[k.strip()] = v.strip()
    if cli_overrides:
        unknown = set(cli_overrides) - set(DEFAULTS)
        if unknown:
            print(f"[ERROR] Unknown parameter(s): {', '.join(sorted(unknown))} "
                  f"(expected: {', '.join(DEFAULTS)}).")
            sys.exit(1)
        print(f"[INFO] CLI overrides: {cli_overrides}")
        params.update(cli_overrides)

    age, driver_tick, rep = params["AGE"], params["DRIVER_TICK"], params["REP"]

    slim_path = os.path.join(script_dir, SLIM_FILE)
    if not os.path.exists(slim_path):
        print(f"[ERROR] Script not found: {slim_path}")
        sys.exit(1)

    out_dir = os.path.join(script_dir, OUTPUT_SUBDIR)
    os.makedirs(out_dir, exist_ok=True)  # SLiM will not create this itself

    # Must exactly match the Eidos string the .slim file builds internally:
    # "CaseStudyOutputs/Sim4_CHIP" + AGE + "_rep" + REP + ".trees"
    stem       = f"Sim4_CHIP{age}_rep{rep}"
    trees_path = os.path.join(out_dir, stem + ".trees")
    log_path   = os.path.join(out_dir, "log.txt")

    d_args   = build_d_args(params)
    full_cmd = f"{slim_exe} " + " ".join(d_args) + f" {SLIM_FILE}"

    def _log(status, reason=None, slim_output="", outputs=None):
        with open(log_path, "a", encoding="utf-8") as lf:
            lf.write(f"\n=== casestudy.py Run: {datetime.datetime.now()} ===\n")
            lf.write(f"Script      : {SLIM_FILE}  (AGE={age}, DRIVER_TICK={driver_tick}, REP={rep})\n")
            lf.write(f"SLiM command: {full_cmd}\n")
            lf.write(f"Status      : {status}\n")
            if reason:
                lf.write(f"Reason      : {reason}\n")
            if slim_output:
                lf.write("SLiM output:\n")
                lf.write(slim_output)
            if outputs:
                lf.write("\nOutputs:\n")
                for line in outputs:
                    lf.write(line + "\n")

    print(f"\n{'='*60}")
    print(f"  Script : {SLIM_FILE}")
    print(f"  AGE={age}  DRIVER_TICK={driver_tick}  REP={rep}")
    print(f"  Output : {trees_path}")
    print(f"{'='*60}")
    print(f"[INFO] Equivalent terminal command (run from inside CaseStudy/):")
    print(f"    {full_cmd}\n")

    ok, slim_output = su.run_slim(slim_exe, slim_path, extra_args=d_args)
    if not ok:
        print(f"[SKIP] Skipping Newick/PNG conversion due to SLiM error.")
        _log("FAILED", reason="SLiM exited with a non-zero return code — see SLiM output above/below.",
             slim_output=slim_output)
        return

    if not os.path.exists(trees_path):
        print(f"[ERROR] SLiM reported success but the expected output wasn't found: {trees_path}")
        _log("FAILED", reason=f"Expected .trees file not found: {trees_path}", slim_output=slim_output)
        return

    # -- Post-processing: .trees -> .nwk -> horizontal PNGs (labels + no labels) --
    tree_opts    = su.load_tree_options(script_dir)
    sim_end_tick = int(age)

    dir_newick = os.path.join(out_dir, "newick")
    dir_png    = os.path.join(out_dir, "pngTree")
    dir_snv    = os.path.join(out_dir, "snv")
    os.makedirs(dir_newick, exist_ok=True)
    os.makedirs(dir_png, exist_ok=True)
    os.makedirs(dir_snv, exist_ok=True)
    nwk_path = os.path.join(dir_newick, stem + ".nwk")
    snv_path = os.path.join(dir_snv, stem + "_snv.csv")

    horizontal_only = ["horizontal_labels", "horizontal_no_labels"]

    if not su.convert_to_newick(trees_path, nwk_path, mutation_position=MUTATION_POSITION,
                                 tree_opts=tree_opts, snv_csv_path=snv_path,
                                 sim_end_tick=sim_end_tick):
        _log("FAILED",
             reason="Newick conversion raised an error:\n" + (su._last_newick_error or "(no traceback captured)"),
             slim_output=slim_output)
        return

    r_ok = su.plot_png_r(nwk_path, dir_png, stem, SCRIPT_NAME, "nonWF", script_dir,
                          mode="Snapshot", variants=horizontal_only)
    if not r_ok:
        print("[INFO] Rscript unavailable — using matplotlib.")
        su.plot_png(nwk_path, dir_png, stem, SCRIPT_NAME, "nonWF",
                    mode="Snapshot", variants=horizontal_only)

    outputs = [
        f"  {trees_path}",
        f"  {nwk_path}",
        f"  {snv_path}",
        f"  {dir_png}/{stem}_horizontal_labels.png",
        f"  {dir_png}/{stem}_horizontal_no_labels.png",
    ]
    _log("OK", slim_output=slim_output, outputs=outputs)

    print(f"\n{'='*60}")
    print(f"  Run complete.")
    print(f"  Tree   : {trees_path}")
    print(f"  Newick : {nwk_path}")
    print(f"  SNV CSV: {snv_path}")
    print(f"  PNGs   : {dir_png}/{stem}_horizontal_labels.png / _no_labels.png")
    print(f"  Log    : {log_path}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
