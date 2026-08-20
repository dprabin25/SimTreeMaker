"""
simtreemaker.py — runs TG.slim directly for no-migration. For migration
runs, this file itself builds the migration variant (in memory, from
TG.slim + config.txt — no separate generator module) and runs it from a
short-lived temporary .slim file that is deleted immediately after SLiM
finishes. Nothing named TG_Migration.slim (or similar) is ever saved
permanently in this folder.

TG.slim is a plain no-migration backbone with no migrationEnabled flag and
no p2/p3 code at all — there's no way to turn migration on by typing extra
-d flags at TG.slim itself. Migration is Python-only:

  1. Directly from the terminal, supplying every -d constant yourself —
     this always means no migration, since TG.slim has no migration branch
     to turn on:

       slim -d initialPopSize=1 -d targetSize=150 -d mutationRate=1e-8 \
            -d driverMutationRate=1e-8 -d m2_dominance=0.5 -d m2_effect=0.6 \
            -d simulationEndTick=200 -d treesPath="'TG.trees'" \
            TG.slim

  2. Through this script, which reads slim_config.txt (for the SLiM
     executable path) and config.txt (for the parameter values) instead of
     you typing the -d flags by hand:

       python simtreemaker.py

     config.txt's migrationEnabled=F is the default — plain
     `python simtreemaker.py` runs TG.slim, no migration. To run migration
     for a single invocation without editing config.txt, add `migration`
     on the command line — this builds the migration variant fresh from
     the current TG.slim + config.txt, writes it to a temporary .slim
     file, runs that instead of TG.slim, then deletes the temp file:

       python simtreemaker.py migration

     Any other config.txt key can be overridden the same way, e.g. to also
     move where p2 shows up for just this run:

       python simtreemaker.py migration Pop2Time=30

     This script prints the equivalent terminal command it built before
     running it, so you can see exactly what it's doing (including the
     path of the temporary migration script, when applicable).

After SLiM finishes, this script converts the resulting .trees file to
Newick and writes a per-cell SNV profile CSV from that same sampled tree
(both via slim_utilities.py), and renders horizontal tree PNGs — with
labels and without — using R+ape if available, falling back to matplotlib
otherwise. Vertical PNGs aren't generated here (that machinery still lives
in slim_utilities.py if needed later). Post-processing samples from
whatever populations actually exist in the .trees file, so migration runs
(p1+p2+p3+p4) and no-migration runs (p1 only) both work without any change
here. How many cells get sampled per population is controlled by
config.txt's Pop1SampleCells..Pop4SampleCells (build_sample_cells()) —
leave any/all unset to fall back to the pooled snapshotSamples cap in
TreeOptions/Snapshot.txt instead.

Every run also appends an entry to TG/log.txt — the SLiM command used,
SLiM's own console output, and the resulting output paths (or, on failure,
the reason) — so a run can be audited after the fact without needing to
keep terminal scrollback around.
"""

import os
import re
import sys
import tempfile
import datetime

import slim_utilities as su

SCRIPT_NAME = "TG"
SOURCE_SLIM = f"{SCRIPT_NAME}.slim"   # "TG.slim"


# ── Migration script building (in memory — no separate generator module,
#    no permanent output file; see build_migration_slim_text()) ───────────

MARKER = (
    "    ////////////////////////////////\n"
    "    //Add migration setting if any//\n"
    "    ////////////////////////////////\n"
)

_POP_TEMPLATE = """    if (sim.cycle >= Pop{n}Time) {{
        // Pop{n}Origin is an integer population ID (1 = p1, 2 = p2, ...) —
        // written in config.txt as a readable "PopM" reference (e.g.
        // Pop{n}Origin = Pop2) and converted to that integer by
        // simtreemaker.py before being passed here via -d. Looked up
        // dynamically every generation so p{n} can originate from any
        // population that exists by then, not just p1. If the configured
        // origin doesn't exist yet (e.g. its own PopTime hasn't happened),
        // this is empty and migration into p{n} is simply skipped that
        // generation.
        source{n} = sim.subpopulations[sim.subpopulations.id == Pop{n}Origin];
        if (sim.cycle == Pop{n}Time) {{
            // Founding pull: seed p{n} from its configured source population.
            sim.addSubpop("p{n}", 0);
            if (size(source{n}) > 0) {{
                nMig = min(Pop{n}IniPopSize, source{n}.individualCount);
                migrants = sample(source{n}.individuals, nMig);
                p{n}.takeMigrants(migrants);
            }}
        }} else if (size(source{n}) > 0) {{
            // Ongoing migration: nonWF models have no automatic
            // per-generation migration, so each subsequent generation
            // draws its own migrant count from the source population,
            // sized by Pop{n}MigrationRate (per-individual probability of
            // migrating this generation).
            nOngoing = rbinom(1, source{n}.individualCount, Pop{n}MigrationRate);
            if (nOngoing > 0)
                p{n}.takeMigrants(sample(source{n}.individuals, nOngoing));
        }}
        n{n} = p{n}.individualCount;
        if (n{n} > 0)
            p{n}.fitnessScaling = Pop{n}TargetSize / n{n};
    }}
"""


def build_migration_block(pop_ids):
    return (
        "    ////////////////////////////////\n"
        "    //Add migration setting if any//\n"
        "    ////////////////////////////////\n"
        + "".join(_POP_TEMPLATE.format(n=n) for n in pop_ids)
    )


HEADER_OLD = """// TG — no-migration, single-population TumorGrowth backbone.
//
// This file has NO migration branching at all — there is no
// migrationEnabled flag, no p2/p3 code, nothing conditional. It is the
// backbone that generate_migration_slim.py reads to build a migration
// variant (TG_Migration.slim) by inserting Eidos code at the
// "//Add migration setting if any//" marker in early() below. That
// marker is otherwise inert — three comment lines, no behavior change —
// so running TG.slim directly always means no migration, full stop.
//
// Migration is Python-only: there is no supported way to turn migration
// on by hand-typing extra -d flags at this file. Use simtreemaker.py:
//
//   python simtreemaker.py              — no migration (runs TG.slim)
//   python simtreemaker.py migration    — migration (generates/refreshes
//                                          TG_Migration.slim from this file,
//                                          then runs that instead)
//
// Growth model: every individual produces exactly one cloned offspring per
// generation (reproduction() below) and is NOT explicitly killed afterward
// — population size is regulated entirely by nonWF density-dependent
// viability selection via fitnessScaling (early() below), not by a forced
// per-age-cohort kill. There is no guaranteed founder driver mutation
// either — all m2 driver mutations arise spontaneously via
// driverMutationRate in the g2 region (positions 50001-99999); nothing is
// manually placed at generation 1.
//
// Required -d constants:
//   initialPopSize, targetSize, mutationRate, driverMutationRate,
//   m2_dominance, m2_effect, simulationEndTick, treesPath
//
//   Direct terminal run:
//     slim -d initialPopSize=1 -d targetSize=150 -d mutationRate=1e-8 \\
//          -d driverMutationRate=1e-8 -d m2_dominance=0.5 -d m2_effect=0.6 \\
//          -d simulationEndTick=200 -d treesPath="'TG.trees'" \\
//          TG.slim
//
// Running this file in SLiMgui (or via `slim` with no -d flags) will fail
// with "undefined identifier" errors — that's expected; supply all eight
// -d flags listed above."""


def build_header(pop_ids):
    pop_list = ", ".join(f"p{n}" for n in pop_ids)
    example_flags = " \\\n".join(
        f"//        -d Pop{n}Time=50 -d Pop{n}IniPopSize=5 -d Pop{n}TargetSize=100 "
        f"-d Pop{n}MigrationRate=0.01 -d Pop{n}Origin=1"
        for n in pop_ids
    )
    n_flags = 8 + 5 * len(pop_ids)
    return f"""// TG_Migration variant — built in memory by simtreemaker.py from TG.slim
// + config.txt, and written to a short-lived temporary .slim file for this
// run only. That temp file is deleted automatically as soon as SLiM
// finishes (success or failure) — nothing named TG_Migration.slim (or
// similar) is ever saved permanently in this folder. If you're reading
// this comment, you're looking at that temporary file while it still
// exists; don't hand-edit it, since it won't be here after the run and
// any edits would be lost. To change the migration model, edit TG.slim
// (the backbone) or the Pop* template inside simtreemaker.py, not this
// file.
//
// Adds {len(pop_ids)} population(s), {pop_list}: each is seeded by a one-time
// migrant pull from its configured PopNOrigin at a configurable
// generation, then receives further migrants every generation after that
// (this is a nonWF model, so migration is drawn manually every generation
// via rbinom() rather than through setMigrationRates() — that method only
// drives automatic migration in WF models, and calling it here would
// raise a runtime error). Note: below, the -d flags show PopNOrigin as a
// plain integer — that's what SLiM actually receives; config.txt itself
// uses a readable "PopM" reference instead (e.g. Pop4Origin = Pop2),
// which simtreemaker.py converts to this integer before invoking SLiM.
//
// Required -d constants: everything TG.slim needs, plus PopNTime,
// PopNIniPopSize, PopNTargetSize, PopNMigrationRate, and PopNOrigin
// (1 = p1, 2 = p2, ...) for N in {", ".join(str(n) for n in pop_ids)}:
//
//   slim -d initialPopSize=1 -d targetSize=150 -d mutationRate=1e-8 \\
//        -d driverMutationRate=1e-8 -d m2_dominance=0.5 -d m2_effect=0.6 \\
//        -d simulationEndTick=200 -d treesPath="'TG.trees'" \\
{example_flags} \\
//        <this temporary file>
//
// Running this file directly in SLiMgui (or via `slim` with no -d flags)
// will fail with "undefined identifier" errors — that's expected; supply
// all {n_flags} -d flags listed above."""


def insert_migration_block(text, pop_ids):
    """Replace the inert 3-line marker inside early() with the real
    migration Eidos code. Raises ValueError if the marker isn't found, so
    a drifted/hand-edited base file fails loudly instead of silently
    producing a script with no migration logic in it."""
    if MARKER not in text:
        raise ValueError(
            f"Marker not found in {SOURCE_SLIM} — expected exactly:\n{MARKER!r}\n"
            "(base script may have been edited since this generator was written)"
        )
    return text.replace(MARKER, build_migration_block(pop_ids), 1)


def rewrite_header(text, pop_ids):
    if HEADER_OLD not in text:
        raise ValueError(f"Header comment block not found/changed in {SOURCE_SLIM}.")
    return text.replace(HEADER_OLD, build_header(pop_ids), 1)


def build_migration_slim_text(script_dir, pop_ids):
    """
    Build the full migration-variant Eidos script text in memory, from the
    current TG.slim + the given pop_ids. Returns a string — the caller is
    responsible for writing it to a temporary file and cleaning that file
    up afterward (see main()). Explicit encoding="utf-8" matters here:
    TG.slim's header comment contains em-dash characters, and on Windows
    the platform default text encoding is not UTF-8 — opening without
    encoding="utf-8" would silently corrupt those characters and make the
    HEADER_OLD match in rewrite_header() fail.
    """
    src_path = os.path.join(script_dir, SOURCE_SLIM)
    with open(src_path, encoding="utf-8") as f:
        text = f.read()
    text = rewrite_header(text, pop_ids)
    text = insert_migration_block(text, pop_ids)
    return text


def detect_migration_pops(params):
    """
    Scan config.txt's parsed params for any PopNTime key (N >= 2) and
    return the sorted list of population numbers found — this is what
    makes the number of migration populations config-driven instead of
    hardcoded. Add a 5th population by adding a Pop5Time/Pop5IniPopSize/
    Pop5TargetSize/Pop5MigrationRate block to config.txt; no code changes
    needed anywhere in this file.
    """
    pops = set()
    for key in params:
        m = re.fullmatch(r"Pop(\d+)Time", key)
        if m:
            n = int(m.group(1))
            if n >= 2:
                pops.add(n)
    return sorted(pops)


def parse_pop_origin(value, n):
    """
    Parse a PopNOrigin config value into the plain integer population id
    the migration script actually needs. config.txt writes this as a
    readable population reference — "Pop1", "pop1", or the shorthand "p1"
    all work (case-insensitive, "op" optional) — or a bare integer (e.g.
    Pop4Origin = 2) if preferred. Exits with an error message if nothing
    matches, since a silently-wrong origin id would just make p{n} draw
    migrants from the wrong (or no) population.
    """
    value = value.strip()
    m = re.fullmatch(r"[Pp](?:op)?(\d+)", value)
    if m:
        return int(m.group(1))
    if re.fullmatch(r"\d+", value):
        return int(value)
    print(f"[ERROR] Pop{n}Origin={value!r} is not a valid population reference "
          f"(expected e.g. 'Pop1', 'p1', or a bare integer like 1).")
    sys.exit(1)


def build_d_args(params, trees_path, migration_enabled, pop_ids):
    """
    Build the -d flag list, from config.txt's parsed params dict — the
    same flags you'd type by hand if running SLiM directly from the
    terminal. No migrationEnabled flag is ever passed — TG.slim doesn't
    have that symbol at all, and the migration variant doesn't need it
    either (its migration code always runs). The PopN* flags (for each n
    in pop_ids) are only added when migration_enabled is True. PopNOrigin
    is parsed from config.txt's readable "PopM" form into the plain
    integer SLiM needs (see parse_pop_origin()) before being passed as -d.
    """
    d_args = []

    def _d(name, value):
        d_args.extend(["-d", f"{name}={value}"])

    target_val = params["targetSize"].split(",")[0].strip()

    _d("initialPopSize", params["initialPopSize"])
    _d("targetSize", target_val)
    _d("mutationRate", params["mutationRate"])
    _d("driverMutationRate", params["driverMutationRate"])
    _d("m2_dominance", params["m2_dominance"])
    _d("m2_effect", params["m2_effect"])
    _d("simulationEndTick", params["simulationEndTick"])
    _d("treesPath", f"'{trees_path}'")

    if migration_enabled:
        for n in pop_ids:
            for suffix in ("Time", "IniPopSize", "TargetSize", "MigrationRate", "Origin"):
                key = f"Pop{n}{suffix}"
                if key not in params:
                    print(f"[ERROR] migration requested but config.txt is missing {key}.")
                    sys.exit(1)
                if suffix == "Origin":
                    _d(key, parse_pop_origin(params[key], n))
                else:
                    _d(key, params[key])

    return d_args


def build_sample_cells(params):
    """
    Build the {population_id: count} dict for per-population Newick/SNV
    sampling, from any Pop<N>SampleCells key in config.txt — not limited
    to p1-p4, so a Pop5SampleCells (or higher) works the same way.
    Populations without a SampleCells key just aren't included — the ones
    that don't exist for this run (e.g. p3+ when migration is off, or
    simply not configured) contribute nothing when sampling, so it's safe
    to leave any of them out. Returns None (falls back to the pooled
    snapshotSamples cap in TreeOptions/Snapshot.txt) if none are set.
    """
    sample_cells = {}
    for key, value in params.items():
        m = re.fullmatch(r"Pop(\d+)SampleCells", key)
        if m:
            sample_cells[int(m.group(1))] = int(value)
    return sample_cells or None


def main():
    # This file lives in TumorGrowth/, alongside the .slim model, config.txt,
    # slim_config.txt, and the run outputs — everything for this scenario is
    # self-contained in one folder.
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # One-time cleanup: older versions of this script (and the now-retired
    # generate_migration_slim.py) wrote a permanent TG_Migration.slim into
    # this folder. This version never writes one — migration runs build
    # their script in memory and use an auto-deleted temp file instead — so
    # any TG_Migration.slim still sitting here is a leftover from that old
    # approach. Removed automatically on every run so it doesn't need to be
    # cleaned up by hand; every run's own log.txt (see out_dir below) still
    # records exactly what was built and run, so nothing is lost by not
    # keeping this file around.
    stale_migration_path = os.path.join(script_dir, "TG_Migration.slim")
    if os.path.exists(stale_migration_path):
        os.remove(stale_migration_path)
        print(f"[INFO] Removed stale leftover from the old approach: {stale_migration_path}")

    config   = su.load_config(os.path.join(script_dir, "slim_config.txt"))
    slim_exe = config.get("SLIM_EXE", "").strip()
    if not slim_exe:
        print("[ERROR] slim_config.txt must define SLIM_EXE.")
        sys.exit(1)
    if sys.platform == "win32" and not slim_exe.lower().endswith(".exe"):
        slim_exe += ".exe"

    args     = sys.argv[1:]
    cfg_name = next((a for a in args if a.endswith(".txt")), "config.txt")

    txt_path = os.path.join(script_dir, cfg_name)
    params   = su.load_config(txt_path)

    # Terminal overrides for this run only — config.txt itself, and its
    # migrationEnabled=F default, are never modified. `migration` is
    # shorthand for migrationEnabled=T; any other KEY=VALUE overrides that
    # one config.txt key (e.g. Pop2Time=30) for this run.
    cli_overrides = {}
    if "migration" in args:
        cli_overrides["migrationEnabled"] = "T"
    for a in args:
        if a == "migration" or a.endswith(".txt"):
            continue
        if "=" in a:
            k, _, v = a.partition("=")
            cli_overrides[k.strip()] = v.strip()
    if cli_overrides:
        print(f"[INFO] CLI overrides: {cli_overrides}")
        params.update(cli_overrides)

    migration_enabled = params.get("migrationEnabled", "F").strip().upper() in ("T", "TRUE", "1")
    pop_ids = detect_migration_pops(params) if migration_enabled else []

    # Set below only for migration runs — the temp file this points to is
    # always removed in the `finally` block near the end of this function,
    # whether the run below succeeds, fails, or raises.
    temp_slim_path = None

    if migration_enabled:
        if not pop_ids:
            print("[ERROR] migrationEnabled=T but config.txt has no PopNTime keys (N >= 2).")
            sys.exit(1)
        # Built fresh from the current TG.slim every run, so the migration
        # variant can never silently drift out of sync with the backbone.
        migration_text = build_migration_slim_text(script_dir, pop_ids)
        tmp = tempfile.NamedTemporaryFile(
            mode="w", suffix=".slim", prefix="TG_Migration_",
            delete=False, encoding="utf-8",
        )
        tmp.write(migration_text)
        tmp.close()
        slim_path = tmp.name
        slim_name = os.path.basename(slim_path)
        temp_slim_path = slim_path
        print(f"[INFO] Migration script built in memory from TG.slim + config.txt.")
        print(f"[INFO] Running it from a temporary file (deleted after this run): {slim_path}")
    else:
        slim_path = os.path.join(script_dir, SOURCE_SLIM)
        slim_name = SOURCE_SLIM
        if not os.path.exists(slim_path):
            print(f"[ERROR] Script not found: {slim_path}")
            sys.exit(1)

    # Separate output folders per mode (TG/Migration vs TG/NoMigration) — both
    # modes previously wrote to the same TG/ folder with the same filenames,
    # so running one after the other silently overwrote the first run's
    # trees/newick/snv/PNGs. log.txt still lives inside each mode's own
    # folder too, so each mode's run history is self-contained.
    out_dir    = os.path.join(script_dir, SCRIPT_NAME, "Migration" if migration_enabled else "NoMigration")
    os.makedirs(out_dir, exist_ok=True)
    stem       = os.path.splitext(params.get("treeOutputFile", SCRIPT_NAME))[0] or SCRIPT_NAME
    trees_path = os.path.join(out_dir, stem + ".trees").replace("\\", "/")
    log_path   = os.path.join(out_dir, "log.txt")

    d_args   = build_d_args(params, trees_path, migration_enabled, pop_ids)
    full_cmd = f"{slim_exe} " + " ".join(d_args) + f" {slim_name}"

    def _log(status, reason=None, slim_output="", outputs=None):
        with open(log_path, "a", encoding="utf-8") as lf:
            lf.write(f"\n=== simtreemaker.py Run: {datetime.datetime.now()} ===\n")
            lf.write(f"Script      : {slim_name}  (config: {cfg_name}, migration: {migration_enabled})\n")
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

    try:
        print(f"\n{'='*60}")
        print(f"  Script : {slim_name}  (config: {cfg_name}, migration: {migration_enabled})")
        print(f"  Output : {trees_path}")
        print(f"{'='*60}")
        print(f"[INFO] Equivalent terminal command:")
        print(f"    {full_cmd}\n")

        ok, slim_output = su.run_slim(slim_exe, slim_path, extra_args=d_args)
        if not ok:
            print(f"[SKIP] Skipping Newick/PNG conversion due to SLiM error.")
            _log("FAILED", reason="SLiM exited with a non-zero return code — see SLiM output above/below.",
                 slim_output=slim_output)
            return

        # -- Post-processing: .trees -> .nwk -> horizontal PNGs (labels + no labels) --
        tree_opts = su.load_tree_options(script_dir)
        mut_pos   = int(params.get("mutationPosition", 50000))

        dir_newick = os.path.join(out_dir, "newick")
        dir_png    = os.path.join(out_dir, "pngTree")
        dir_snv    = os.path.join(out_dir, "snv")
        os.makedirs(dir_newick, exist_ok=True)
        os.makedirs(dir_png, exist_ok=True)
        os.makedirs(dir_snv, exist_ok=True)
        nwk_path = os.path.join(dir_newick, stem + ".nwk")
        snv_path = os.path.join(dir_snv, stem + "_snv.csv")

        horizontal_only = ["horizontal_labels", "horizontal_no_labels"]
        sim_end_tick    = int(params["simulationEndTick"])
        sample_cells    = build_sample_cells(params)

        if not su.convert_to_newick(trees_path, nwk_path, mutation_position=mut_pos, tree_opts=tree_opts,
                                     snv_csv_path=snv_path, sim_end_tick=sim_end_tick,
                                     sample_cells=sample_cells):
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

    finally:
        # Runs whether the try block above returned early, succeeded, or
        # raised — the temporary migration script never lingers on disk.
        # Removed silently (no console line) — its path is already in the
        # "[INFO] Running it from a temporary file..." line printed above,
        # and in log.txt via full_cmd, if it's ever needed for reference.
        if temp_slim_path and os.path.exists(temp_slim_path):
            os.remove(temp_slim_path)


if __name__ == "__main__":
    main()
