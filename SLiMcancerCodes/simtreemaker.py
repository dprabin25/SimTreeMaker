"""
simtreemaker.py — runs TG.slim directly for no-migration. For migration
runs, this file itself builds the migration variant (in memory, from
TG.slim + config.txt — no separate generator module) and runs it from a
short-lived temporary .slim file that is deleted immediately after SLiM
finishes. Nothing named TG_Migration.slim (or similar) is ever saved
permanently in this folder.

TG.slim is a plain no-migration backbone with no migration switch and
no p2/p3 code at all — there's no way to turn migration on by typing extra
-d flags at TG.slim itself. Migration is Python-only:

  1. Directly from the terminal, supplying every -d constant yourself —
     this always means no migration, since TG.slim has no migration branch
     to turn on:

       slim -d initialPopSize=1 -d targetSize=150 \
            -d mutationRate=1e-8 -d driverMutationRate=1e-8 \
            -d m1Start=0 -d m1End=49999 -d m2Start=50000 -d m2End=99999 \
            -d m2_dominance=0.5 -d m2_effect=0.6 \
            -d simulationEndTick=200 -d treesPath="'TG.trees'" TG.slim

  2. Through this script, which reads slim_config.txt (for the SLiM
     executable path) and config.txt (for the parameter values) instead of
     you typing the -d flags by hand:

       python simtreemaker.py

     The mode is chosen on the command line only. Plain
     `python simtreemaker.py` (or `nomigration`) runs TG.slim, p1 only.
     `migration` builds the migration variant fresh from the current
     TG.slim + config.txt's PopN blocks, writes it to a temporary .slim
     file, runs that instead of TG.slim, then deletes the temp file:

       python simtreemaker.py nomigration
       python simtreemaker.py migration

     Add `noise` to either mode to finish with add_noise.py on the clean
     SNV FASTA, using config.txt's Sequencing noise settings (output in
     <outputDir>/TG/<mode>/noise/):

       python simtreemaker.py migration noise

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
slim_utilities.load_tree_options() instead.

Every run also appends an entry to <outputDir>/TG/<mode>/log.txt — the SLiM command used,
SLiM's own console output, and the resulting output paths (or, on failure,
the reason) — so a run can be audited after the fact without needing to
keep terminal scrollback around.
"""

import os
import re
import sys
import csv
import shutil
import tempfile
import datetime

import slim_utilities as su
import add_noise

SCRIPT_NAME = "TG"
SOURCE_SLIM = f"{SCRIPT_NAME}.slim"   # "TG.slim"
TREATMENT_SLIM = f"{SCRIPT_NAME}_Treatment.slim"   # used with `nomigration treatment`


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


BASE_FLAGS = (
    "//   slim -d initialPopSize=1 -d targetSize=150 \\\n"
    "//        -d mutationRate=1e-8 -d driverMutationRate=1e-8 \\\n"
    "//        -d m1Start=0 -d m1End=49999 -d m2Start=50000 -d m2End=99999 \\\n"
    "//        -d m2_dominance=0.5 -d m2_effect=0.6 \\\n"
    "//        -d simulationEndTick=200 -d treesPath=\"'TG.trees'\" \\"
)
N_BASE_FLAGS = 12


def build_header(pop_ids, source=SOURCE_SLIM):
    pop_list = ", ".join(f"p{n}" for n in pop_ids)
    example_flags = " \\\n".join(
        f"//        -d Pop{n}Time=50 -d Pop{n}IniPopSize=5 -d Pop{n}TargetSize=100 "
        f"-d Pop{n}MigrationRate=0.01 -d Pop{n}Origin=1"
        for n in pop_ids
    )
    n_flags = N_BASE_FLAGS + 5 * len(pop_ids)
    return f"""// TG_Migration variant — built in memory by simtreemaker.py from {source}
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
{BASE_FLAGS}
{example_flags} \\
//        <this temporary file>
//
// Running this file directly in SLiMgui (or via `slim` with no -d flags)
// will fail with "undefined identifier" errors — that's expected; supply
// all {n_flags} -d flags listed above."""


def insert_migration_block(text, pop_ids, source=SOURCE_SLIM):
    """Replace the inert 3-line marker inside early() with the real
    migration Eidos code. Raises ValueError if the marker isn't found, so
    a drifted/hand-edited base file fails loudly instead of silently
    producing a script with no migration logic in it."""
    if MARKER not in text:
        raise ValueError(
            f"Marker not found in {source} — expected exactly:\n{MARKER!r}\n"
            "(base script may have been edited since this generator was written)"
        )
    return text.replace(MARKER, build_migration_block(pop_ids), 1)


def rewrite_header(text, pop_ids, source=SOURCE_SLIM):
    """Replace the backbone's leading comment block (everything after the
    "//SLiM version" line, up to initialize()) with the migration header.
    Matching on structure rather than exact text means editing the
    backbone's header comments never breaks migration builds."""
    first_nl = text.find("\n")
    init_at  = text.find("\ninitialize()")
    if not text.startswith("//SLiM version") or init_at < 0:
        raise ValueError(f"{source} must start with '//SLiM version' and contain initialize().")
    return text[:first_nl + 1] + "//\n" + build_header(pop_ids, source) + "\n" + text[init_at:]


def build_migration_slim_text(script_dir, pop_ids, source=SOURCE_SLIM):
    """
    Build the full migration-variant Eidos script text in memory, from the
    current backbone (TG.slim) +
    the given pop_ids. Returns a string — the caller is responsible for
    writing it to a temporary file and cleaning that file up afterward (see
    main()). Read as UTF-8 explicitly, since the Windows default encoding
    is not UTF-8.
    """
    src_path = os.path.join(script_dir, source)
    with open(src_path, encoding="utf-8") as f:
        text = f.read()
    text = rewrite_header(text, pop_ids, source)
    text = insert_migration_block(text, pop_ids, source)
    return text


POP_COLUMNS = ("Origin", "Time", "IniPopSize", "TargetSize", "MigrationRate", "SampleCells")


def expand_pop_rows(params):
    """
    Expand config.txt's migration table rows (`Pop2 = Pop1 50 5 100 0.01 100`)
    into the individual PopN<Column> keys the rest of this script uses
    (Pop2Origin, Pop2Time, ...), in POP_COLUMNS order. Done before CLI
    overrides are applied, so `Pop2Time=30` on the command line still
    overrides a single cell of the table.
    """
    for key in [k for k in params if re.fullmatch(r"Pop\d+", k)]:
        values = params.pop(key).split()
        if len(values) != len(POP_COLUMNS):
            print(f"[ERROR] {key} row has {len(values)} value(s); expected {len(POP_COLUMNS)}: "
                  f"{' '.join(POP_COLUMNS)}.")
            sys.exit(1)
        for column, value in zip(POP_COLUMNS, values):
            params[f"{key}{column}"] = value


def detect_migration_pops(params):
    """
    Scan config.txt's parsed params for any PopNTime key (N >= 2) and
    return the sorted list of population numbers found — this is what
    makes the number of migration populations config-driven instead of
    hardcoded. Add a 5th population by adding a `Pop5 = ...` row to
    config.txt's migration table; no code changes needed anywhere in this
    file.
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


def parse_regions(params, names=("m1", "m2")):
    """
    Read and validate the genome layout (m1Start/m1End/m2Start/m2End, plus
    m3Start/m3End for treatment runs) from config.txt. Returns
    {"m1": (start, end), ...}, inclusive. Exits with a clear message on a
    bad layout so SLiM never starts with one. The .slim scripts repeat
    these checks for direct terminal runs.
    """
    regions = {}
    for m in names:
        bounds = []
        for edge in ("Start", "End"):
            key = f"{m}{edge}"
            raw = params.get(key)
            if raw is None:
                print(f"[ERROR] config.txt is missing {key}.")
                sys.exit(1)
            if not re.fullmatch(r"\d+", raw):
                print(f"[ERROR] {key}={raw!r} must be a whole number >= 0.")
                sys.exit(1)
            bounds.append(int(raw))
        if bounds[0] > bounds[1]:
            print(f"[ERROR] {m}Start ({bounds[0]}) must be <= {m}End ({bounds[1]}).")
            sys.exit(1)
        regions[m] = tuple(bounds)
    for i, p in enumerate(names):
        for q in names[i + 1:]:
            (a1, b1), (a2, b2) = regions[p], regions[q]
            if a1 <= b2 and a2 <= b1:
                print(f"[ERROR] {p} [{a1}-{b1}] and {q} [{a2}-{b2}] overlap.")
                sys.exit(1)
    return regions


def parse_treatment(params):
    """
    Read and validate config.txt's Treatment section (only with
    `nomigration treatment`). Returns the extra -d flags TG_Treatment.slim
    needs, as {name: value}.
    """
    def need(key):
        if key not in params:
            print(f"[ERROR] 'treatment' requested but config.txt is missing {key}.")
            sys.exit(1)
        return params[key].strip()

    end_tick = int(params["simulationEndTick"])
    tick = need("treatmentTick")
    if not re.fullmatch(r"\d+", tick) or not 1 < int(tick) < end_tick:
        print(f"[ERROR] treatmentTick={tick!r} must be a whole number between 2 and "
              f"simulationEndTick - 1 ({end_tick - 1}).")
        sys.exit(1)
    try:
        kill = add_noise.check_rate(need("drugKill"), "drugKill")
        cost = add_noise.check_rate(need("resistanceCost"), "resistanceCost")
        rate = float(need("resistanceMutationRate"))
        if rate < 0:
            raise ValueError("resistanceMutationRate must be >= 0.")
    except ValueError as e:
        print(f"[ERROR] {e}")
        sys.exit(1)

    return {"treatmentTick": tick, "drugKill": kill, "resistanceCost": cost,
            "resistanceMutationRate": params["resistanceMutationRate"].strip(),
            "m3Start": params["m3Start"], "m3End": params["m3End"]}


NOISE_KEYS = {          # config.txt key -> add_noise() argument
    "missingRate": "missing",
    "fpRate":      "fp",
    "fnRate":      "fn",
    "cellSpread":  "cell_spread",
    "siteSpread":  "site_spread",
}


def parse_noise(params):
    """Read and validate config.txt's Sequencing noise section (only when the
    `noise` keyword is given). Checked before SLiM runs so a bad value never
    costs a full simulation."""
    noise = {}
    for key, arg in NOISE_KEYS.items():
        if key not in params:
            print(f"[ERROR] 'noise' requested but config.txt is missing {key}.")
            sys.exit(1)
        try:
            noise[arg] = add_noise.check_rate(params[key], key,
                                              inclusive=not key.endswith("Spread"))
        except ValueError as e:
            print(f"[ERROR] {e}")
            sys.exit(1)
    seed = params.get("noiseSeed", "").strip()
    noise["seed"] = int(seed) if seed else None
    return noise


def control_matches_clean(control_fasta, clean_fasta):
    """True if the zero-noise control holds exactly the clean cells (same
    names, order and sequences), plus only the added Normal record."""
    clean_names, clean = add_noise.read_fasta(clean_fasta)
    ctrl_names, ctrl   = add_noise.read_fasta(control_fasta)
    return (ctrl_names == clean_names + ["Normal"]
            and all(ctrl[n] == clean[n] for n in clean_names))


def summarise_trajectory(traj_csv, dose_tick):
    """Plain-language summary of a treatment run's trajectory CSV: size at
    the dose, recovery time, and resistant share before/after/at the end."""
    if not os.path.exists(traj_csv):
        return ["trajectory file missing"]
    total, res = {}, {}
    with open(traj_csv) as fh:
        for r in csv.DictReader(fh):
            g = int(r["generation"])
            total[g], res[g] = int(r["total"]), int(r["resistant"])
    gens = sorted(total)
    pct  = lambda g: 100.0 * res[g] / total[g] if total.get(g) else 0.0
    before, at = total.get(dose_tick - 1), total.get(dose_tick)
    if before is None or at is None:
        return ["the run ended before the dose"]
    back = next((g - dose_tick for g in gens if g > dose_tick and total[g] >= 0.95 * before), None)
    end  = gens[-1]
    return [
        f"dose at generation {dose_tick}: {before} -> {at} cells",
        f"back to pre-dose size after {back} generations" if back else
        f"not back to pre-dose size by generation {end}",
        f"resistant cells: {pct(dose_tick - 1):.1f}% before -> {pct(dose_tick):.1f}% after the dose "
        f"-> {pct(end):.1f}% at generation {end}",
    ]


def warn_config(params, regions):
    """Warn (non-fatal) when a region is so small, or its rate so low, that
    fewer than one mutation is expected over the whole run."""
    target = int(params["targetSize"].split(",")[0])
    ticks  = int(params["simulationEndTick"])
    # ~target new cells per tick, 2 genome copies each
    copies = 2 * target * ticks
    for m, rate_key in (("m1", "mutationRate"), ("m2", "driverMutationRate")):
        rate = float(params[rate_key])
        length = regions[m][1] - regions[m][0] + 1
        expected = copies * length * rate
        if rate > 0 and expected < 1:
            print(f"[WARN] {m}: ~{expected:.2g} {m} mutations expected over the run "
                  f"({length} bp x {rate_key}={rate:g}); raise the rate or widen the region.")


def build_d_args(params, trees_path, migration_enabled, pop_ids):
    """
    Build the -d flag list, from config.txt's parsed params dict — the
    same flags you'd type by hand if running SLiM directly from the
    terminal. No migration switch is ever passed — TG.slim doesn't
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
    for key in ("m1Start", "m1End", "m2Start", "m2End"):
        _d(key, params[key])
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
    snapshotSamples cap in slim_utilities.load_tree_options()) if none are set.
    """
    sample_cells = {}
    for key, value in params.items():
        m = re.fullmatch(r"Pop(\d+)SampleCells", key)
        if m:
            sample_cells[int(m.group(1))] = int(value)
    return sample_cells or None


def main():
    # This file lives in SLiMcancerCodes/, alongside the .slim models,
    # config.txt and slim_config.txt; run outputs go to config.txt's outputDir.
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
    if not os.path.isfile(slim_exe) and not shutil.which(slim_exe):
        print(f"[ERROR] SLiM not found at SLIM_EXE = {slim_exe}")
        print("        Set SLIM_EXE in slim_config.txt to your SLiM program; find it with")
        print("        `which slim` (macOS/Linux) or `where slim` (Windows) in the SimTreeMaker environment.")
        sys.exit(1)

    args     = sys.argv[1:]
    cfg_name = next((a for a in args if a.endswith(".txt")), "config.txt")

    txt_path = os.path.join(script_dir, cfg_name)
    params   = su.load_config(txt_path)
    expand_pop_rows(params)

    # Run mode: `migration` or `nomigration` (the default). Any KEY=VALUE
    # overrides that one config.txt key (e.g. Pop2Time=30) for this run
    # only — config.txt itself is never modified.
    if "migration" in args and "nomigration" in args:
        print("[ERROR] Use either 'migration' or 'nomigration', not both.")
        sys.exit(1)
    migration_enabled = "migration" in args
    if "noise" in args and "nonoise" in args:
        print("[ERROR] Use either 'noise' or 'nonoise', not both.")
        sys.exit(1)
    if "treatment" in args and "notreatment" in args:
        print("[ERROR] Use either 'treatment' or 'notreatment', not both.")
        sys.exit(1)
    if "treatment" in args and migration_enabled:
        print("[ERROR] 'treatment' is for a single population: use `python simtreemaker.py nomigration treatment`.")
        sys.exit(1)

    cli_overrides = {}
    for a in args:
        if a in ("migration", "nomigration", "noise", "nonoise", "treatment", "notreatment") or a.endswith(".txt"):
            continue
        if "=" in a:
            k, _, v = a.partition("=")
            cli_overrides[k.strip()] = v.strip()
        else:
            print(f"[ERROR] Unknown argument {a!r} (expected migration, nomigration, noise, nonoise, treatment, "
                  f"notreatment, a .txt config file, or KEY=VALUE).")
            sys.exit(1)
    if cli_overrides:
        print(f"[INFO] CLI overrides: {cli_overrides}")
        params.update(cli_overrides)

    # Noise: the `noise` / `nonoise` keyword wins; otherwise config.txt's
    # `noise = on/off` decides (both modes).
    noise_switch = params.get("noise", "off").strip().lower()
    if noise_switch not in ("on", "off"):
        print(f"[ERROR] noise={noise_switch!r} must be 'on' or 'off'.")
        sys.exit(1)
    if "noise" in args:
        noise_enabled = True
    elif "nonoise" in args:
        noise_enabled = False
    else:
        noise_enabled = noise_switch == "on"

    # Treatment: the `treatment` / `notreatment` keyword wins; otherwise
    # config.txt's `treatment = on/off` decides, for no-migration runs only.
    switch = params.get("treatment", "off").strip().lower()
    if switch not in ("on", "off"):
        print(f"[ERROR] treatment={switch!r} must be 'on' or 'off'.")
        sys.exit(1)
    if "treatment" in args:
        treatment_enabled = True
    elif "notreatment" in args:
        treatment_enabled = False
    else:
        treatment_enabled = switch == "on" and not migration_enabled
        if switch == "on" and migration_enabled:
            print("[INFO] treatment = on is ignored for migration runs (treatment is no-migration only).")
    # Treatment runs use their own (larger) carrying capacity, so normal runs
    # and the manuscript batch keep targetSize unchanged.
    if treatment_enabled and params.get("treatmentTargetSize", "").strip():
        params["targetSize"] = params["treatmentTargetSize"].strip()
        print(f"[INFO] Treatment run: targetSize = treatmentTargetSize = {params['targetSize']}")

    pop_ids = detect_migration_pops(params) if migration_enabled else []
    regions = parse_regions(params, ("m1", "m2", "m3") if treatment_enabled else ("m1", "m2"))
    warn_config(params, regions)
    noise   = parse_noise(params) if noise_enabled else None
    treatment = parse_treatment(params) if treatment_enabled else None
    backbone  = TREATMENT_SLIM if treatment_enabled else SOURCE_SLIM

    # Set below only for migration runs — the temp file this points to is
    # always removed in the `finally` block near the end of this function,
    # whether the run below succeeds, fails, or raises.
    temp_slim_path = None

    if migration_enabled:
        if not pop_ids:
            print("[ERROR] 'migration' requested but config.txt has no PopNTime keys (N >= 2).")
            sys.exit(1)
        end_tick = int(params["simulationEndTick"])
        late = [n for n in pop_ids if int(params[f"Pop{n}Time"]) >= end_tick]
        if late:
            print(f"[ERROR] Pop{late[0]}Time ({params[f'Pop{late[0]}Time']}) must be earlier than "
                  f"simulationEndTick ({end_tick}), or that population is never created.")
            sys.exit(1)
        # Built fresh from the current TG.slim every run, so the migration
        # variant can never silently drift out of sync with the backbone.
        migration_text = build_migration_slim_text(script_dir, pop_ids, backbone)
        tmp = tempfile.NamedTemporaryFile(
            mode="w", suffix=".slim", prefix="TG_Migration_",
            delete=False, encoding="utf-8",
        )
        tmp.write(migration_text)
        tmp.close()
        slim_path = tmp.name
        slim_name = os.path.basename(slim_path)
        temp_slim_path = slim_path
        print(f"[INFO] Migration script built in memory from {backbone} + config.txt.")
        print(f"[INFO] Running it from a temporary file (deleted after this run): {slim_path}")
    else:
        slim_path = os.path.join(script_dir, backbone)
        slim_name = backbone
        if not os.path.exists(slim_path):
            print(f"[ERROR] Script not found: {slim_path}")
            sys.exit(1)

    # Separate output folders per mode (TG/Migration vs TG/NoMigration) — both
    # modes previously wrote to the same TG/ folder with the same filenames,
    # so running one after the other silently overwrote the first run's
    # trees/newick/snv/PNGs. log.txt still lives inside each mode's own
    # folder too, so each mode's run history is self-contained.
    # outputDir (config.txt) is relative to this script's folder unless absolute;
    # left blank, outputs go next to the script as before.
    mode_dir   = "Treatment" if treatment_enabled else ("Migration" if migration_enabled else "NoMigration")
    out_root   = os.path.join(script_dir, params.get("outputDir", "").strip())
    out_dir    = os.path.normpath(os.path.join(out_root, SCRIPT_NAME, mode_dir))
    os.makedirs(out_dir, exist_ok=True)
    stem       = os.path.splitext(params.get("treeOutputFile", SCRIPT_NAME))[0] or SCRIPT_NAME
    trees_path = os.path.join(out_dir, stem + ".trees").replace("\\", "/")
    log_path   = os.path.join(out_dir, "log.txt")
    dir_traj   = os.path.join(out_dir, "trajectory")
    traj_csv   = os.path.join(dir_traj, stem + "_trajectory.csv")

    d_args   = build_d_args(params, trees_path, migration_enabled, pop_ids)
    if treatment:
        os.makedirs(dir_traj, exist_ok=True)
        for name, value in treatment.items():
            d_args.extend(["-d", f"{name}={value}"])
        d_args.extend(["-d", f"trajectoryPath='{traj_csv.replace(chr(92), '/')}'"])
    sim_seed = params.get("simSeed", "").strip()
    if sim_seed:
        if not re.fullmatch(r"\d+", sim_seed):
            print(f"[ERROR] simSeed={sim_seed!r} must be a whole number >= 0 (or empty for random).")
            sys.exit(1)
        d_args = ["-s", sim_seed] + d_args
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
        # SLiM prints the seed it used ("// Initial random seed:" + number);
        # read it back so random-seed runs can still be reproduced exactly.
        m = re.search(r"Initial random seed:\s*(\d+)", slim_output)
        used_seed = m.group(1) if m else (sim_seed or "unknown")
        if not ok:
            print(f"[SKIP] Skipping Newick/PNG conversion due to SLiM error.")
            _log("FAILED", reason="SLiM exited with a non-zero return code — see SLiM output above/below.",
                 slim_output=slim_output)
            return

        treatment_summary = []
        if treatment:
            treatment_summary = summarise_trajectory(traj_csv, int(treatment["treatmentTick"]))
            su.plot_trajectory_r(traj_csv, os.path.join(dir_traj, stem + "_trajectory.png"),
                                 int(treatment["treatmentTick"]), script_dir)
            erad = re.search(r"TUMOR ERADICATED at generation (\d+)", slim_output)
            if erad:
                # No cells left to sample: clear this mode's old tree/SNV/FASTA
                # outputs so they can't be mistaken for this run's.
                for sub in ("newick", "pngTree", "snv", "fasta", "noise"):
                    shutil.rmtree(os.path.join(out_dir, sub), ignore_errors=True)
                if os.path.exists(trees_path):
                    os.remove(trees_path)
                reason = f"Tumor eradicated at generation {erad.group(1)}; no cells to sample."
                print(f"\n[RESULT] {reason}")
                _log("ERADICATED", reason=reason, slim_output=slim_output,
                     outputs=[f"  {traj_csv}", *(f"  {l}" for l in treatment_summary)])
                print(f"  Trajectory: {traj_csv}")
                return

        # -- Post-processing: .trees -> .nwk -> horizontal PNGs (labels + no labels) --
        tree_opts = su.load_tree_options(script_dir)
        # Seed the cell sampling with the same seed SLiM used, so one seed
        # reproduces the whole run (tumor + which cells are sampled).
        if used_seed.isdigit():
            tree_opts["seed"] = used_seed

        dir_newick = os.path.join(out_dir, "newick")
        dir_png    = os.path.join(out_dir, "pngTree")
        dir_snv    = os.path.join(out_dir, "snv")
        os.makedirs(dir_newick, exist_ok=True)
        os.makedirs(dir_png, exist_ok=True)
        os.makedirs(dir_snv, exist_ok=True)
        nwk_path = os.path.join(dir_newick, stem + ".nwk")
        snv_path = os.path.join(dir_snv, stem + "_snv.csv")
        dir_fasta    = os.path.join(out_dir, "fasta")
        os.makedirs(dir_fasta, exist_ok=True)
        fasta_prefix = os.path.join(dir_fasta, stem)
        fasta_seed   = int(params.get("fastaSeed", 12345))
        write_full   = params.get("writeFullFasta", "F").strip().upper() in ("T", "TRUE", "1")

        horizontal_only = ["horizontal_labels", "horizontal_no_labels"]
        sim_end_tick    = int(params["simulationEndTick"])
        sample_cells    = build_sample_cells(params)

        if not su.convert_to_newick(trees_path, nwk_path, m2_range=regions["m2"], tree_opts=tree_opts,
                                     snv_csv_path=snv_path, sim_end_tick=sim_end_tick,
                                     sample_cells=sample_cells,
                                     fasta_prefix=fasta_prefix, fasta_seed=fasta_seed,
                                     write_full_fasta=write_full):
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
            f"  SLiM seed: {used_seed}",
            f"  {trees_path}",
            f"  {nwk_path}",
            f"  {snv_path}",
            f"  {fasta_prefix}_parent.fasta",
            *([f"  {fasta_prefix}_full.fasta"] if write_full else []),
            f"  {fasta_prefix}_snv.fasta",
            f"  {fasta_prefix}_snv_positions.tsv",
            f"  {dir_png}/{stem}_horizontal_labels.png",
            f"  {dir_png}/{stem}_horizontal_no_labels.png",
        ]

        # -- Optional last step: sequencing noise on the clean SNV FASTA --
        # Old noise files describe the previous run's cells, so clear them
        # every run (like every other output here, they are regenerated).
        dir_noise = os.path.join(out_dir, "noise")
        if os.path.isdir(dir_noise):
            shutil.rmtree(dir_noise)
        if noise:
            # Two outputs: control (every rate 0, so the same step, format and
            # Normal record, but no noise) and test (config.txt's rates).
            # Analyse both the same way; any difference is due to the noise.
            clean_fasta = fasta_prefix + "_snv.fasta"
            control = dict(noise, missing=0.0, fp=0.0, fn=0.0, cell_spread=0.0, site_spread=0.0)
            try:
                print("[noise] -- control (no noise) --")
                control_paths, control_stats = add_noise.add_noise(
                    clean_fasta, normal_name="Normal", out_dir=os.path.join(dir_noise, "control"), **control)
                print("[noise] -- test (config.txt rates) --")
                test_paths, test_stats = add_noise.add_noise(
                    clean_fasta, normal_name="Normal", out_dir=os.path.join(dir_noise, "test"), **noise)
            except ValueError as e:
                print(f"[ERROR] Noise step failed: {e}")
                _log("FAILED", reason=f"Noise step failed: {e}", slim_output=slim_output, outputs=outputs)
                return
            if not control_matches_clean(control_paths[0], clean_fasta):
                print("[ERROR] Noise control check FAILED: control FASTA differs from the clean SNV FASTA.")
                _log("FAILED", reason="Noise control FASTA differs from the clean SNV FASTA.",
                     slim_output=slim_output, outputs=outputs)
                return
            print("[noise] Control check passed: control FASTA = clean SNV FASTA (+ Normal).")
            outputs.extend(f"  {p}" for p in (*control_paths, *test_paths))
            noise_summary = [
                f"control check PASSED (control = clean SNV FASTA + Normal)",
                *(f"{name:<8}missing {s['missing']:>6}  FN {s['FN']:>5}  FP {s['FP']:>5}   "
                  f"({s['cells']} cells x {s['sites']} sites)"
                  for name, s in (("control", control_stats), ("test", test_stats))),
            ]
            outputs.extend(["", "Noise:", *(f"  {line}" for line in noise_summary)])
        if treatment:
            outputs.extend([f"  {traj_csv}", "", "Treatment:", *(f"  {l}" for l in treatment_summary)])
        _log("OK", slim_output=slim_output, outputs=outputs)

        print(f"\n{'='*60}")
        print(f"  Run complete.")
        print(f"  Seed   : {used_seed}  ({'simSeed' if sim_seed else 'random; set simSeed=' + used_seed + ' to repeat'})")
        print(f"  Tree   : {trees_path}")
        print(f"  Newick : {nwk_path}")
        print(f"  SNV CSV: {snv_path}")
        print(f"  FASTA  : {dir_fasta}  (parent, {'full, ' if write_full else ''}snv, snv_positions)")
        print(f"  PNGs   : {dir_png}/{stem}_horizontal_labels.png / _no_labels.png")
        if noise:
            print(f"  Noise  : {dir_noise}  (control/ = no noise, test/ = config.txt rates)")
            for line in noise_summary:
                print(f"           {line}")
        if treatment:
            print(f"  Treated: {dir_traj}  (trajectory CSV + PNG)")
            for line in treatment_summary:
                print(f"           {line}")
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
