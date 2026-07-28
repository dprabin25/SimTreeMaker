"""
simtreemaker.py  —  SLiM Cancer Simulation Runner
--------------------------------------------------
Run a specific model:
    python simtreemaker.py MutationSpread_NWF               # default: Lineage
    python simtreemaker.py MutationSpread_NWF Snapshot      # final-gen cells only
    python simtreemaker.py MutationSpread_NWF Lineage       # all generations, recency bias
    python simtreemaker.py MutationSpread_WF
    python simtreemaker.py ClonalGrowth_NWF
    python simtreemaker.py ClonalGrowth_WF
    python simtreemaker.py Metastasis_NWF
    python simtreemaker.py Metastasis_WF
    python simtreemaker.py CaseStudy
    python simtreemaker.py Tree

Tree sampling is controlled by TreeOptions/Snapshot.txt or TreeOptions/Lineage.txt.

Steps per model:
  1. Read slim_config.txt       (SLiM exe path)
  2. Read Options/<Name>.txt    (simulation parameters)
  3. Read TreeOptions/<Mode>.txt (tree sampling parameters)
  4. Generate + run a .slim script
  5. Convert .trees to .nwk     (requires: tskit)
  6. Save four PNG plots         (requires: biopython matplotlib)
  7. Write log.txt to output folder

Dependencies:
    conda create -n SimTreeMaker -c conda-forge python=3.11 biopython=1.85 matplotlib pyslim=1.1.1 tskit=1.0.3 r-base r-ape -y
"""

import sys
import os
import glob
import subprocess
import textwrap
import datetime


# -- 1. Parse slim_config.txt -------------------------------------------------

def load_config(config_path):
    config = {}
    if not os.path.exists(config_path):
        print(f"[ERROR] Config file not found: {config_path}")
        sys.exit(1)
    with open(config_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                key, _, value = line.partition("=")
                config[key.strip()] = value.strip()
    return config


# -- 2. Load model and tree-options .txt files --------------------------------

def load_model_txt(txt_path):
    """Read a key = value parameter file. Returns a list with one model dict."""
    if not os.path.exists(txt_path):
        print(f"[ERROR] Parameter file not found: {txt_path}")
        sys.exit(1)
    params = {}
    with open(txt_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                key, _, value = line.partition("=")
                value = value.split("#")[0]   # strip inline comments
                params[key.strip()] = value.strip()
    return [params]


def parse_valid_pops(vp_raw):
    """Parse valid_pops string → list of ints, or None for 'all'."""
    if vp_raw.strip().lower() == "all":
        return None
    return [int(x.strip()) for x in vp_raw.split(",") if x.strip()]


def load_tree_options(script_dir, mode="Lineage"):
    """
    Load TreeOptions/Snapshot.txt or TreeOptions/Lineage.txt.
    Returns a dict of parameters plus 'mode' key.
    """
    filename = f"{mode}.txt"
    path     = os.path.join(script_dir, "TreeOptions", filename)
    opts     = {"mode": mode}
    if not os.path.exists(path):
        print(f"[WARN] {path} not found — using defaults for {mode}.")
        return opts
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                key, _, value = line.partition("=")
                value = value.split("#")[0].strip()
                opts[key.strip()] = value
    print(f"[INFO] Tree options loaded: TreeOptions/{filename}")
    return opts


# -- 3. Generate SLiM scripts -------------------------------------------------

def generate_mutspread_nwf_script(p, trees_path):
    """MutationSpread nonWF — matches Sample 1 NWF.
    Driver (m2) injected into a newborn at mutationIntroTick.
    Cells older than maxAge are killed each generation.
    """
    convert = "F" if p["m2_convertToSubstitution"].strip() == "F" else "T"
    clones  = ("subpop.addCloned(individual); " * int(p["cloneCount"])).strip()
    return textwrap.dedent(f"""\
        initialize() {{
            initializeSLiMModelType("nonWF");
            initializeTreeSeq();
            initializeMutationRate({p['mutationRate']});
            initializeMutationType("m1", {p['m1_dominance']}, "f", {p['m1_effect']});
            initializeMutationType("m2", {p['m2_dominance']}, "f", {p['m2_effect']});
            m2.convertToSubstitution = {convert};
            initializeGenomicElementType("g1", m1, {p['m1_proportion']});
            initializeGenomicElement(g1, 0, {p['chromosomeEnd']});
            initializeRecombinationRate({p['recombinationRate']});
        }}

        1 early() {{
            sim.addSubpop("p1", {p['initialPopSize']});
        }}

        reproduction() {{
            {clones}
        }}

        early() {{
            targetSize = {p['targetSize']};
            popSize = p1.individualCount;
            p1.fitnessScaling = targetSize / popSize;
            for (ind in p1.individuals) {{
                if (ind.age >= {p['maxAge']}) {{
                    ind.fitnessScaling = 0.0;
                    sim.treeSeqRememberIndividuals(sim.subpopulations.individuals);
                }}
            }}
        }}

        {p['mutationIntroTick']} late() {{
            inds = p1.individuals;
            newborns = inds[inds.age == 0];
            if (size(newborns) > 0) {{
                ind = sample(newborns, 1);
                hapIndex = sample(0:1, 1);
                genome = ind.haplosomes[hapIndex];
                genome.addNewDrawnMutation(m2, {p['mutationPosition']});
            }}
        }}

        {p['simulationEndTick']} late() {{
            sim.treeSeqOutput("{trees_path}");
            catn("SLiM done. Tree written to: {trees_path}");
            sim.simulationFinished();
        }}
    """)


def generate_mutspread_wf_script(p, trees_path):
    """MutationSpread WF — matches Sample 1 WF.
    Driver (m2) seeded at mutationIntroTick, then manually spread
    each generation from mutationSpreadStart to mutationSpreadEnd.
    """
    convert = "F" if p["m2_convertToSubstitution"].strip() == "F" else "T"
    return textwrap.dedent(f"""\
        initialize() {{
            initializeTreeSeq();
            initializeMutationRate({p['mutationRate']});
            initializeMutationType("m1", {p['m1_dominance']}, "f", {p['m1_effect']});
            initializeMutationType("m2", {p['m2_dominance']}, "f", {p['m2_effect']});
            m2.convertToSubstitution = {convert};
            initializeGenomicElementType("g1", m1, {p['m1_proportion']});
            initializeGenomicElement(g1, 0, {p['chromosomeEnd']});
            initializeRecombinationRate({p['recombinationRate']});
        }}

        1 early() {{
            sim.addSubpop("p1", {p['targetSize']});
        }}

        {p['mutationIntroTick']} late() {{
            targetInd = sample(p1.individuals, 1);
            targetGenome = sample(targetInd.haplosomes, 2);
            m = targetGenome.addNewDrawnMutation(m2, {p['mutationPosition']});
            defineGlobal("seedMutation", m);
        }}

        {p['mutationSpreadStart']}:{p['mutationSpreadEnd']} late() {{
            targetInd = sample(p1.individuals, 1);
            targetGenome = sample(targetInd.haplosomes, 2);
            targetGenome.addMutations(seedMutation);
        }}

        {p['simulationEndTick']} late() {{
            sim.treeSeqOutput("{trees_path}");
            catn("SLiM done. Tree written to: {trees_path}");
            sim.simulationFinished();
        }}
    """)


def generate_clonalgrowth_nwf_script(p, trees_path):
    """ClonalGrowth nonWF — matches Sample 2 NWF.
    High m1_proportion, near-immortal cells (large maxAge).
    Driver injected directly (no newborn check) at mutationIntroTick.
    """
    convert = "F" if p["m2_convertToSubstitution"].strip() == "F" else "T"
    clones  = ("subpop.addCloned(individual); " * int(p["cloneCount"])).strip()
    return textwrap.dedent(f"""\
        initialize() {{
            initializeSLiMModelType("nonWF");
            initializeTreeSeq();
            initializeMutationRate({p['mutationRate']});
            initializeMutationType("m1", {p['m1_dominance']}, "f", {p['m1_effect']});
            m1.convertToSubstitution = F;
            initializeMutationType("m2", {p['m2_dominance']}, "f", {p['m2_effect']});
            m2.convertToSubstitution = {convert};
            initializeGenomicElementType("g1", m1, {p['m1_proportion']});
            initializeGenomicElement(g1, 0, {p['chromosomeEnd']});
            initializeRecombinationRate({p['recombinationRate']});
        }}

        1 early() {{
            sim.addSubpop("p1", {p['initialPopSize']});
        }}

        reproduction() {{
            {clones}
        }}

        early() {{
            targetSize = {p['targetSize']};
            popSize = p1.individualCount;
            p1.fitnessScaling = targetSize / popSize;
            for (ind in p1.individuals) {{
                if (ind.age >= {p['maxAge']}) {{
                    ind.fitnessScaling = 0.0;
                    sim.treeSeqRememberIndividuals(sim.subpopulations.individuals);
                }}
            }}
        }}

        {p['mutationIntroTick']} late() {{
            targetInd = sample(p1.individuals, 1);
            targetGenome = sample(targetInd.haplosomes, 2);
            targetGenome.addNewDrawnMutation(m2, {p['mutationPosition']});
        }}

        {p['simulationEndTick']} late() {{
            sim.treeSeqOutput("{trees_path}");
            catn("SLiM done. Tree written to: {trees_path}");
            sim.simulationFinished();
        }}
    """)


def generate_clonalgrowth_wf_script(p, trees_path):
    """ClonalGrowth WF — matches Sample 2 WF.
    Exponential growth of p1, population split to p2 at populationSplitTick,
    two driver mutations each with their own spread window.
    """
    convert = "F" if p["m2_convertToSubstitution"].strip() == "F" else "T"
    return textwrap.dedent(f"""\
        initialize() {{
            initializeTreeSeq();
            initializeMutationRate({p['mutationRate']});
            initializeMutationType("m1", {p['m1_dominance']}, "f", {p['m1_effect']});
            initializeMutationType("m2", {p['m2_dominance']}, "f", {p['m2_effect']});
            m2.convertToSubstitution = {convert};
            initializeGenomicElementType("g1", m1, {p['m1_proportion']});
            initializeGenomicElement(g1, 0, {p['chromosomeEnd']});
            initializeRecombinationRate({p['recombinationRate']});
            defineGlobal("growthAcc", {p['growthAcc']});
            defineGlobal("growthRate", {p['growthRate']});
        }}

        1 early() {{
            sim.addSubpop("p1", {p['initialPopSize']});
        }}

        {p['mutationIntroTick']} late() {{
            targetInd = sample(p1.individuals, 1);
            targetGenome = sample(targetInd.haplosomes, 2);
            m = targetGenome.addNewDrawnMutation(m2, {p['mutationPosition']});
            defineGlobal("seedMutation", m);
        }}

        {p['mutationSpreadStart']}:{p['mutationSpreadEnd']} late() {{
            targetInd = sample(p1.individuals, 1);
            targetGenome = sample(targetInd.haplosomes, 2);
            targetGenome.addMutations(seedMutation);
        }}

        {p['growthStartTick']}:{p['growthEndTick']} early() {{
            growthAcc = growthAcc + (p1.individualCount * (growthRate - 1));
            if (growthAcc >= 1.0) {{
                increase = asInteger(growthAcc);
                p1.setSubpopulationSize(p1.individualCount + increase);
                growthAcc = growthAcc - increase;
            }}
        }}

        {p['populationSplitTick']} early() {{
            sim.addSubpopSplit("p2", {p['splitPopSize']}, p1);
        }}

        {p['secondMutationTick']} late() {{
            targetInd = sample(p1.individuals, 1);
            targetGenome = sample(targetInd.haplosomes, 2);
            mm = targetGenome.addNewDrawnMutation(m2, {p['secondMutationPosition']});
            defineGlobal("seedMutation2", mm);
        }}

        {p['secondSpreadStart']}:{p['secondSpreadEnd']} late() {{
            targetInd = sample(p1.individuals, 1);
            targetGenome = sample(targetInd.haplosomes, 2);
            targetGenome.addMutations(seedMutation2);
        }}

        {p['pop2GrowthStartTick']}:{p['pop2GrowthEndTick']} early() {{
            growthAcc = growthAcc + (p2.individualCount * (growthRate - 1));
            if (growthAcc >= 1.0) {{
                increase = asInteger(growthAcc);
                p2.setSubpopulationSize(p2.individualCount + increase);
                growthAcc = growthAcc - increase;
            }}
        }}

        {p['simulationEndTick']} late() {{
            sim.treeSeqOutput("{trees_path}");
            catn("SLiM done. Tree written to: {trees_path}");
            sim.simulationFinished();
        }}
    """)


def generate_metastasis_nwf_script(p, trees_path):
    """Metastasis nonWF — matches Sample 3 NWF.
    Two populations: p1 (primary tumour) and p2 (metastatic site, starts empty).
    Driver introduced in p1, then cells migrate to p2 at migrationTick.
    """
    convert = "F" if p["m2_convertToSubstitution"].strip() == "F" else "T"
    clones  = ("subpop.addCloned(individual); " * int(p["cloneCount"])).strip()
    return textwrap.dedent(f"""\
        initialize() {{
            initializeSLiMModelType("nonWF");
            initializeTreeSeq();
            initializeMutationRate({p['mutationRate']});
            initializeMutationType("m1", {p['m1_dominance']}, "f", {p['m1_effect']});
            initializeMutationType("m2", {p['m2_dominance']}, "f", {p['m2_effect']});
            m2.convertToSubstitution = {convert};
            initializeGenomicElementType("g1", m1, {p['m1_proportion']});
            initializeGenomicElement(g1, 0, {p['chromosomeEnd']});
            initializeRecombinationRate({p['recombinationRate']});
        }}

        1 early() {{
            sim.addSubpop("p1", {p['initialPopSize']});
            sim.addSubpop("p2", {p['initialPop2Size']});
        }}

        reproduction() {{
            {clones}
        }}

        early() {{
            targetSize = {p['targetSize']};
            popSize = p1.individualCount;
            p1.fitnessScaling = targetSize / popSize;
            for (ind in p1.individuals) {{
                if (ind.age >= {p['maxAge']}) {{
                    ind.fitnessScaling = 0.0;
                    sim.treeSeqRememberIndividuals(sim.subpopulations.individuals);
                }}
            }}
        }}

        early() {{
            targetSize = {p['targetSize']};
            popSize = p2.individualCount;
            p2.fitnessScaling = (popSize > 0) ? targetSize / popSize else 1.0;
            for (ind in p2.individuals) {{
                if (ind.age >= {p['maxAge']}) {{
                    ind.fitnessScaling = 0.0;
                }}
            }}
        }}

        {p['mutationIntroTick']} late() {{
            inds = p1.individuals;
            newborns = inds[inds.age == 0];
            if (size(newborns) > 0) {{
                ind = sample(newborns, 1);
                hapIndex = sample(0:1, 1);
                genome = ind.haplosomes[hapIndex];
                genome.addNewDrawnMutation(m2, {p['mutationPosition']});
            }}
        }}

        {p['migrationTick']} early() {{
            if (p1.individualCount > 0) {{
                migrants = sample(p1.individuals, {p['migrantCount']});
                p2.takeMigrants(migrants);
            }}
        }}

        {p['simulationEndTick']} late() {{
            sim.treeSeqOutput("{trees_path}");
            catn("SLiM done. Tree written to: {trees_path}");
            sim.simulationFinished();
        }}
    """)


def generate_metastasis_wf_script(p, trees_path):
    """Metastasis WF — matches Sample 3 WF.
    Two populations from generation 1. Driver seeded and spread in p1.
    Migration rate from p1 → p2 turned on at migrationStartTick.
    """
    convert = "F" if p["m2_convertToSubstitution"].strip() == "F" else "T"
    return textwrap.dedent(f"""\
        initialize() {{
            initializeTreeSeq();
            initializeMutationRate({p['mutationRate']});
            initializeMutationType("m1", {p['m1_dominance']}, "f", {p['m1_effect']});
            initializeMutationType("m2", {p['m2_dominance']}, "f", {p['m2_effect']});
            m2.convertToSubstitution = {convert};
            initializeGenomicElementType("g1", m1, {p['m1_proportion']});
            initializeGenomicElement(g1, 0, {p['chromosomeEnd']});
            initializeRecombinationRate({p['recombinationRate']});
        }}

        1 early() {{
            sim.addSubpop("p1", {p['targetSize']});
            sim.addSubpop("p2", {p['initialPop2Size']});
        }}

        {p['mutationIntroTick']} late() {{
            targetInd = sample(p1.individuals, 1);
            targetGenome = sample(targetInd.haplosomes, 2);
            m = targetGenome.addNewDrawnMutation(m2, {p['mutationPosition']});
            defineGlobal("seedMutation", m);
        }}

        {p['mutationSpreadStart']}:{p['mutationSpreadEnd']} late() {{
            targetInd = sample(p1.individuals, 1);
            targetGenome = sample(targetInd.haplosomes, 2);
            targetGenome.addMutations(seedMutation);
        }}

        {p['migrationStartTick']}:{p['migrationEndTick']} early() {{
            p2.setMigrationRates(p1, {p['migrationRate']});
        }}

        {p['simulationEndTick']} late() {{
            sim.treeSeqOutput("{trees_path}");
            catn("SLiM done. Tree written to: {trees_path}");
            sim.simulationFinished();
        }}
    """)


def generate_script(p, trees_path):
    """Dispatch to the correct model-specific generator based on ModelName."""
    name = p.get("ModelName", "")
    if "MutationSpread" in name:
        if p["modelType"].strip() == "WF":
            return generate_mutspread_wf_script(p, trees_path)
        return generate_mutspread_nwf_script(p, trees_path)
    if "ClonalGrowth" in name:
        if p["modelType"].strip() == "WF":
            return generate_clonalgrowth_wf_script(p, trees_path)
        return generate_clonalgrowth_nwf_script(p, trees_path)
    if "Metastasis" in name:
        if p["modelType"].strip() == "WF":
            return generate_metastasis_wf_script(p, trees_path)
        return generate_metastasis_nwf_script(p, trees_path)
    # fallback: treat as MutationSpread
    if p["modelType"].strip() == "WF":
        return generate_mutspread_wf_script(p, trees_path)
    return generate_mutspread_nwf_script(p, trees_path)


# -- 5. Run SLiM --------------------------------------------------------------

def run_slim(slim_exe, script_path, extra_args=None):
    """Run SLiM, streaming output live to terminal and capturing it for the log.
    extra_args: list of additional flags e.g. ['-d', 'AGE=100', '-d', 'REP=1']
    """
    cmd = [slim_exe] + (extra_args or []) + [script_path]
    print(f"[INFO] Running SLiM: {' '.join(cmd)}")
    process = subprocess.Popen(
        cmd,  # -l 0 suppresses init chatter
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    output_lines = []
    for line in process.stdout:
        print(line, end="", flush=True)
        output_lines.append(line)
    process.wait()
    slim_output = "".join(output_lines)
    if process.returncode != 0:
        print(f"[ERROR] SLiM exited with code {process.returncode}")
        return False, slim_output
    print("[INFO] SLiM run complete.")
    return True, slim_output


# -- 5. Convert .trees to .nwk ------------------------------------------------

def convert_to_newick(trees_path, nwk_path, mutation_position=10000, tree_opts=None):
    """
    Route to Snapshot or Lineage newick builder based on tree_opts['mode'].
    tree_opts is loaded from TreeOptions/Snapshot.txt or TreeOptions/Lineage.txt.
    """
    if tree_opts is None:
        tree_opts = {"mode": "Lineage"}

    try:
        from slim_newick import build_newick, build_newick_snapshot
    except ImportError:
        print("[ERROR] slim_newick.py not found in the same folder as simtreemaker.py")
        return False

    mode = tree_opts.get("mode", "Lineage")

    try:
        seed_raw = tree_opts.get("seed", "none").strip().lower()
        seed     = None if seed_raw == "none" else int(seed_raw)

        if mode == "Snapshot":
            valid_pops       = parse_valid_pops(tree_opts.get("valid_pops", "all"))
            snapshot_samples = int(tree_opts.get("snapshotSamples", 400))
            stats = build_newick_snapshot(
                trees_path        = trees_path,
                nwk_path          = nwk_path,
                mutation_position = mutation_position,
                valid_pops        = valid_pops,
                snapshot_samples  = snapshot_samples,
                keep_unary        = False,
                seed              = seed,
            )
        else:  # Lineage
            valid_pops  = parse_valid_pops(tree_opts.get("valid_pops", "all"))
            max_per_gen = int(tree_opts.get("max_samples_per_generation", 5))
            min_per_gen = int(tree_opts.get("min_samples_per_generation", 1))
            stats = build_newick(
                trees_path                 = trees_path,
                nwk_path                   = nwk_path,
                mutation_position          = mutation_position,
                valid_pops                 = valid_pops,
                max_samples_per_generation = max_per_gen,
                min_samples_per_generation = min_per_gen,
                keep_unary                 = True,
                seed                       = seed,
            )

        print(f"[INFO] {stats['n_sampled']} nodes sampled | {stats['n_roots']} root(s)")
        return True

    except Exception as e:
        print(f"[ERROR] Newick conversion failed: {e}")
        return False


# -- 6. PNG plots -------------------------------------------------------------

import shutil
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def find_rscript():
    """Return path to Rscript if available, else None."""
    return shutil.which("Rscript")


def plot_png_r(nwk_path, dir_png, stem, model_name, model_type, script_dir, mode="Lineage"):
    """Render tree PNGs using R (ape). Returns True if successful."""
    rscript = find_rscript()
    if not rscript:
        return False
    r_script = os.path.join(script_dir, "plot_tree.R")
    if not os.path.exists(r_script):
        print("[WARN] plot_tree.R not found alongside simtreemaker.py")
        return False
    result = subprocess.run(
        [rscript, r_script, nwk_path, dir_png, stem, model_name, model_type, mode],
        capture_output=True, text=True
    )
    if result.stdout:
        print(result.stdout.rstrip())
    if result.stderr:
        print(result.stderr.rstrip())
    return result.returncode == 0


def _layout(clade, depths, counter, node_y):
    if not clade.clades:
        y = counter[0]
        node_y[id(clade)] = y
        counter[0] += 1
        return y
    child_ys = [_layout(c, depths, counter, node_y) for c in clade.clades]
    y = (min(child_ys) + max(child_ys)) / 2.0
    node_y[id(clade)] = y
    return y


def _draw_clade(clade, depths, max_depth, node_y, ax, show_labels, lw, color, orientation="horizontal"):
    depth = depths.get(clade, 0) / max_depth
    pos   = node_y[id(clade)]

    if orientation == "horizontal":
        x, y = depth, pos
        if clade.clades:
            child_ys = [node_y[id(c)] for c in clade.clades]
            ax.plot([x, x], [min(child_ys), max(child_ys)], color=color, lw=lw, solid_capstyle="round")
            for child in clade.clades:
                cx = depths.get(child, 0) / max_depth
                cy = node_y[id(child)]
                ax.plot([x, cx], [cy, cy], color=color, lw=lw, solid_capstyle="round")
                _draw_clade(child, depths, max_depth, node_y, ax, show_labels, lw, color, orientation)
        else:
            if show_labels:
                ax.text(x + 0.008, y, clade.name or "", va="center", ha="left",
                        fontsize=5.5, fontfamily="monospace", color="#222222")
    else:
        x, y = pos, depth
        if clade.clades:
            child_xs = [node_y[id(c)] for c in clade.clades]
            ax.plot([min(child_xs), max(child_xs)], [y, y], color=color, lw=lw, solid_capstyle="round")
            for child in clade.clades:
                cx = node_y[id(child)]
                cy = depths.get(child, 0) / max_depth
                ax.plot([cx, cx], [y, cy], color=color, lw=lw, solid_capstyle="round")
                _draw_clade(child, depths, max_depth, node_y, ax, show_labels, lw, color, orientation)
        else:
            if show_labels:
                ax.text(x, y + 0.018, clade.name or "", va="top", ha="center",
                        fontsize=5.5, fontfamily="monospace", color="#222222", rotation=90)


def _render_tree(bio_tree, ax, show_labels, title, subtitle, color="#2c5f8a", orientation="horizontal"):
    depths = bio_tree.depths(unit_branch_lengths=False)
    max_d  = max(depths.values()) if depths else 1
    if max_d == 0:
        depths = bio_tree.depths(unit_branch_lengths=True)
        max_d  = max(depths.values()) if depths else 1

    counter = [0]
    node_y  = {}
    _layout(bio_tree.root, depths, counter, node_y)
    n_tips = counter[0]

    _draw_clade(bio_tree.root, depths, max_d, node_y, ax,
                show_labels=show_labels, lw=0.7, color=color, orientation=orientation)

    pad = 0.3
    if orientation == "horizontal":
        ax.set_xlim(-0.02, 1.35 if show_labels else 1.05)
        ax.set_ylim(-pad, n_tips - 1 + pad)
        ax.invert_yaxis()
    else:
        ax.set_xlim(-pad, n_tips - 1 + pad)
        ax.set_ylim(-0.05, 1.22 if show_labels else 1.05)
        ax.invert_yaxis()

    ax.axis("off")
    ax.set_title(title, fontsize=10, fontweight="bold", pad=8, color="#111111")
    if subtitle:
        ax.text(0.5, -0.02, subtitle, transform=ax.transAxes,
                ha="center", fontsize=6.5, color="#888888", style="italic")


def plot_png(nwk_path, dir_png, stem, model_name, model_type, mode="Lineage"):
    try:
        from Bio import Phylo
    except ImportError:
        print("[ERROR] Missing biopython.")
        return

    trees = list(Phylo.parse(nwk_path, "newick"))
    if not trees:
        print("[ERROR] No trees found in Newick file.")
        return

    tree    = max(trees, key=lambda t: t.count_terminals())
    n_tips  = tree.count_terminals()
    n_trees = len(trees)

    subtitle = (
        f"Largest of {n_trees} subtrees · {n_tips} tips  (increase simulationEndTick for full coalescence)"
        if n_trees > 1 else f"{n_tips} tips"
    )
    tag = f"{model_name}  [{model_type}]  |  {mode}"
    h   = max(5, n_tips * 0.16)
    w   = max(8, n_tips * 0.16)

    variants = [
        ("horizontal", True,  (13, h), f"{stem}_horizontal_labels.png"),
        ("horizontal", False, (10, h), f"{stem}_horizontal_no_labels.png"),
        ("vertical",   True,  (w, 10), f"{stem}_vertical_labels.png"),
        ("vertical",   False, (w, 8),  f"{stem}_vertical_no_labels.png"),
    ]

    for orientation, show_labels, figsize, fname in variants:
        fig, ax = plt.subplots(figsize=figsize, dpi=150)
        fig.patch.set_facecolor("white")
        _render_tree(tree, ax,
                     show_labels=show_labels,
                     title=f"{tag} - {orientation.capitalize()} | {'With' if show_labels else 'No'} Labels",
                     subtitle=subtitle,
                     orientation=orientation)
        out = os.path.join(dir_png, fname)
        fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        print(f"[INFO] Saved: {out}")


# -- 7. Run one model ---------------------------------------------------------

def write_log(log_path, model, summary_lines, slim_cmd="", slim_output="", slim_script=""):
    with open(log_path, "a", encoding="utf-8") as lf:
        lf.write(f"\n=== Run: {datetime.datetime.now()} ===\n")
        lf.write(f"Model : {model.get('ModelName')} [{model.get('modelType')}]\n\n")
        if slim_cmd:
            slim_exe_part = slim_cmd.split()[0]
            script_name   = os.path.basename(slim_cmd.split()[-1])
            d_flags       = " ".join(f"-d {k}={v}" for k, v in model.items())
            lf.write(f"Command:\n  {slim_exe_part} {script_name} {d_flags}\n\n")
        if slim_output:
            lf.write("SLiM output:\n")
            lf.write(slim_output)
            lf.write("\n")
        lf.write("Outputs:\n")
        for line in summary_lines:
            lf.write(line + "\n")


def run_model(model, slim_exe, work_dir, script_dir, tree_opts=None, tree_only=False):
    if tree_opts is None:
        tree_opts = {"mode": "Lineage"}

    model_name = model["ModelName"]
    model_type = model["modelType"]
    stem       = os.path.splitext(model["treeOutputFile"])[0]
    mode       = tree_opts.get("mode", "Lineage")

    base       = os.path.join(work_dir, model_name, mode)
    dir_tree   = os.path.join(base, "tree")
    dir_newick = os.path.join(base, "newick")
    dir_png    = os.path.join(base, "pngTree")
    for d in (base, dir_tree, dir_newick, dir_png):
        os.makedirs(d, exist_ok=True)

    trees_path = os.path.join(dir_tree,   stem + ".tree").replace("\\", "/")
    nwk_path   = os.path.join(dir_newick, stem + ".nwk")
    log_path   = os.path.join(base, "log.txt")

    print(f"\n{'='*60}")
    print(f"  Model       : {model_name}  [{model_type}]")
    print(f"  Tree mode   : {mode}")
    print(f"  Step        : {'Tree only (skipping SLiM)' if tree_only else 'Full run'}")
    print(f"  Folder      : {base}")
    print(f"{'='*60}")

    slim_cmd    = ""
    slim_output = ""
    slim_script = ""

    if not tree_only:
        # Generate .slim and run SLiM
        sim_options_dir  = os.path.join(script_dir, "SimOptions")
        slim_script_path = os.path.join(sim_options_dir, f"{model_name}.slim")
        slim_script = generate_script(model, trees_path)
        with open(slim_script_path, "w") as f:
            f.write(slim_script)

        slim_script_rel  = os.path.relpath(slim_script_path, script_dir).replace("\\", "/")
        slim_cmd = f"{slim_exe} {slim_script_rel}"
        print(f"[INFO] SLiM script : {slim_script_rel}")
        print(f"[INFO] SLiM command: {slim_cmd}")
        print(f"\nParameters:")
        max_k = max(len(k) for k in model)
        for k, v in model.items():
            print(f"  {k:<{max_k}}  =  {v}")
        if tree_opts:
            print(f"\nTree options ({mode}):")
            to = {k: v for k, v in tree_opts.items() if k != "mode"}
            max_k2 = max(len(k) for k in to) if to else 1
            for k, v in to.items():
                print(f"  {k:<{max_k2}}  =  {v}")
        print()

        ok, slim_output = run_slim(slim_exe, slim_script_path)
        if not ok:
            print(f"[SKIP] Skipping post-processing for {stem} due to SLiM error.")
            return
    else:
        if not os.path.exists(trees_path):
            print(f"[ERROR] No .tree file found at: {trees_path}")
            print(f"  Run without 'tree' flag first to generate it.")
            return
        print(f"[INFO] Using existing: {trees_path}")

    mut_pos = int(model.get("mutationPosition", 10000))
    if not convert_to_newick(trees_path, nwk_path,
                             mutation_position=mut_pos,
                             tree_opts=tree_opts):
        return

    r_ok = plot_png_r(nwk_path, dir_png, stem, model_name, model_type, script_dir, mode)
    if not r_ok:
        print("[INFO] Rscript unavailable — using matplotlib.")
        plot_png(nwk_path, dir_png, stem, model_name, model_type, mode)

    rel = os.path.relpath(base, work_dir).replace("\\", "/")
    summary = [
        f"  {rel}/tree/    {stem}.tree",
        f"  {rel}/newick/  {stem}.nwk",
        f"  {rel}/pngTree/ {stem}_horizontal_labels.png / _no_labels.png",
        f"  {rel}/pngTree/ {stem}_vertical_labels.png   / _no_labels.png",
        f"  {rel}/log.txt",
    ]
    print(f"\n[DONE] Outputs in: {base}")
    for line in summary:
        print(line)

    write_log(log_path, model, summary, slim_cmd=slim_cmd, slim_output=slim_output, slim_script=slim_script)


# -- 8. CaseStudy runner ------------------------------------------------------

import re as _re

def _extract_slim_defs(slim_file):
    """
    Parse the first line of a .slim script for inline -d constants.
    Format:  //slim -d AGE=100 -d REP=1
    Returns a list like ['-d', 'AGE=100', '-d', 'REP=1'].
    """
    params = []
    try:
        with open(slim_file) as f:
            first = f.readline().strip()
        if first.startswith("//slim") or first.startswith("// slim"):
            tokens = first.lstrip("/").split()
            i = 0
            while i < len(tokens):
                if tokens[i] == "-d" and i + 1 < len(tokens):
                    params += ["-d", tokens[i + 1]]
                    i += 2
                else:
                    i += 1
    except Exception:
        pass
    return params


def _parse_tree_output_path(slim_file, work_dir):
    """
    Read the .slim script and extract the path from sim.treeSeqOutput("...").
    Returns absolute path, or None if not found.
    """
    try:
        with open(slim_file, encoding="utf-8", errors="ignore") as f:
            content = f.read()
        m = _re.search(r'treeSeqOutput\s*\(\s*["\']([^"\']+)["\']', content)
        if m:
            raw = m.group(1)
            # resolve relative to work_dir (where SLiM will be run from)
            return os.path.normpath(os.path.join(work_dir, raw))
    except Exception:
        pass
    return None


def run_casestudy(work_dir, slim_exe, script_dir, tree_opts=None, tree_only=False, cli_defs=None):
    """
    Run every .slim file in CaseStudy/.
    SLiM writes its .tree wherever the script says.
    Python then applies Snapshot or Lineage on that .tree file.
    Outputs go to CaseStudyOutputs/<stem>/<Mode>/
    """
    if tree_opts is None:
        tree_opts = {"mode": "Snapshot"}

    mode   = tree_opts.get("mode", "Snapshot")
    cs_dir = os.path.join(work_dir, "CaseStudy")

    if not os.path.isdir(cs_dir):
        print(f"[ERROR] CaseStudy folder not found: {cs_dir}")
        return

    slim_files = sorted(glob.glob(os.path.join(cs_dir, "*.slim")))
    if not slim_files:
        print(f"[ERROR] No .slim files found in {cs_dir}")
        return

    # Create CaseStudyOutputs/ immediately so it always exists
    cs_out = os.path.join(work_dir, "CaseStudyOutputs")
    os.makedirs(cs_out, exist_ok=True)

    print(f"\n[CaseStudy] Found {len(slim_files)} script(s) in {cs_dir}")
    print(f"[CaseStudy] Tree mode : {mode}")
    print(f"[CaseStudy] Outputs   : CaseStudyOutputs/")
    print(f"[CaseStudy] Step      : {'Tree only (skipping SLiM)' if tree_only else 'Full run'}\n")

    def _find_tree_files(root):
        """Return set of all .tree / .trees files under root."""
        found = set()
        for dirpath, _, filenames in os.walk(root):
            for fn in filenames:
                if fn.endswith(".tree") or fn.endswith(".trees"):
                    found.add(os.path.join(dirpath, fn))
        return found

    for slim_file in slim_files:
        stem = os.path.splitext(os.path.basename(slim_file))[0]

        # Merge -d defs: script first-line defaults, overridden by CLI
        script_defs = {}
        for pair in _extract_slim_defs(slim_file):
            if "=" in pair:
                k, _, v = pair.partition("=")
                script_defs[k.upper()] = v
        merged_defs = {**script_defs, **(cli_defs or {})}

        # Build extra_args list for subprocess
        extra_args = []
        for k, v in merged_defs.items():
            extra_args += ["-d", f"{k}={v}"]

        # Include def values in output folder name for traceability
        def_suffix = "_".join(f"{k}{v}" for k, v in merged_defs.items())
        named_stem = f"{stem}_{def_suffix}" if def_suffix else stem

        slim_cmd    = f"{slim_exe} CaseStudy/{stem}.slim"
        slim_output = ""

        # Output folders (created before SLiM so they're always ready)
        base       = os.path.join(cs_out, named_stem, mode)
        dir_newick = os.path.join(base, "newick")
        dir_png    = os.path.join(base, "pngTree")
        for d in (base, dir_newick, dir_png):
            os.makedirs(d, exist_ok=True)
        log_path = os.path.join(base, "log.txt")

        print(f"\n{'='*60}")
        print(f"  Script  : CaseStudy/{stem}.slim")
        if merged_defs:
            for k, v in merged_defs.items():
                print(f"  -d {k:<8}: {v}")
        print(f"  Outputs : CaseStudyOutputs/{named_stem}/{mode}/")
        print(f"{'='*60}")

        # ── Run SLiM ──────────────────────────────────────────────────────
        if not tree_only:
            defs_str = " ".join(f"-d {k}={v}" for k, v in merged_defs.items())
            print(f"[INFO] SLiM command: {slim_exe} {defs_str + ' ' if defs_str else ''}CaseStudy/{stem}.slim")

            slim_start = datetime.datetime.now().timestamp()

            ok, slim_output = run_slim(slim_exe, slim_file, extra_args=extra_args)
            if not ok:
                print(f"[SKIP] SLiM failed for {stem}.slim")
                continue

            # Find the .tree file written/updated during this SLiM run
            all_trees = _find_tree_files(work_dir)
            recent    = [f for f in all_trees if os.path.getmtime(f) >= slim_start]
            if not recent:
                # fallback: most recently modified tree file overall
                recent = sorted(all_trees, key=os.path.getmtime, reverse=True)
            if not recent:
                print(f"[ERROR] SLiM ran but no .tree file was found.")
                continue
            trees_path = recent[0]
            print(f"[INFO] Tree written: {os.path.relpath(trees_path, work_dir).replace(chr(92), '/')}")

        else:
            # tree_only: find the most recent .tree in work_dir matching stem name
            candidates = [f for f in _find_tree_files(work_dir)
                          if stem.lower() in os.path.basename(f).lower()]
            if not candidates:
                # fallback: any .tree anywhere
                candidates = sorted(_find_tree_files(work_dir),
                                    key=os.path.getmtime, reverse=True)
            if not candidates:
                print(f"[ERROR] No .tree file found. Run without 'tree' flag first.")
                continue
            trees_path = candidates[0]
            print(f"[INFO] Using existing: {os.path.relpath(trees_path, work_dir).replace(chr(92), '/')}")

        nwk_path = os.path.join(dir_newick, f"{stem}.nwk")

        # ── Newick conversion ─────────────────────────────────────────────
        if not convert_to_newick(trees_path, nwk_path,
                                 mutation_position=10000,
                                 tree_opts=tree_opts):
            continue

        # ── PNG plots ─────────────────────────────────────────────────────
        r_ok = plot_png_r(nwk_path, dir_png, stem, stem, "CaseStudy", script_dir, mode)
        if not r_ok:
            plot_png(nwk_path, dir_png, stem, stem, mode)

        rel = os.path.relpath(base, work_dir).replace("\\", "/")
        tree_rel = os.path.relpath(trees_path, work_dir).replace("\\", "/")
        summary = [
            f"  {tree_rel}  (.tree — written by SLiM script)",
            f"  {rel}/newick/  {stem}.nwk",
            f"  {rel}/pngTree/ {stem}_horizontal_labels.png / _no_labels.png",
            f"  {rel}/pngTree/ {stem}_vertical_labels.png   / _no_labels.png",
            f"  {rel}/log.txt",
        ]
        print(f"\n[DONE] {stem}")
        for line in summary:
            print(line)

        with open(log_path, "a", encoding="utf-8") as lf:
            d_flags  = " ".join(f"-d {k}={v}" for k, v in merged_defs.items())
            full_cmd = f"{slim_exe} {stem}.slim{' ' + d_flags if d_flags else ''}"
            lf.write(f"\n=== CaseStudy Run: {datetime.datetime.now()} ===\n")
            lf.write(f"Script      : CaseStudy/{stem}.slim\n")
            lf.write(f"SLiM command: {full_cmd}\n")
            lf.write(f"Tree mode   : {mode}\n")
            if slim_output:
                lf.write("SLiM output:\n")
                lf.write(slim_output)
            lf.write("\nOutputs:\n")
            for line in summary:
                lf.write(line + "\n")

    print(f"\n{'='*60}")
    print(f"  CaseStudy complete.")
    print(f"{'='*60}")


# -- 9. ReadyTrees ------------------------------------------------------------

def run_readytrees(script_dir, max_samples_per_generation=5, min_samples_per_generation=1):
    in_dir  = os.path.join(script_dir, "ReadyTrees")
    out_dir = os.path.join(script_dir, "ReadyTreesOutputs")

    if not os.path.isdir(in_dir):
        os.makedirs(in_dir, exist_ok=True)
        print(f"[ReadyTrees] Created input folder: {in_dir}")
        print(f"[ReadyTrees] Drop your .trees files there, then re-run.")
        return

    tree_files = sorted(
        glob.glob(os.path.join(in_dir, "*.trees")) +
        glob.glob(os.path.join(in_dir, "*.tree"))
    )
    if not tree_files:
        print(f"[ReadyTrees] No .trees files found in {in_dir}")
        return

    print(f"\n[ReadyTrees] Found {len(tree_files)} file(s) in {in_dir}")
    print(f"[ReadyTrees] Sampling: max {max_samples_per_generation}/gen, min {min_samples_per_generation}/gen (recency bias)")
    print(f"[ReadyTrees] Outputs -> {out_dir}\n")

    for trees_path in tree_files:
        stem     = os.path.splitext(os.path.basename(trees_path))[0]
        run_out  = os.path.join(out_dir, f"{stem}Output")
        dir_nwk  = os.path.join(run_out, "newick")
        dir_png  = os.path.join(run_out, "pngTree")
        log_path = os.path.join(run_out, "log.txt")
        for d in (dir_nwk, dir_png):
            os.makedirs(d, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"  File: {os.path.basename(trees_path)}")
        print(f"  Out : {run_out}")
        print(f"{'='*60}")

        nwk_path = os.path.join(dir_nwk, f"{stem}.nwk")
        if not convert_to_newick(
            trees_path, nwk_path,
            max_samples_per_generation = max_samples_per_generation,
            min_samples_per_generation = min_samples_per_generation,
        ):
            continue

        plot_png(nwk_path, dir_png, stem, stem, "ReadyTree")

        summary = [
            f"  newick/  {stem}.nwk",
            f"  pngTree/ {stem}_horizontal_labels.png / _no_labels.png",
            f"           {stem}_vertical_labels.png   / _no_labels.png",
            f"  log.txt",
        ]
        print(f"[DONE] {stem}")
        for line in summary:
            print(line)

        with open(log_path, "a", encoding="utf-8") as lf:
            lf.write(f"\n=== ReadyTrees Run: {datetime.datetime.now()} ===\n")
            lf.write(f"Input: {os.path.basename(trees_path)}\n")
            lf.write(f"max_samples_per_generation: {max_samples_per_generation}\n")
            lf.write(f"min_samples_per_generation: {min_samples_per_generation}\n")
            for line in summary:
                lf.write(line + "\n")

    print(f"\n{'='*60}")
    print(f"  ReadyTrees complete. Outputs in: {out_dir}")
    print(f"{'='*60}")


# -- Main ---------------------------------------------------------------------

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    config     = load_config(os.path.join(script_dir, "slim_config.txt"))
    slim_exe   = config.get("SLIM_EXE", "").strip()

    if not slim_exe:
        print("[ERROR] slim_config.txt must define SLIM_EXE.")
        sys.exit(1)

    # Ensure .exe on Windows
    if sys.platform == "win32" and not slim_exe.lower().endswith(".exe"):
        slim_exe += ".exe"

    work_dir = script_dir

    if len(sys.argv) < 2:
        print("Usage: python simtreemaker.py <ModelName|CaseStudy|Tree> [Snapshot|Lineage] [tree]")
        print("")
        print("  ModelName   — run a built-in model (reads SimOptions/<Name>.txt)")
        print("  CaseStudy   — run all .slim scripts in the CaseStudy/ folder")
        print("  Tree        — build tree from .tree files in ReadyTrees/")
        print("")
        print("  Snapshot    — sample cells at final generation only  (default)")
        print("  Lineage     — sample cells across all generations with recency bias")
        print("  tree        — skip SLiM, rebuild tree from existing .tree file")
        print("")
        print("  Examples:")
        print("    python simtreemaker.py MutationSpread_NWF")
        print("    python simtreemaker.py MutationSpread_NWF Snapshot")
        print("    python simtreemaker.py MutationSpread_NWF Lineage")
        print("    python simtreemaker.py MutationSpread_NWF Snapshot tree")
        print("    python simtreemaker.py CaseStudy Snapshot")
        print("    python simtreemaker.py CaseStudy Snapshot AGE=100 REP=1")
        print("    python simtreemaker.py CaseStudy Lineage  AGE=50  REP=3")
        print("    python simtreemaker.py CaseStudy Snapshot tree")
        print("    python simtreemaker.py Tree")
        sys.exit(1)

    arg       = sys.argv[1]
    args      = sys.argv[2:]
    tree_mode = next((a for a in args if a in ("Snapshot", "Lineage")), "Snapshot")
    tree_only = "tree" in args

    # Parse KEY=VALUE overrides from command line (e.g. AGE=100 REP=2)
    cli_defs = {}
    for a in args:
        if "=" in a and a not in ("Snapshot", "Lineage") and not a.startswith("-"):
            k, _, v = a.partition("=")
            cli_defs[k.strip().upper()] = v.strip()

    if arg == "CaseStudy":
        run_casestudy(work_dir, slim_exe, script_dir,
                      tree_opts=load_tree_options(script_dir, tree_mode),
                      tree_only=tree_only,
                      cli_defs=cli_defs)
        return

    if arg == "Tree":
        run_readytrees(script_dir)
        return

    tree_opts   = load_tree_options(script_dir, tree_mode)
    options_dir = os.path.join(work_dir, "SimOptions")
    txt_name    = arg if arg.endswith(".txt") else arg + ".txt"
    txt_path    = os.path.join(options_dir, txt_name)
    if not os.path.exists(txt_path):
        print(f"[ERROR] Parameter file not found: {txt_path}")
        print(f"  Place your .txt parameter files in: SimOptions/")
        sys.exit(1)

    print(f"[INFO] Using: {txt_path}")
    models = load_model_txt(txt_path)
    print(f"[INFO] {os.path.basename(txt_path)}: {len(models)} model(s)")
    for model in models:
        run_model(model, slim_exe, work_dir, script_dir, tree_opts, tree_only=tree_only)

    print(f"\n{'='*60}")
    print(f"  All runs complete. Outputs in: {work_dir}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
