"""
simtreemaker.py  —  SLiM Cancer Simulation Runner
----------------------------------------------
Run a specific model by CSV name:
    python simtreemaker.py MutationSpread
    python simtreemaker.py ClonalGrowth
    python simtreemaker.py Metastasis
    python simtreemaker.py CaseStudy
    python simtreemaker.py Tree

Steps per model row:
  1. Read slim_config.txt  (SLiM exe path + working directory)
  2. Read models.csv       (parameters per model)
  3. Generate + run a .slim script
  4. Convert .trees → .nwk  (requires: pip install tskit)
  5. Save two PNG plots     (requires: pip install biopython matplotlib)

Dependencies:
    pip install tskit biopython matplotlib
"""

import sys
import os
import csv
import glob
import subprocess
import textwrap


# ── 1. Parse slim_config.txt ─────────────────────────────────────────────────

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


# ── 2. Load models from models.csv ───────────────────────────────────────────

def load_all_models(models_path):
    if not os.path.exists(models_path):
        print(f"[ERROR] Models file not found: {models_path}")
        sys.exit(1)
    rows = []
    with open(models_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({k: v.strip() for k, v in row.items()})
    return rows


# ── 3. Generate SLiM scripts (separate templates for WF / nonWF) ─────────────

def generate_nonwf_script(p, trees_path):
    convert = "F" if p["m2_convertToSubstitution"] == "F" else "T"
    return textwrap.dedent(f"""\
        initialize() {{
            initializeSLiMModelType("nonWF");
            initializeTreeSeq();
            initializeMutationRate({p['mutationRate']});
            initializeMutationType("m1", {p['m1_dominance']}, "f", {p['m1_effect']});
            initializeMutationType("m2", {p['m2_dominance']}, "f", {p['m2_effect']});
            m2.convertToSubstitution = {convert};
            initializeGenomicElementType("g1", m1, {p['g1_mutationFraction']});
            initializeGenomicElement(g1, 0, {p['chromosomeEnd']});
            initializeRecombinationRate({p['recombinationRate']});
        }}

        1 early() {{
            sim.addSubpop("p1", {p['initialPopSize']});
        }}

        reproduction() {{
            {('subpop.addCloned(individual); ' * int(p['cloneCount'])).strip()}
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
            sim.simulationFinished();
        }}
    """)


def generate_wf_script(p, trees_path):
    """WF models: fixed pop size (= targetSize), no reproduction() block, no age-based death."""
    convert = "F" if p["m2_convertToSubstitution"] == "F" else "T"
    return textwrap.dedent(f"""\
        initialize() {{
            initializeSLiMModelType("WF");
            initializeTreeSeq();
            initializeMutationRate({p['mutationRate']});
            initializeMutationType("m1", {p['m1_dominance']}, "f", {p['m1_effect']});
            initializeMutationType("m2", {p['m2_dominance']}, "f", {p['m2_effect']});
            m2.convertToSubstitution = {convert};
            initializeGenomicElementType("g1", m1, {p['g1_mutationFraction']});
            initializeGenomicElement(g1, 0, {p['chromosomeEnd']});
            initializeRecombinationRate({p['recombinationRate']});
        }}

        1 early() {{
            // WF: use targetSize if set, otherwise initialPopSize
            sim.addSubpop("p1", {p['targetSize'] if int(p['targetSize']) > 0 else p['initialPopSize']});
        }}

        {p['mutationIntroTick']} late() {{
            // WF: ind.age not available - sample any individual directly
            ind = sample(p1.individuals, 1);
            genome = ind.haplosomes[sample(0:1, 1)];
            genome.addNewDrawnMutation(m2, {p['mutationPosition']});
        }}

        {p['simulationEndTick']} late() {{
            sim.treeSeqOutput("{trees_path}");
            sim.simulationFinished();
        }}
    """)


def generate_script(p, trees_path):
    if p["modelType"].strip() == "WF":
        return generate_wf_script(p, trees_path)
    return generate_nonwf_script(p, trees_path)


# ── 4. Run SLiM ──────────────────────────────────────────────────────────────

def run_slim(slim_exe, script_path):
    print(f"[INFO] Running SLiM: {slim_exe} {script_path}")
    result = subprocess.run([slim_exe, script_path], capture_output=True, text=True)
    if result.stdout:
        print(result.stdout)
    if result.stderr:
        print(result.stderr)
    if result.returncode != 0:
        print(f"[ERROR] SLiM exited with code {result.returncode}")
        return False
    print("[INFO] SLiM run complete.")
    return True


# ── 5. Convert .trees → .nwk (generation-sampled, mutation-labelled) ─────────

def convert_to_newick(trees_path, nwk_path, mutation_position=10000):
    try:
        from slim_newick import build_newick
    except ImportError:
        print("[ERROR] slim_newick.py not found in the same folder as simtreemaker.py")
        return False
    try:
        stats = build_newick(
            trees_path=trees_path,
            nwk_path=nwk_path,
            mutation_position=mutation_position,
            max_samples=400,
        )
        print(f"[INFO] {stats['n_sampled']} nodes sampled | "
              f"{stats['n_roots']} root(s)")
        return True
    except Exception as e:
        print(f"[ERROR] Newick conversion failed: {e}")
        return False


# ── 6. PNG plots — custom rectangular cladogram renderer ─────────────────────

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def _layout(clade, depths, counter, node_y):
    """Recursively assign y-positions. Tips get sequential integers; internals get midpoint."""
    if not clade.clades:                        # terminal
        y = counter[0]
        node_y[id(clade)] = y
        counter[0] += 1
        return y
    child_ys = [_layout(c, depths, counter, node_y) for c in clade.clades]
    y = (min(child_ys) + max(child_ys)) / 2.0
    node_y[id(clade)] = y
    return y


def _draw_clade(clade, depths, max_depth, node_y, ax, show_labels, lw, color,
                orientation="horizontal"):
    """Draw lines for this clade recursively. orientation: 'horizontal' or 'vertical'."""
    depth = depths.get(clade, 0) / max_depth
    pos   = node_y[id(clade)]          # sequential tip position

    if orientation == "horizontal":
        x, y = depth, pos
        if clade.clades:
            child_ys = [node_y[id(c)] for c in clade.clades]
            ax.plot([x, x], [min(child_ys), max(child_ys)],
                    color=color, lw=lw, solid_capstyle="round")
            for child in clade.clades:
                cx = depths.get(child, 0) / max_depth
                cy = node_y[id(child)]
                ax.plot([x, cx], [cy, cy],
                        color=color, lw=lw, solid_capstyle="round")
                _draw_clade(child, depths, max_depth, node_y, ax,
                            show_labels, lw, color, orientation)
        else:
            if show_labels:
                ax.text(x + 0.008, y, clade.name or "",
                        va="center", ha="left", fontsize=5.5,
                        fontfamily="monospace", color="#222222")
    else:  # vertical — root at top, tips at bottom
        x, y = pos, depth             # swap: sequential → x-axis, depth → y-axis
        if clade.clades:
            child_xs = [node_y[id(c)] for c in clade.clades]
            ax.plot([min(child_xs), max(child_xs)], [y, y],
                    color=color, lw=lw, solid_capstyle="round")
            for child in clade.clades:
                cx = node_y[id(child)]
                cy = depths.get(child, 0) / max_depth
                ax.plot([cx, cx], [y, cy],
                        color=color, lw=lw, solid_capstyle="round")
                _draw_clade(child, depths, max_depth, node_y, ax,
                            show_labels, lw, color, orientation)
        else:
            if show_labels:
                ax.text(x, y + 0.018, clade.name or "",
                        va="top", ha="center", fontsize=5.5,
                        fontfamily="monospace", color="#222222",
                        rotation=90)


def _render_tree(bio_tree, ax, show_labels, title, subtitle, color="#2c5f8a",
                 orientation="horizontal"):
    """Full render of one Bio.Phylo tree onto ax."""
    depths   = bio_tree.depths(unit_branch_lengths=True)
    max_d    = max(depths.values()) if depths else 1
    counter  = [0]
    node_y   = {}
    _layout(bio_tree.root, depths, counter, node_y)
    n_tips = counter[0]

    _draw_clade(bio_tree.root, depths, max_d, node_y, ax,
                show_labels=show_labels, lw=0.7, color=color,
                orientation=orientation)

    pad = 0.3
    if orientation == "horizontal":
        ax.set_xlim(-0.02, 1.35 if show_labels else 1.05)
        ax.set_ylim(-pad, n_tips - 1 + pad)
        ax.invert_yaxis()
    else:  # vertical
        ax.set_xlim(-pad, n_tips - 1 + pad)
        ax.set_ylim(-0.05, 1.22 if show_labels else 1.05)
        ax.invert_yaxis()   # depth=0 (root) at visual top

    ax.axis("off")
    ax.set_title(title, fontsize=10, fontweight="bold", pad=8, color="#111111")
    if subtitle:
        ax.text(0.5, -0.02, subtitle, transform=ax.transAxes,
                ha="center", fontsize=6.5, color="#888888", style="italic")


def plot_png(nwk_path, dir_png, stem, model_name, model_type):
    """
    Saves 4 PNGs into dir_png (pngTree/):
      stem_horizontal_labels.png
      stem_horizontal_no_labels.png
      stem_vertical_labels.png
      stem_vertical_no_labels.png
    """
    try:
        from Bio import Phylo
    except ImportError:
        print("[ERROR] Missing biopython. Run: pip install biopython")
        return

    trees = list(Phylo.parse(nwk_path, "newick"))
    if not trees:
        print("[ERROR] No trees found in Newick file.")
        return

    tree    = max(trees, key=lambda t: t.count_terminals())
    n_tips  = tree.count_terminals()
    n_trees = len(trees)

    subtitle = (
        f"Largest of {n_trees} subtrees · {n_tips} tips  "
        f"(increase simulationEndTick for full coalescence)"
        if n_trees > 1 else f"{n_tips} tips"
    )
    tag    = f"{stem}  [{model_type}]"
    h      = max(5, n_tips * 0.16)   # height for horizontal plots
    w      = max(8, n_tips * 0.16)   # width  for vertical plots

    variants = [
        ("horizontal", True,  (13, h), f"{stem}_horizontal_labels.png"),
        ("horizontal", False, (10, h), f"{stem}_horizontal_no_labels.png"),
        ("vertical",   True,  (w, 10), f"{stem}_vertical_labels.png"),
        ("vertical",   False, (w, 8),  f"{stem}_vertical_no_labels.png"),
    ]

    for orientation, show_labels, figsize, fname in variants:
        fig, ax = plt.subplots(figsize=figsize, dpi=150)
        fig.patch.set_facecolor("white")
        label_str = "With Labels" if show_labels else "No Labels"
        orient_str = orientation.capitalize()
        _render_tree(tree, ax,
                     show_labels=show_labels,
                     title=f"{tag} - {orient_str} | {label_str}",
                     subtitle=subtitle,
                     orientation=orientation)
        out = os.path.join(dir_png, fname)
        fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        print(f"[INFO] Saved: {out}")




# ── 8. Run one model row ──────────────────────────────────────────────────────

def run_model(model, slim_exe, work_dir, script_dir):
    model_name = model["ModelName"]
    model_type = model["modelType"]
    stem       = os.path.splitext(model["treeOutputFile"])[0]   # e.g. "S1"

    # Subfolder layout: work_dir/<ModelName>/<stem>/scripts|tree|newick|pngTree
    base          = os.path.join(work_dir, model_name, stem)
    dir_scripts   = os.path.join(base, "scripts")
    dir_tree      = os.path.join(base, "tree")
    dir_newick    = os.path.join(base, "newick")
    dir_png       = os.path.join(base, "pngTree")
    for d in (dir_scripts, dir_tree, dir_newick, dir_png):
        os.makedirs(d, exist_ok=True)

    slim_script_path = os.path.join(dir_scripts, f"{stem}.slim")
    trees_path       = os.path.join(dir_tree,    stem + ".tree").replace("\\", "/")
    nwk_path         = os.path.join(dir_newick,  stem + ".nwk")

    print(f"\n{'='*60}")
    print(f"  Model : {model_name}  [{model_type}]")
    print(f"  Folder: {base}")
    print(f"{'='*60}")

    # Generate + run SLiM
    script = generate_script(model, trees_path)
    with open(slim_script_path, "w") as f:
        f.write(script)
    print(f"[INFO] Script written: {slim_script_path}")

    if not run_slim(slim_exe, slim_script_path):
        print(f"[SKIP] Skipping post-processing for {stem} due to SLiM error.")
        return

    # Newick conversion
    mut_pos = int(model.get("mutationPosition", 10000))
    if not convert_to_newick(trees_path, nwk_path, mutation_position=mut_pos):
        return

    # PNG plots — all 4 variants into pngTree/
    plot_png(nwk_path, dir_png, stem, model_name, model_type)

    print(f"\n[DONE] Outputs in: {base}")
    print(f"  scripts/  {stem}.slim")
    print(f"  tree/     {stem}.tree")
    print(f"  newick/   {stem}.nwk")
    print(f"  pngTree/  {stem}_horizontal_labels.png / _no_labels.png")
    print(f"            {stem}_vertical_labels.png   / _no_labels.png")
    print(f"  pngTree/          {stem}_vertical_no_labels.png")
    print(f"  {stem}_tree.html  <- open in browser")


# ── 8. Main ──────────────────────────────────────────────────────────────────

# ── CaseStudy runner ─────────────────────────────────────────────────────────

def _extract_slim_params(slim_file):
    """
    Parse the first comment line for default -d arguments.
    e.g.  //slim -d AGE=100 -d REP=1 CHIP2Env.slim
    Returns a list like ["-d", "AGE=100", "-d", "REP=1"]
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


def run_casestudy(work_dir, slim_exe):
    import shutil

    cs_dir = os.path.join(work_dir, "CaseStudy")
    if not os.path.isdir(cs_dir):
        print(f"[ERROR] CaseStudy folder not found: {cs_dir}")
        return

    slim_files = sorted(glob.glob(os.path.join(cs_dir, "*.slim")))
    if not slim_files:
        print(f"[ERROR] No .slim files found in {cs_dir}")
        return

    # ── Locate / promote Slim2tree3.py ───────────────────────────────────────
    slim2tree_py  = os.path.join(cs_dir, "Slim2tree3.py")
    slim2tree_txt = os.path.join(cs_dir, "Slim2tree3.py.txt")
    if not os.path.exists(slim2tree_py) and os.path.exists(slim2tree_txt):
        shutil.copy(slim2tree_txt, slim2tree_py)
        print(f"[INFO] Copied Slim2tree3.py.txt → Slim2tree3.py")
    use_slim2tree = os.path.exists(slim2tree_py)

    # ── Output structure: CaseStudy/CaseStudyTrees/ ───────────────────────────
    out_base   = os.path.join(cs_dir, "CaseStudyTrees")
    dir_tree   = os.path.join(out_base, "tree")
    dir_newick = os.path.join(out_base, "newick")
    dir_png = os.path.join(out_base, "pngTree")
    for d in (dir_tree, dir_newick, dir_png):
        os.makedirs(d, exist_ok=True)

    print(f"\n[CaseStudy] Found {len(slim_files)} script(s) in {cs_dir}")
    print(f"[CaseStudy] Outputs → {out_base}")
    print(f"[CaseStudy] Mode: {'Slim2tree3.py' if use_slim2tree else 'slim_newick.py'}\n")

    # ── Add slim_exe directory to PATH so os.system('slim ...') works ────────
    slim_dir = os.path.dirname(os.path.abspath(slim_exe))
    env = os.environ.copy()
    env["PATH"] = slim_dir + os.pathsep + env.get("PATH", "")

    for slim_file in slim_files:
        stem = os.path.splitext(os.path.basename(slim_file))[0]  # e.g. CHIP2EnvN

        print(f"\n{'='*60}")
        print(f"  Script: {stem}.slim")
        print(f"{'='*60}")

        # Extract AGE / REP from first comment line (e.g. //slim -d AGE=100 -d REP=1 ...)
        params = _extract_slim_params(slim_file)
        age = "100"; rep = "1"
        i = 0
        while i < len(params):
            if params[i] == "-d" and i + 1 < len(params):
                kv = params[i + 1]
                if kv.startswith("AGE="):
                    age = kv.split("=", 1)[1]
                elif kv.startswith("REP="):
                    rep = kv.split("=", 1)[1]
                i += 2
            else:
                i += 1

        nwk_path   = os.path.join(dir_newick, f"{stem}.nwk")
        trees_path = os.path.join(dir_tree,   f"{stem}.tree")

        if use_slim2tree:
            # ── Slim2tree3.py handles SLiM + newick in one call ──────────────
            # It calls: slim -d AGE=<age> -d REP=<rep> CHIP2EnvN.slim
            # and writes: CHIPage_{age}_rep{rep}.trees  +  .nwk  in cs_dir
            cmd = [sys.executable, slim2tree_py, age, rep]
            print(f"[INFO] Running: python Slim2tree3.py {age} {rep}")
            res = subprocess.run(cmd, cwd=cs_dir, capture_output=True,
                                 text=True, env=env)
            if res.stdout: print(res.stdout)
            if res.stderr: print(res.stderr)
            if res.returncode != 0:
                print(f"[ERROR] Slim2tree3.py failed (code {res.returncode})")
                continue

            # Expected output filenames (matching Slim2tree3.py conventions)
            src_trees = os.path.join(cs_dir, f"CHIPage_{age}_rep{rep}.trees")
            src_nwk   = os.path.join(cs_dir, f"CHIPage_{age}_rep{rep}.trees.nwk")

            if not os.path.exists(src_trees):
                print(f"[ERROR] Expected tree file not found: {src_trees}")
                continue
            shutil.move(src_trees, trees_path)
            print(f"[INFO] Tree → {trees_path}")

            if not os.path.exists(src_nwk):
                print(f"[ERROR] Expected newick file not found: {src_nwk}")
                continue
            shutil.move(src_nwk, nwk_path)
            print(f"[INFO] Newick → {nwk_path}")

        else:
            # ── Fallback: run SLiM directly + slim_newick for newick ─────────
            before = set(glob.glob(os.path.join(dir_tree, "*.trees")) +
                         glob.glob(os.path.join(dir_tree, "*.tree")))
            cmd = [slim_exe] + params + [os.path.abspath(slim_file)]
            print(f"[INFO] Running: {' '.join(cmd)}")
            res = subprocess.run(cmd, cwd=dir_tree, capture_output=True,
                                 text=True, env=env)
            if res.stdout: print(res.stdout)
            if res.stderr: print(res.stderr)
            if res.returncode != 0:
                print(f"[ERROR] SLiM failed for {stem} (code {res.returncode})")
                continue
            after     = set(glob.glob(os.path.join(dir_tree, "*.trees")) +
                            glob.glob(os.path.join(dir_tree, "*.tree")))
            new_files = after - before
            if not new_files:
                print(f"[ERROR] No tree file produced by {stem}.slim")
                continue
            produced = sorted(new_files)[0]
            if produced != trees_path:
                os.replace(produced, trees_path)
            print(f"[INFO] Tree → {trees_path}")
            if not convert_to_newick(trees_path, nwk_path):
                continue

        # ── PNG plots — all 4 variants into pngTree/ ─────────────────────────
        plot_png(nwk_path, dir_png, stem, stem, "CaseStudy")

        print(f"[DONE] {stem}")
        print(f"  tree/    {stem}.tree")
        print(f"  newick/  {stem}.nwk")
        print(f"  pngTree/ {stem}_horizontal_labels.png / _no_labels.png / _vertical_*.png")

    print(f"\n{'='*60}")
    print(f"  CaseStudy complete. Outputs in: {out_base}")
    print(f"{'='*60}")


# ── ReadyTrees: process pre-existing .trees files ────────────────────────────

def run_readytrees(script_dir):
    """
    Scans ReadyTrees/ (next to simtreemaker.py) for *.trees files.
    For each file, converts to Newick and generates 4 PNGs.
    Outputs go to ReadyTreesOutputs/{stem}Output/newick/ and pngTree/.
    """
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
        print(f"[ReadyTrees] Drop .trees files there and re-run.")
        return

    print(f"\n[ReadyTrees] Found {len(tree_files)} file(s) in {in_dir}")
    print(f"[ReadyTrees] Outputs → {out_dir}\n")

    for trees_path in tree_files:
        stem = os.path.splitext(os.path.basename(trees_path))[0]
        run_out  = os.path.join(out_dir, f"{stem}Output")
        dir_nwk  = os.path.join(run_out, "newick")
        dir_png  = os.path.join(run_out, "pngTree")
        for d in (dir_nwk, dir_png):
            os.makedirs(d, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"  File: {os.path.basename(trees_path)}")
        print(f"  Out : {run_out}")
        print(f"{'='*60}")

        nwk_path = os.path.join(dir_nwk, f"{stem}.nwk")
        if not convert_to_newick(trees_path, nwk_path):
            continue

        plot_png(nwk_path, dir_png, stem, stem, "ReadyTree")

        print(f"[DONE] {stem}")
        print(f"  newick/  {stem}.nwk")
        print(f"  pngTree/ {stem}_horizontal_labels.png / _no_labels.png")
        print(f"           {stem}_vertical_labels.png   / _no_labels.png")

    print(f"\n{'='*60}")
    print(f"  ReadyTrees complete. Outputs in: {out_dir}")
    print(f"{'='*60}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    script_dir  = os.path.dirname(os.path.abspath(__file__))
    config      = load_config(os.path.join(script_dir, "slim_config.txt"))

    slim_exe = config.get("SLIM_EXE")

    if not slim_exe:
        print("[ERROR] slim_config.txt must define SLIM_EXE.")
        sys.exit(1)

    # WORK_DIR is always the folder containing simtreemaker.py — makes it portable
    work_dir = script_dir
    os.makedirs(work_dir, exist_ok=True)

    # ── CaseStudy mode ────────────────────────────────────────────────────────
    if len(sys.argv) >= 2 and sys.argv[1] == "CaseStudy":
        run_casestudy(work_dir, slim_exe)
        return

    # ── ReadyTrees mode ───────────────────────────────────────────────────────
    if len(sys.argv) >= 2 and sys.argv[1] == "Tree":
        run_readytrees(script_dir)
        return

    # ── Options/ folder: named CSV files ─────────────────────────────────────
    options_dir = os.path.join(work_dir, "Options")

    if len(sys.argv) >= 2:
        # e.g.  python simtreemaker.py MutationSpread
        #        python simtreemaker.py MutationSpread.csv
        arg      = sys.argv[1]
        csv_name = arg if arg.endswith(".csv") else arg + ".csv"
        csv_path = os.path.join(options_dir, csv_name)
        if not os.path.exists(csv_path):
            print(f"[ERROR] CSV not found: {csv_path}")
            print(f"  Place your CSV files in: {options_dir}")
            sys.exit(1)
        csv_files = [csv_path]
        print(f"[INFO] Using: {csv_path}")
    else:
        print("Usage: python simtreemaker.py <CsvName|CaseStudy|Tree>")
        print("  Examples: MutationSpread, ClonalGrowth, Metastasis")
        sys.exit(1)

    for csv_path in csv_files:
        models = load_all_models(csv_path)
        print(f"\n[INFO] {os.path.basename(csv_path)}: {len(models)} row(s)")
        for model in models:
            run_model(model, slim_exe, work_dir, script_dir)

    print(f"\n{'='*60}")
    print(f"  All runs complete. Outputs in: {work_dir}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
