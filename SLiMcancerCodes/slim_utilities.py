"""
slim_utilities.py — shared plumbing used by both simtreemaker.py (TumorGrowth_NoMigration)
and casestudy.py (arbitrary CaseStudy/ scripts):

  - load_config()          simple key = value text file parser
  - load_tree_options()    fixed Snapshot sampling defaults
  - run_slim()             invoke SLiM as a subprocess, streaming output live
  - convert_to_newick()    .trees -> .nwk (+ optional SNV CSV), Snapshot sampling
  - mutation_type_label() / collect_snv_rows() / write_snv_csv()
                           per-cell SNV profiling, called from convert_to_newick()
  - plot_png_r() / plot_png()   tree PNG rendering (R+ape via plot_tree.R, or
                                 matplotlib fallback)

Kept in its own module (rather than living inside simtreemaker.py) so both
runner scripts can share it without either one depending on the other.
Everything that used to be split across sim_utils.py / slim_newick.py /
snvprofiling.py now lives here in one file — the Newick-building and SNV
steps still run within the same sampling pass (inside convert_to_newick),
so the .nwk and the SNV CSV always describe the exact same sampled cells.
"""

import os
import sys
import csv
import random
import shutil
import subprocess

import numpy as np
import tskit

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ── Config file loading ─────────────────────────────────────────────────────

def load_config(config_path):
    """Read a simple `key = value` text file into a dict. Inline '#' comments
    (after the value) are stripped, so this works for both slim_config.txt
    and config.txt."""
    if not os.path.exists(config_path):
        print(f"[ERROR] Config file not found: {config_path}")
        sys.exit(1)
    config = {}
    with open(config_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                key, _, value = line.partition("=")
                value = value.split("#")[0]  # strip inline comments
                config[key.strip()] = value.strip()
    return config


def parse_valid_pops(vp_raw):
    """Parse valid_pops string → list of ints, or None for 'all'."""
    if vp_raw.strip().lower() == "all":
        return None
    return [int(x.strip()) for x in vp_raw.split(",") if x.strip()]


def load_tree_options(script_dir):
    """
    Sampling options for the Newick/SNV snapshot. Every method in this
    project samples the same way — a snapshot at the simulation's final
    generation — so these are fixed defaults rather than something read
    from a config file (there used to be an optional TreeOptions/Snapshot.txt
    for this; it's no longer used). To change sampling behavior, edit the
    values below directly.
    """
    return {
        "mode": "Snapshot",
        "valid_pops": "all",
        "snapshotSamples": "400",
        "seed": "none",
    }


# ── Run SLiM ─────────────────────────────────────────────────────────────────

def run_slim(slim_exe, script_path, extra_args=None):
    """Run SLiM, streaming output live to terminal and capturing it for the log.
    extra_args: list of additional flags e.g. ['-d', 'AGE=100', '-d', 'REP=1']
    """
    cmd = [slim_exe] + (extra_args or []) + [script_path]
    print(f"[INFO] Running SLiM: {' '.join(cmd)}")
    process = subprocess.Popen(
        cmd,
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


# ── Mutation typing + SNV profiling ─────────────────────────────────────────
# (formerly snvprofiling.py)

def mutation_type_label(mut, site_position, m2_range):
    """
    Determine whether a tskit Mutation is SLiM type m1 or m2.

    Prefers the mutation's own SLiM metadata (mutation_type id: 1 = m1,
    2 = m2). This is required under the multi-driver design, where m2 can
    arise spontaneously at any position in its genomic element's range.
    Falls back to checking whether the position lies inside m2_range
    (an inclusive (start, end) tuple, from config.txt's m2Start/m2End)
    only if SLiM metadata isn't present (e.g. a .trees file written
    without it).
    """
    try:
        meta = mut.metadata
        if meta and meta.get("mutation_list"):
            type_id = meta["mutation_list"][0].get("mutation_type")
            if type_id in (1, 2, 3):                  # m3 = resistance (TG_Treatment.slim)
                return f"m{type_id}"
    except Exception:
        pass
    if m2_range is not None and m2_range[0] <= int(site_position) <= m2_range[1]:
        return "m2"
    return "m1"


def collect_snv_rows(simplified_ts, tree, node_labels, node_to_population,
                      m2_range, sim_end_tick=None):
    """
    For each sampled leaf, return one row per mutation it carries (its own
    plus everything inherited from its ancestors): cell_id, population,
    mutation_type, genomic position, and how long ago the mutation arose.

    If sim_end_tick is given, an additional origin_tick column is filled in
    as the forward-time generation the mutation arose (sim_end_tick minus
    the tree-sequence's "time ago" for that mutation). Without it, only the
    raw time-ago value is reported.
    """
    mutation_info = {}  # mut.id -> (type, position, time_ago)
    for site in simplified_ts.sites():
        for mut in site.mutations:
            mut_type = mutation_type_label(mut, site.position, m2_range)
            mutation_info[mut.id] = (mut_type, site.position, mut.time)

    mutations_by_node = {n: [] for n in tree.nodes()}
    for mut in simplified_ts.mutations():
        mutations_by_node[mut.node].append(mutation_info.get(mut.id, ("unknown", -1, -1.0)))

    rows = []
    for node in tree.leaves():
        cell_id = node_labels.get(node, f"node{node}")
        pop     = node_to_population.get(node, -1)
        seen    = set()
        carried = []
        cur = node
        while cur != tskit.NULL:
            for mut_type, position, time_ago in mutations_by_node.get(cur, []):
                key = (mut_type, position)
                if key in seen:
                    continue
                seen.add(key)
                carried.append((mut_type, position, time_ago))
            cur = tree.parent(cur)

        if not carried:
            rows.append({
                "cell_id": cell_id, "population": pop, "mutation_type": "none",
                "position": "", "generations_before_sampling": "", "origin_tick": "",
            })
            continue

        for mut_type, position, time_ago in carried:
            origin_tick = round(sim_end_tick - time_ago) if sim_end_tick is not None else ""
            rows.append({
                "cell_id": cell_id,
                "population": pop,
                "mutation_type": mut_type,
                "position": int(position),
                "generations_before_sampling": round(time_ago, 3),
                "origin_tick": origin_tick,
            })
    return rows


# ── FASTA output (baseline, no dropout) ─────────────────────────────────────
#
#   parent genome  : the founder cell's starting genome (ts.sequence_length bp)
#   full FASTA     : each sampled cell = parent genome + every SNV it carries
#   SNV FASTA      : the same cells, only the SNV columns (exactly the
#                    positions in the SNV CSV, ascending), + _positions.tsv
#
# Headers are the Newick tip labels (gen0_ind<i>|pop<N>|<muts>): the first
# "|" field is the SNV CSV's cell_id, and the full header matches the .nwk.
#
# TG.slim is not nucleotide-based: the .trees file records where/when each
# mutation arose and which cells carry it, but not the base it changed to.
# The parent genome is therefore drawn at random (uniform A/C/G/T) and each
# SNV's ALT base is a random base != the parent base at that position. Both
# use fastaSeed, so the same .trees + seed always gives identical FASTA.

_BASES = "ACGT"


def _write_records(path, headers, rows_u8, line_width=60):
    import numpy as np
    lut = np.frombuffer(_BASES.encode(), dtype=np.uint8)
    with open(path, "w", newline="\n") as fh:
        for h, codes in zip(headers, rows_u8):
            seq = lut[codes].tobytes().decode()
            fh.write(f">{h}\n")
            for i in range(0, len(seq), line_width):
                fh.write(seq[i:i + line_width] + "\n")


def write_fasta_set(simplified_ts, tip_labels, out_prefix, m2_range,
                    fasta_seed=12345, write_full=True):
    """
    Writes, using out_prefix (e.g. .../fasta/TumorGrowth):
        <prefix>_parent.fasta         starting genome
        <prefix>_full.fasta           full genome per sampled cell (if write_full)
        <prefix>_snv.fasta            SNV columns only, per sampled cell
        <prefix>_snv_positions.tsv    column -> position, m1/m2, REF, ALT
    Returns dict with the in-memory SNV matrix for downstream steps.
    """
    import numpy as np

    rng     = np.random.default_rng(fasta_seed)
    seq_len = int(simplified_ts.sequence_length)
    parent  = rng.integers(0, 4, size=seq_len, dtype=np.uint8)
    samples = list(simplified_ts.samples())
    headers = [tip_labels.get(n, f"node{n}") for n in samples]

    cols, info = [], []
    for var in simplified_ts.variants(isolated_as_missing=False):
        if not var.site.mutations:
            continue
        pos      = int(var.site.position)
        ref_base = int(parent[pos])
        alt_pool = [b for b in range(4) if b != ref_base]
        rng.shuffle(alt_pool)
        # Allele 0 = SLiM ancestral state (= parent base); each derived allele
        # (stacked mutations at one site) gets its own non-REF base.
        allele_base = [ref_base] + [alt_pool[i % 3] for i in range(len(var.alleles) - 1)]
        cols.append(np.array(allele_base, dtype=np.uint8)[var.genotypes])
        mtype = mutation_type_label(var.site.mutations[0], var.site.position, m2_range)
        info.append((pos, mtype, _BASES[ref_base],
                     ",".join(_BASES[b] for b in allele_base[1:])))

    snv = (np.stack(cols, axis=1) if cols
           else np.zeros((len(samples), 0), dtype=np.uint8))
    positions = np.array([i[0] for i in info], dtype=np.int64)

    parent_path = out_prefix + "_parent.fasta"
    _write_records(parent_path, [f"parent length={seq_len} fastaSeed={fasta_seed}"], [parent])
    print(f"[fasta] Saved: {parent_path}  ({seq_len} bp)")

    full_path = None
    if write_full:
        full_path = out_prefix + "_full.fasta"
        lut = np.frombuffer(_BASES.encode(), dtype=np.uint8)
        with open(full_path, "w", newline="\n") as fh:
            for row, h in enumerate(headers):
                genome = parent.copy()
                genome[positions] = snv[row]          # parent + this cell's SNVs
                seq = lut[genome].tobytes().decode()
                fh.write(f">{h}\n")
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i + 60] + "\n")
        print(f"[fasta] Saved: {full_path}  ({len(samples)} cells x {seq_len} bp)")

    snv_path = out_prefix + "_snv.fasta"
    _write_records(snv_path, headers, snv)
    pos_path = out_prefix + "_snv_positions.tsv"
    with open(pos_path, "w", newline="\n") as fh:
        fh.write("column\tposition\tmutation_type\tref\talt\n")
        for i, (pos, mtype, ref, alt) in enumerate(info, 1):
            fh.write(f"{i}\t{pos}\t{mtype}\t{ref}\t{alt}\n")
    print(f"[fasta] Saved: {snv_path}  ({len(samples)} cells x {len(info)} SNV column(s))")
    print(f"[fasta] Saved: {pos_path}")

    return {"headers": headers, "snv": snv, "info": info,
            "parent": parent_path, "full": full_path, "snv_fasta": snv_path,
            "positions": pos_path}


def write_snv_csv(rows, csv_path):
    fieldnames = ["cell_id", "population", "mutation_type", "position",
                  "generations_before_sampling", "origin_tick"]
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"[snv] Saved: {csv_path}  ({len(rows)} row(s))")


# ── Newick building (Snapshot sampling strategy) ────────────────────────────
# (formerly slim_newick.py)
#
# Samples only cells alive at the final generation (like a biopsy).
# Sampling defaults come from load_tree_options().
# Install: conda install -c conda-forge pyslim=1.1.1 tskit=1.0.3

def _path_length_to_root(tree, node):
    length = 0.0
    while node != tskit.NULL:
        parent = tree.parent(node)
        if parent == tskit.NULL:
            break
        length += tree.branch_length(node)
        node = parent
    return length


def _build_tip_labels(simplified_ts, tree, node_labels, node_to_population, m2_range):
    """Build tip labels: gen{T}_ind{i}|pop{N}|{mutations}"""
    # Map each mutation to m1 or m2 (see mutation_type_label above)
    mutation_type_map = {}
    for site in simplified_ts.sites():
        for mut in site.mutations:
            mutation_type_map[mut.id] = mutation_type_label(mut, site.position, m2_range)

    # Collect mutations per node
    mutations_by_node = {n: [] for n in tree.nodes()}
    for mut in simplified_ts.mutations():
        mutations_by_node[mut.node].append(mutation_type_map.get(mut.id, "unknown"))

    # Build full tip label for each leaf
    simplified_labels = {}
    for node in tree.leaves():
        inherited = set()
        cur = node
        while cur != tskit.NULL:
            inherited.update(mutations_by_node.get(cur, []))
            cur = tree.parent(cur)
        base = node_labels.get(node, f"node{node}")
        pop  = node_to_population.get(node, -1)
        muts = "+".join(sorted(inherited)) if inherited else "no_mut"
        simplified_labels[node] = f"{base}|pop{pop}|{muts}"

    return simplified_labels


def _print_root_distances(tree, simplified_labels):
    print("\n[newick] Root distances for sampled nodes:")
    for node in tree.leaves():
        label = simplified_labels.get(node, f"node{node}")
        dist  = _path_length_to_root(tree, node)
        print(f"  {label} → root path length: {dist:.2f}")


def _write_newick(tree, simplified_labels, nwk_path):
    """Join roots into one Newick string and write to file."""
    roots = list(tree.roots)
    if len(roots) == 1:
        newick = tree.as_newick(
            root        = roots[0],
            node_labels = simplified_labels,
            include_branch_lengths = True,
        )
    else:
        subtrees = [
            tree.as_newick(
                root        = r,
                node_labels = simplified_labels,
                include_branch_lengths = True,
            ).rstrip(";")
            for r in roots
        ]
        newick = "(" + ",".join(subtrees) + ");"

    with open(nwk_path, "w") as f:
        f.write(newick)
    return newick, len(roots)


def _merge_genome_copies(ts, sampled, keep_unary=False):
    """
    Return (merged_ts, sample node ids) where every mutation carried on a
    sampled cell's second genome copy is moved onto the matching node of its
    first copy, so tracing copy 0 alone gives each cell's full mutation set.

    Reproduction is clonal with no recombination and copies a cell's two
    genome copies to its clone in order, so the copy-0 and copy-1 genealogies
    of the sampled cells are identical: the copy-0 counterpart of a copy-1
    node is the copy-0 MRCA of the same cells, at the same time. Only reading
    copy 0 (as before) silently dropped about half of all mutations.
    """
    k  = len(sampled)
    n0 = [ind.nodes[0] for ind in sampled]
    n1 = [ind.nodes[1] for ind in sampled]
    both = ts.simplify(samples=n0 + n1, keep_unary=keep_unary)   # samples: 0..k-1 copy 0, k..2k-1 copy 1
    tree = both.first()
    tables = both.dump_tables()
    nodes  = tables.mutations.node.copy()
    counterpart = {}
    for j, x in enumerate(nodes):
        x = int(x)
        if x not in counterpart:
            cells = [s - k for s in tree.samples(x)]
            if cells and min(cells) >= 0:                         # a copy-1 lineage
                counterpart[x] = cells[0] if len(cells) == 1 else tree.mrca(*cells)
            else:
                counterpart[x] = x
        nodes[j] = counterpart[x]
    tables.mutations.node   = nodes
    tables.mutations.parent = np.full(len(nodes), -1, dtype=np.int32)
    tables.sort()
    tables.build_index()
    tables.compute_mutation_parents()
    return tables.tree_sequence(), list(range(k))


def build_newick_snapshot(
    trees_path,
    nwk_path,
    m2_range           = None,
    valid_pops         = None,     # list of SLiM pop IDs, or None = all
    snapshot_samples   = 400,      # max cells drawn from final generation (pooled)
    sample_cells       = None,     # dict {pop_id: count} — overrides snapshot_samples
                                    # with independent per-population sample sizes
    keep_unary         = False,    # F = compressed branches | T = real branch lengths
    seed               = None,     # random seed for reproducibility
    snv_csv_path       = None,     # if given, also write a per-cell SNV profile CSV
    sim_end_tick       = None,     # generation the run stopped at (for origin_tick)
    fasta_prefix       = None,     # if given, write parent/full/SNV FASTA as <prefix>_*.fasta
    fasta_seed         = 12345,    # seed for the parent genome + ALT bases
    write_full_fasta   = False,    # also write the full-length per-cell FASTA
):
    """
    Snapshot strategy.
    Only individuals alive at time=0 (final generation) are sampled.
    Either up to snapshot_samples are drawn randomly across all valid_pops
    (pooled), or, if sample_cells is given, that many are drawn independently
    from each population.

    If snv_csv_path is given, the per-cell SNV profile CSV is built from this
    same sampled/simplified tree sequence — not a separate re-sampling pass —
    so it always matches exactly the cells written to nwk_path.
    """
    if seed is not None:
        random.seed(seed)
        print(f"[newick] Random seed: {seed}")
    print(f"[newick] Loading: {trees_path}")
    ts = tskit.load(trees_path)
    print(f"[newick] {ts.num_individuals} individuals  |  "
          f"{ts.num_trees} tree(s)  |  {ts.num_sites} site(s)")

    # If the caller didn't supply sim_end_tick (e.g. CaseStudy scripts, which
    # have no config.txt simulationEndTick to read), fall back to the actual
    # final tick SLiM recorded in the tree sequence's own metadata — this is
    # the true forward-time tick treeSeqOutput() was called at, so origin_tick
    # in the SNV CSV can still be computed instead of coming out blank.
    if sim_end_tick is None:
        try:
            sim_end_tick = ts.metadata["SLiM"]["tick"]
        except Exception:
            pass

    # Filter to alive individuals at final gen
    alive_inds = [
        ind for ind in ts.individuals()
        if ind.time == 0
        and len(ind.nodes) > 0
        and (valid_pops is None or ind.population in valid_pops)
    ]

    if not alive_inds:
        raise ValueError(
            "No individuals alive at final generation. "
            "Check valid_pops in load_tree_options()."
        )

    print(f"[newick] Alive at final generation: {len(alive_inds)}")

    if sample_cells:
        sampled = []
        by_pop = {}
        for ind in alive_inds:
            by_pop.setdefault(ind.population, []).append(ind)
        for pop_id, count in sample_cells.items():
            pool = by_pop.get(pop_id, [])
            sampled.extend(random.sample(pool, count) if len(pool) > count else pool)
    else:
        sampled = (
            random.sample(alive_inds, snapshot_samples)
            if len(alive_inds) > snapshot_samples
            else alive_inds
        )

    print(f"[newick] Sampled {len(sampled)} cells from final generation (Snapshot).")

    # Cells are diploid and a new mutation hits one of the two genome copies,
    # so a cell's mutations are those on either copy. Merge them onto copy 0
    # (see _merge_genome_copies); its tree is the cell tree.
    merged_ts, sampled_nodes = _merge_genome_copies(ts, sampled, keep_unary)
    node_labels        = {n: f"gen0_ind{i}" for i, n in enumerate(sampled_nodes)}
    node_to_population = {n: ind.population for n, ind in zip(sampled_nodes, sampled)}

    simplified_ts, node_map = merged_ts.simplify(samples=sampled_nodes, keep_unary=keep_unary,
                                                  map_nodes=True)
    tree = simplified_ts.first()

    # simplify() reassigns node IDs — translate the label/population dicts
    # (keyed by pre-simplification IDs) to the new IDs before looking them up.
    node_labels        = {int(node_map[old]): lbl  for old, lbl  in node_labels.items()}
    node_to_population = {int(node_map[old]): pop  for old, pop  in node_to_population.items()}

    simplified_labels = _build_tip_labels(
        simplified_ts, tree, node_labels, node_to_population, m2_range
    )
    _print_root_distances(tree, simplified_labels)

    newick, n_roots = _write_newick(tree, simplified_labels, nwk_path)
    print(f"\n[newick] Saved: {nwk_path}  ({n_roots} root(s))")

    if snv_csv_path:
        rows = collect_snv_rows(simplified_ts, tree, node_labels, node_to_population,
                                 m2_range, sim_end_tick=sim_end_tick)
        write_snv_csv(rows, snv_csv_path)

    if fasta_prefix:
        # Same simplified_ts and tip labels as the .nwk and SNV CSV, so the
        # FASTA describes exactly the same cells and SNVs.
        write_fasta_set(simplified_ts, simplified_labels, fasta_prefix,
                        m2_range, fasta_seed=fasta_seed,
                        write_full=write_full_fasta)

    return {"n_sampled": len(sampled_nodes), "n_roots": n_roots, "newick": newick}


# ── Convert .trees to .nwk (+ optional SNV CSV) ─────────────────────────────

# Set to the exception text whenever convert_to_newick() returns False, so
# callers that only check the boolean can still recover the real reason
# afterward for logging — e.g. casestudy.py's log.txt entries.
_last_newick_error = None


def convert_to_newick(trees_path, nwk_path, m2_range=None, tree_opts=None,
                       sample_cells=None, snv_csv_path=None, sim_end_tick=None,
                       fasta_prefix=None, fasta_seed=12345, write_full_fasta=False,
                       mutation_position=None):
    """
    mutation_position : older name, still accepted (CaseStudy/casestudy.py
                   and Trees/process_trees.py pass it): a single driver
                   position, used as m2_range = (pos, pos) when m2_range
                   isn't given. Only matters for .trees files without SLiM
                   mutation metadata.

    Build a Newick tree using the Snapshot strategy (samples cells alive at
    the final generation). tree_opts comes from load_tree_options().

    sample_cells : optional dict {population_id: count} for independent
                   per-population sampling — overrides the pooled
                   snapshotSamples cap when given.
    snv_csv_path : optional path — if given, also write a per-cell SNV
                   profile CSV alongside the Newick file.
    sim_end_tick : generation the run stopped at, used to convert mutation
                   ages into forward-time origin ticks in the SNV CSV.
    """
    global _last_newick_error
    if tree_opts is None:
        tree_opts = {"mode": "Snapshot"}
    if m2_range is None and mutation_position is not None:
        m2_range = (int(mutation_position), int(mutation_position))

    try:
        seed_raw = tree_opts.get("seed", "none").strip().lower()
        seed     = None if seed_raw == "none" else int(seed_raw)

        valid_pops       = parse_valid_pops(tree_opts.get("valid_pops", "all"))
        snapshot_samples = int(tree_opts.get("snapshotSamples", 400))
        stats = build_newick_snapshot(
            trees_path        = trees_path,
            nwk_path          = nwk_path,
            m2_range = m2_range,
            valid_pops        = valid_pops,
            snapshot_samples  = snapshot_samples,
            sample_cells      = sample_cells,
            keep_unary        = False,
            seed              = seed,
            snv_csv_path      = snv_csv_path,
            sim_end_tick      = sim_end_tick,
            fasta_prefix      = fasta_prefix,
            fasta_seed        = fasta_seed,
            write_full_fasta  = write_full_fasta,
        )

        print(f"[INFO] {stats['n_sampled']} nodes sampled | {stats['n_roots']} root(s)")
        _last_newick_error = None
        return True

    except Exception as e:
        import traceback
        _last_newick_error = traceback.format_exc()
        print(f"[ERROR] Newick conversion failed: {e}")
        return False


# ── PNG plots ────────────────────────────────────────────────────────────────

def plot_trajectory_r(traj_csv, out_png, dose_tick, script_dir):
    """Draw a treatment run's cells-per-generation trajectory (plot_trajectory.R).
    Returns True if the PNG was written."""
    rscript = find_rscript()
    r_script = os.path.join(script_dir, "plot_trajectory.R")
    if not rscript or not os.path.exists(r_script) or not os.path.exists(traj_csv):
        print("[WARN] Trajectory plot skipped (Rscript, plot_trajectory.R or the CSV not found).")
        return False
    result = subprocess.run([rscript, r_script, traj_csv, out_png, str(dose_tick)],
                            capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stderr.rstrip())
        print("[WARN] Trajectory plot failed.")
        return False
    print(f"[trajectory] Saved: {out_png}")
    return True


def find_rscript():
    """Return path to Rscript if available, else None."""
    return shutil.which("Rscript")


def plot_png_r(nwk_path, dir_png, stem, model_name, model_type, script_dir, mode="Snapshot", variants=None):
    """
    Render tree PNGs using R (ape). Returns True if successful.

    variants : optional list of keys to restrict which PNGs get rendered —
               any of "horizontal_labels", "horizontal_no_labels",
               "vertical_labels", "vertical_no_labels". None (default)
               renders all four.
    """
    rscript = find_rscript()
    if not rscript:
        return False
    r_script = os.path.join(script_dir, "plot_tree.R")
    if not os.path.exists(r_script):
        print("[WARN] plot_tree.R not found alongside slim_utilities.py")
        return False
    cmd = [rscript, r_script, nwk_path, dir_png, stem, model_name, model_type, mode]
    if variants is not None:
        cmd.append(",".join(variants))
    result = subprocess.run(cmd, capture_output=True, text=True)
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


def plot_png(nwk_path, dir_png, stem, model_name, model_type, mode="Snapshot", variants=None):
    """
    variants : optional list of keys to restrict which PNGs get rendered —
               any of "horizontal_labels", "horizontal_no_labels",
               "vertical_labels", "vertical_no_labels". None (default)
               renders all four.
    """
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

    all_variants = [
        ("horizontal", True,  (13, h), f"{stem}_horizontal_labels.png",    "horizontal_labels"),
        ("horizontal", False, (10, h), f"{stem}_horizontal_no_labels.png", "horizontal_no_labels"),
        ("vertical",   True,  (w, 10), f"{stem}_vertical_labels.png",     "vertical_labels"),
        ("vertical",   False, (w, 8),  f"{stem}_vertical_no_labels.png",  "vertical_no_labels"),
    ]
    selected = (
        all_variants if variants is None
        else [v for v in all_variants if v[4] in variants]
    )

    for orientation, show_labels, figsize, fname, _key in selected:
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
