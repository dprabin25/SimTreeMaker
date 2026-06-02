# -*- coding: utf-8 -*-
"""
slim_newick.py
Standardized Newick conversion for ALL model types (regular + CaseStudy).

Replaces Slim2tree3.py entirely.

Approach (matching Slim2tree3.py exactly):
  - Find individuals alive at the FINAL generation using pyslim (most accurate)
  - Falls back to tskit built-in, then time==0 filter if pyslim not installed
  - Randomly sample up to max_samples (default 400)
  - Simplify tree sequence
  - Export Newick with tip labels (population + mutation info)

Install:  pip install tskit pyslim biopython
"""

import random


# ── Alive-individual detection (pyslim → tskit → fallback) ───────────────────

def _get_alive_individuals(ts, time=0):
    """
    Return Individual objects alive at `time`.
    Tries pyslim first (most accurate for nonWF), then tskit built-in,
    then filters by ind.time == 0 as last resort.
    """
    # 1. pyslim — handles multi-generation survival correctly
    try:
        import pyslim
        alive_ids  = pyslim.individuals_alive_at(ts, time)
        alive_inds = [ts.individual(i) for i in alive_ids
                      if len(ts.individual(i).nodes) > 0]
        print(f"[newick] pyslim: {len(alive_inds)} individuals alive at time={time}")
        return alive_inds
    except (ImportError, Exception) as e:
        if "pyslim" not in str(type(e).__module__):
            pass  # pyslim not installed — try next method

    # 2. tskit built-in (available in newer tskit for SLiM tree sequences)
    try:
        alive_ids  = ts.individuals_alive_at(time)
        alive_inds = [ts.individual(i) for i in alive_ids
                      if len(ts.individual(i).nodes) > 0]
        print(f"[newick] tskit.individuals_alive_at: {len(alive_inds)} alive")
        return alive_inds
    except AttributeError:
        pass

    # 3. Fallback — individuals born at time==0
    alive_inds = [ind for ind in ts.individuals()
                  if ind.time == time and len(ind.nodes) > 0]
    if not alive_inds:
        # Use most recent generation available
        min_time   = min((ind.time for ind in ts.individuals() if len(ind.nodes) > 0),
                         default=0)
        alive_inds = [ind for ind in ts.individuals()
                      if ind.time == min_time and len(ind.nodes) > 0]
        print(f"[newick] Fallback: most recent generation (time={min_time}), "
              f"{len(alive_inds)} individuals")
    else:
        print(f"[newick] Fallback (time=={time}): {len(alive_inds)} individuals")
    return alive_inds


# ── Main entry point ──────────────────────────────────────────────────────────

def build_newick(
    trees_path,
    nwk_path,
    mutation_position=10000,   # position of driver mutation on chromosome
    max_samples=400,           # max individuals to sample (matches Slim2tree3.py)
):
    """
    Load a SLiM .trees file, sample alive individuals at final generation,
    simplify the tree, label tips, and write a Newick file.
    Returns a summary dict.
    """
    try:
        import tskit
    except ImportError:
        raise ImportError("tskit not installed. Run: pip install tskit")

    print(f"[newick] Loading: {trees_path}")
    ts = tskit.load(trees_path)

    # Count tips across all trees (matching Slim2tree3.py diagnostic)
    n_tips = sum(tree.num_samples() for tree in ts.trees())
    print(f"[newick] Tree sequence: {n_tips} total sample nodes")

    # ── Get alive individuals ─────────────────────────────────────────────────
    alive_inds = _get_alive_individuals(ts, time=0)

    if not alive_inds:
        raise ValueError("No alive individuals found. Check the tree sequence.")

    # ── Sample up to max_samples (matching Slim2tree3.py) ────────────────────
    if len(alive_inds) > max_samples:
        alive_inds = random.sample(alive_inds, max_samples)
        print(f"[newick] Randomly sampled {max_samples} individuals")
    else:
        print(f"[newick] Using all {len(alive_inds)} individuals")

    # First node per individual (haploid), matching Slim2tree3.py
    sample_nodes = [ind.nodes[0] for ind in alive_inds]

    # ── Simplify ──────────────────────────────────────────────────────────────
    ts_simplified = ts.simplify(samples=sample_nodes)
    tree = ts_simplified.first()
    roots = list(tree.roots)
    print(f"[newick] Simplified: {len(sample_nodes)} samples | "
          f"{ts_simplified.num_trees} tree(s) | {len(roots)} root(s)")

    # ── Mutation-type map ─────────────────────────────────────────────────────
    mutation_type_map = {}
    for site in ts_simplified.sites():
        mut_type = "m2" if int(site.position) == int(mutation_position) else "m1"
        for mut in site.mutations:
            mutation_type_map[mut.id] = mut_type

    mutations_by_node = {n: [] for n in tree.nodes()}
    for mut in ts_simplified.mutations():
        mutations_by_node[mut.node].append(
            mutation_type_map.get(mut.id, "unknown")
        )

    # ── Tip labels: ind{i} | pop{N} | mutations ───────────────────────────────
    node_to_pop = {}
    for ind in alive_inds:
        for n in ind.nodes:
            node_to_pop[n] = ind.population

    tip_labels = {}
    for i, node in enumerate(tree.leaves()):
        inherited = set()
        cur = node
        while cur != tskit.NULL:
            inherited.update(mutations_by_node.get(cur, []))
            cur = tree.parent(cur)
        pop  = node_to_pop.get(node, -1)
        muts = "+".join(sorted(inherited)) if inherited else "no_mut"
        tip_labels[node] = f"ind{i}|pop{pop}|{muts}"

    # ── Root distances (diagnostic, matching Slim2tree3.py) ──────────────────
    print("\n[newick] Root distances (sampled tips):")
    for node in tree.leaves():
        print(f"  {tip_labels.get(node, f'node{node}')} -> "
              f"{_path_length(tree, node):.2f}")

    # ── Export Newick ─────────────────────────────────────────────────────────
    if len(roots) == 1:
        newick = tree.as_newick(
            root=roots[0],
            node_labels=tip_labels,
            include_branch_lengths=True,
        )
    else:
        parts = [
            tree.as_newick(
                root=r,
                node_labels=tip_labels,
                include_branch_lengths=True,
            ).rstrip(";")
            for r in roots
        ]
        newick = "(" + ",".join(parts) + ");"

    with open(nwk_path, "w") as f:
        f.write(newick)

    print(f"\n[newick] Saved: {nwk_path}")
    return {
        "n_sampled": len(sample_nodes),
        "n_roots":   len(roots),
        "newick":    newick,
    }


def _path_length(tree, node):
    import tskit
    length = 0.0
    while node != tskit.NULL:
        parent = tree.parent(node)
        if parent == tskit.NULL:
            break
        length += tree.branch_length(node)
        node = parent
    return length
