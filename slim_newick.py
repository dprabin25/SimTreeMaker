# -*- coding: utf-8 -*-
"""
slim_newick.py
Two sampling strategies for SLiM .trees files.

  Snapshot  — samples only cells alive at the final generation (like a biopsy).
              Parameters read from TreeOptions/Snapshot.txt.

  Lineage   — samples cells across ALL generations with recency bias.
              Recent generations get more cells, older ones get fewer.
              Parameters read from TreeOptions/Lineage.txt.

Install: conda install -c conda-forge pyslim=1.1.1 tskit=1.0.3
"""

import tskit
import random


# ── Shared helpers ────────────────────────────────────────────────────────────

def _path_length_to_root(tree, node):
    length = 0.0
    while node != tskit.NULL:
        parent = tree.parent(node)
        if parent == tskit.NULL:
            break
        length += tree.branch_length(node)
        node = parent
    return length


def _build_tip_labels(simplified_ts, tree, node_labels, node_to_population, mutation_position):
    """Build tip labels: gen{T}_ind{i}|pop{N}|{mutations}"""
    # Map each site position to m1 or m2
    mutation_type_map = {}
    for site in simplified_ts.sites():
        mut_type = "m2" if int(site.position) == int(mutation_position) else "m1"
        for mut in site.mutations:
            mutation_type_map[mut.id] = mut_type

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


# ── Snapshot ──────────────────────────────────────────────────────────────────

def build_newick_snapshot(
    trees_path,
    nwk_path,
    mutation_position  = 10000,
    valid_pops         = None,     # list of SLiM pop IDs, or None = all
    snapshot_samples   = 400,      # max cells drawn from final generation
    keep_unary         = False,    # F = compressed branches | T = real branch lengths
    seed               = None,     # random seed for reproducibility
):
    """
    Snapshot strategy.
    Only individuals alive at time=0 (final generation) are sampled.
    Up to snapshot_samples are drawn randomly.
    """
    if seed is not None:
        random.seed(seed)
        print(f"[newick] Random seed: {seed}")
    print(f"[newick] Loading: {trees_path}")
    ts = tskit.load(trees_path)
    print(f"[newick] {ts.num_individuals} individuals  |  "
          f"{ts.num_trees} tree(s)  |  {ts.num_sites} site(s)")

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
            "Check valid_pops in TreeOptions/Snapshot.txt."
        )

    print(f"[newick] Alive at final generation: {len(alive_inds)}")

    sampled = (
        random.sample(alive_inds, snapshot_samples)
        if len(alive_inds) > snapshot_samples
        else alive_inds
    )

    sampled_nodes      = []
    node_labels        = {}
    node_to_population = {}

    for i, ind in enumerate(sampled):
        node_id = ind.nodes[0]
        sampled_nodes.append(node_id)
        node_labels[node_id]        = f"gen0_ind{i}"
        node_to_population[node_id] = ind.population

    print(f"[newick] Sampled {len(sampled_nodes)} cells from final generation (Snapshot).")

    simplified_ts = ts.simplify(samples=sampled_nodes, keep_unary=keep_unary)
    tree          = simplified_ts.first()

    simplified_labels = _build_tip_labels(
        simplified_ts, tree, node_labels, node_to_population, mutation_position
    )
    _print_root_distances(tree, simplified_labels)

    newick, n_roots = _write_newick(tree, simplified_labels, nwk_path)
    print(f"\n[newick] Saved: {nwk_path}  ({n_roots} root(s))")
    return {"n_sampled": len(sampled_nodes), "n_roots": n_roots, "newick": newick}


# ── Lineage ───────────────────────────────────────────────────────────────────

def build_newick(
    trees_path,
    nwk_path,
    mutation_position          = 10000,
    valid_pops                 = None,  # list of SLiM pop IDs, or None = all
    max_samples_per_generation = 5,
    min_samples_per_generation = 1,
    keep_unary                 = True,  # T = real branch lengths | F = compressed
    seed                       = None,  # random seed for reproducibility
):
    """
    Lineage strategy.
    Individuals are sampled across ALL generations with recency bias:
      recency_weight = 1.0 - (gen / max_time)
      num_to_sample  = round(recency_weight * max_samples_per_generation)
      clamped to [min_samples_per_generation, cells available in that generation]
    """
    if seed is not None:
        random.seed(seed)
        print(f"[newick] Random seed: {seed}")
    print(f"[newick] Loading: {trees_path}")
    ts = tskit.load(trees_path)
    print(f"[newick] {ts.num_individuals} individuals  |  "
          f"{ts.num_trees} tree(s)  |  {ts.num_sites} site(s)")

    # Group individuals by generation
    gen_to_inds = {}
    for ind in ts.individuals():
        if valid_pops is not None and ind.population not in valid_pops:
            continue
        if len(ind.nodes) == 0:
            continue
        gen_to_inds.setdefault(ind.time, []).append(ind)

    if not gen_to_inds:
        raise ValueError(
            "No individuals found. Check valid_pops in TreeOptions/Lineage.txt."
        )

    sorted_gens = sorted(gen_to_inds.keys())
    max_time    = max(sorted_gens) if sorted_gens else 0
    print(f"[newick] Generations found: {len(sorted_gens)}  "
          f"(time {int(min(sorted_gens))} – {int(max_time)})")

    sampled_nodes      = []
    node_labels        = {}
    node_to_population = {}

    for gen in sorted_gens:
        inds = gen_to_inds[gen]
        if len(inds) < min_samples_per_generation:
            continue
        recency_weight = 1.0 - (gen / max_time) if max_time > 0 else 1.0
        num_to_sample  = round(recency_weight * max_samples_per_generation)
        num_to_sample  = max(min_samples_per_generation, min(num_to_sample, len(inds)))
        sampled = random.sample(inds, num_to_sample)
        for i, ind in enumerate(sampled):
            node_id = ind.nodes[0]
            sampled_nodes.append(node_id)
            node_labels[node_id]        = f"gen{int(gen)}_ind{i}"
            node_to_population[node_id] = ind.population

    print(f"[newick] Sampled {len(sampled_nodes)} nodes across "
          f"{len(sorted_gens)} generations (Lineage, recency bias).")

    if not sampled_nodes:
        raise ValueError(
            "No nodes sampled. Lower min_samples_per_generation in TreeOptions/Lineage.txt."
        )

    simplified_ts = ts.simplify(samples=sampled_nodes, keep_unary=keep_unary)
    tree          = simplified_ts.first()

    simplified_labels = _build_tip_labels(
        simplified_ts, tree, node_labels, node_to_population, mutation_position
    )
    _print_root_distances(tree, simplified_labels)

    newick, n_roots = _write_newick(tree, simplified_labels, nwk_path)
    print(f"\n[newick] Saved: {nwk_path}  ({n_roots} root(s))")
    return {"n_sampled": len(sampled_nodes), "n_roots": n_roots, "newick": newick}
