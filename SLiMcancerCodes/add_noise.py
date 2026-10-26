"""
add_noise.py — adds single-cell sequencing noise to a clean SNV FASTA.

Runs automatically as the last pipeline step with
`python simtreemaker.py [migration | nomigration] noise` (settings from
config.txt's Sequencing noise section), or on its own on any aligned FASTA:

    python add_noise.py <fasta> <missing> <FP> <FN> <normal name>
                        [--seed N] [--cell-spread S] [--site-spread S]

    python add_noise.py Test.fasta 0.5 0.5 0.5 Nor
    python add_noise.py TG/NoMigration/fasta/TumorGrowth_snv.fasta 0.3 0.01 0.2 Normal --cell-spread 0.5

Noise model (standard single-cell DNA error model), for every site of every
cell except the normal one:

    missing  with that cell/site's missing rate          -> '?'
    otherwise, if the cell is mutated at that site (base != normal):
        false negative with probability FN -> normal base (allelic dropout)
    otherwise (cell has the normal base):
        false positive with probability FP -> the site's ALT base
            (the mutant base seen in other cells there; a random non-normal
            base only if no cell carries one)

FP and FN apply to the sites that are not missing. The normal sequence is
written unchanged, last.

Missing rate: <missing> is the average. With --cell-spread / --site-spread
above 0, each cell / site gets its own rate drawn from a Beta distribution
with that mean; the spread is the SD as a fraction of its maximum possible
value (0 = every cell/site identical, close to 1 = mostly all-or-nothing).
With both, a cell/site rate is cell_rate x site_rate / missing (capped at 1).

If <normal name> is not in the FASTA and the file is a simtreemaker
<name>_snv.fasta, the normal sequence is built from the REF column of the
neighbouring <name>_snv_positions.tsv (the parent genome's base at each SNV
site) and added under that name.

Outputs, in a noise/ folder next to the input (or --out-dir):
    <stem>_Miss_<m>_FP_<fp>_FN_<fn>.fasta        noisy cells + normal
    <stem>_Miss_<m>_FP_<fp>_FN_<fn>_log.tsv      every changed cell/site
    <stem>_Miss_<m>_FP_<fp>_FN_<fn>_cells.tsv    per-cell missing rate (target, observed)
"""

import os
import sys
import argparse

import numpy as np

BASES = np.array(list("ACGT"))


def read_fasta(path):
    names, seqs = [], {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:]
                names.append(name)
                seqs[name] = []
            else:
                seqs[name].append(line)
    return names, {n: "".join(s).upper() for n, s in seqs.items()}


def normal_from_positions(fasta_path):
    """REF bases from <name>_snv_positions.tsv next to a <name>_snv.fasta."""
    if not fasta_path.endswith("_snv.fasta"):
        return None
    tsv = fasta_path[:-len(".fasta")] + "_positions.tsv"
    if not os.path.exists(tsv):
        return None
    with open(tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ref_col = header.index("ref")
        return "".join(line.rstrip("\n").split("\t")[ref_col] for line in fh if line.strip())


def check_rate(value, name, upper=1.0, inclusive=True):
    value = float(value)
    ok = 0 <= value <= upper if inclusive else 0 <= value < upper
    if not ok:
        bound = f"<= {upper}" if inclusive else f"< {upper}"
        raise ValueError(f"{name} must be >= 0 and {bound} (got {value}).")
    return value


def beta_rates(rng, mean, spread, n):
    """n rates with the given mean; spread = SD / sqrt(mean*(1-mean))."""
    if spread == 0 or mean in (0.0, 1.0):
        return np.full(n, mean)
    k = 1.0 / spread ** 2 - 1.0                   # Beta concentration
    return rng.beta(mean * k, (1.0 - mean) * k, n)


def _fmt(x):
    return f"{x:g}"


def add_noise(fasta, missing, fp, fn, normal_name, seed=None, cell_spread=0.0, site_spread=0.0,
              out_dir=None):
    """Write the noisy FASTA + logs to out_dir (default: a noise/ folder next to
    `fasta`). Returns (output paths, counts dict)."""
    missing = check_rate(missing, "missing rate")
    fp      = check_rate(fp, "FP rate")
    fn      = check_rate(fn, "FN rate")
    cell_spread = check_rate(cell_spread, "cell spread", inclusive=False)
    site_spread = check_rate(site_spread, "site spread", inclusive=False)

    names, seqs = read_fasta(fasta)
    if normal_name in seqs:
        normal = seqs[normal_name]
    else:
        normal = normal_from_positions(fasta)
        if normal is None:
            raise ValueError(f"Normal sequence {normal_name!r} not found in {fasta}.")
        print(f"[noise] {normal_name!r} not in FASTA; built it from the REF column of _snv_positions.tsv.")
    cells = [n for n in names if n != normal_name]
    if not cells:
        raise ValueError("No cells to add noise to.")
    lengths = {len(seqs[c]) for c in cells} | {len(normal)}
    if len(lengths) != 1:
        raise ValueError(f"Sequences must all be the same length (found {sorted(lengths)}).")

    rng  = np.random.default_rng(seed)
    mat  = np.array([list(seqs[c]) for c in cells]).reshape(len(cells), -1)   # cells x sites
    norm = np.array(list(normal))
    n_cells, n_sites = mat.shape

    # ALT base per site: the most common non-normal base among the cells.
    alt = np.empty(n_sites, dtype="<U1")
    for j in range(n_sites):
        col = mat[:, j]
        mutant = col[(col != norm[j]) & np.isin(col, BASES)]
        if mutant.size:
            vals, counts = np.unique(mutant, return_counts=True)
            alt[j] = vals[np.argmax(counts)]
        else:
            alt[j] = rng.choice(BASES[BASES != norm[j]])

    # Per-cell x per-site missing probability, mean = `missing`.
    cell_rate = beta_rates(rng, missing, cell_spread, n_cells)
    site_rate = beta_rates(rng, missing, site_spread, n_sites)
    if cell_spread and site_spread:
        p_miss = np.minimum(1.0, np.outer(cell_rate, site_rate) / missing)
    elif site_spread:
        p_miss = np.broadcast_to(site_rate, mat.shape)
    else:
        p_miss = np.broadcast_to(cell_rate[:, None], mat.shape)

    mutated = mat != norm
    is_miss = rng.random(mat.shape) < p_miss
    flip    = rng.random(mat.shape)
    is_fn   = ~is_miss & mutated & (flip < fn)
    is_fp   = ~is_miss & ~mutated & (flip < fp)

    noisy = mat.copy()
    noisy[is_fn]   = np.broadcast_to(norm, mat.shape)[is_fn]
    noisy[is_fp]   = np.broadcast_to(alt, mat.shape)[is_fp]
    noisy[is_miss] = "?"

    out_dir = out_dir or os.path.join(os.path.dirname(os.path.abspath(fasta)), "noise")
    os.makedirs(out_dir, exist_ok=True)
    base = os.path.basename(fasta)
    stem = os.path.join(out_dir, base[:-len(".fasta")] if base.endswith(".fasta") else base)
    tag  ="_".join(["Miss", _fmt(missing), "FP", _fmt(fp), "FN", _fmt(fn)])
    out_fasta = f"{stem}_{tag}.fasta"
    out_log   = f"{stem}_{tag}_log.tsv"
    out_cells = f"{stem}_{tag}_cells.tsv"

    with open(out_fasta, "w", newline="\n") as fh:
        for c, row in zip(cells, noisy):
            fh.write(f">{c}\n{''.join(row)}\n")
        fh.write(f">{normal_name}\n{normal}\n")

    kind = np.full(mat.shape, "", dtype="<U4")
    kind[is_fn], kind[is_fp], kind[is_miss] = "FN", "FP", "MISS"
    with open(out_log, "w", newline="\n") as fh:
        fh.write("cell\tsite\ttype\ttrue\tobserved\n")
        for i, j in zip(*np.nonzero(kind)):
            fh.write(f"{cells[i]}\t{j + 1}\t{kind[i, j]}\t{mat[i, j]}\t{noisy[i, j]}\n")

    observed_frac = is_miss.mean(axis=1)
    with open(out_cells, "w", newline="\n") as fh:
        fh.write("cell\ttarget_missing_rate\tobserved_missing_rate\n")
        for c, target, obs in zip(cells, p_miss.mean(axis=1), observed_frac):
            fh.write(f"{c}\t{target:.4f}\t{obs:.4f}\n")

    seen  = ~is_miss
    stats = {"cells": n_cells, "sites": n_sites, "missing": int(is_miss.sum()),
             "TP": int((seen & mutated & ~is_fn).sum()), "FN": int(is_fn.sum()),
             "TN": int((seen & ~mutated & ~is_fp).sum()), "FP": int(is_fp.sum())}
    print(f"[noise] {n_cells} cells x {n_sites} sites  (seed={seed}, cell spread={_fmt(cell_spread)}, "
          f"site spread={_fmt(site_spread)})")
    print(f"[noise] missing {stats['missing']} (per cell {observed_frac.min():.0%}-{observed_frac.max():.0%})  |  "
          f"TP {stats['TP']}  FN {stats['FN']}  |  TN {stats['TN']}  FP {stats['FP']}")
    for path in (out_fasta, out_log, out_cells):
        print(f"[noise] Saved: {path}")
    return (out_fasta, out_log, out_cells), stats


def main():
    ap = argparse.ArgumentParser(description="Add missing data, false positives and "
                                             "false negatives to a single-cell SNV FASTA.")
    ap.add_argument("fasta")
    ap.add_argument("missing", help="average missing-data rate (0-1)")
    ap.add_argument("fp", help="false-positive rate at unmutated sites (0-1)")
    ap.add_argument("fn", help="false-negative (dropout) rate at mutated sites (0-1)")
    ap.add_argument("normal", help="name of the normal (unmutated) sequence")
    ap.add_argument("--seed", type=int, default=None, help="random seed for reproducible noise")
    ap.add_argument("--cell-spread", type=float, default=0.0,
                    help="variation in missing rate between cells (0 = uniform, <1)")
    ap.add_argument("--site-spread", type=float, default=0.0,
                    help="variation in missing rate between sites (0 = uniform, <1)")
    ap.add_argument("--out-dir", default=None, help="output folder (default: noise/ next to the input)")
    args = ap.parse_args()
    try:
        add_noise(args.fasta, args.missing, args.fp, args.fn, args.normal,
                  seed=args.seed, cell_spread=args.cell_spread, site_spread=args.site_spread,
                  out_dir=args.out_dir)
    except ValueError as e:
        sys.exit(f"[ERROR] {e}")


if __name__ == "__main__":
    main()
