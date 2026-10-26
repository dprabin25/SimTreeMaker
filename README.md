# SimTreeMaker

**Simulate tumors whose true evolutionary history is known, and use them to test single-cell phylogeny and genotyping methods.**

How well does a tumor-phylogeny method work? On real data the true tree is never known, so the question cannot be answered. SimTreeMaker grows a virtual tumor cell by cell in [SLiM](https://messerlab.org/slim/), records every cell division, and gives you:
- the **true cell tree**;
- the **true mutations** of every sampled cell;
- the **same data as a sequencing experiment would deliver**, optionally with realistic single-cell errors.

Run your method on the sequences, compare its answer with the truth, and you know how accurate it is.

---

## Why use SimTreeMaker

- **The right answer is known.** Every sampled cell's ancestry and mutations are recorded exactly.
- **Realistic tumor biology.** Driver mutations with positive selection, clonal expansions, metastatic spread between sites, and drug treatment with resistance.
- **Realistic single-cell data.** Missing data, allelic dropout and false-positive calls, at rates you choose. A noise-free control is produced automatically.
- **No programming needed.** Everything is set in one plain-text file, `config.txt`.
- **Fast and reproducible.** A typical run takes seconds on a laptop. One random seed reproduces a run exactly.

### What you can use it for

- Benchmarking tumor-phylogeny methods on data where the true tree is known.
- Testing how genotype-calling and error-correction methods cope with missing data and sequencing errors.
- Studying how selection, tumor age and size, metastasis or therapy shape cell phylogenies.
- Generating example data for teaching and for developing new methods.

---

## Getting started (3 steps)

**1. Install** SLiM, Python and R in one conda environment:

```
conda create -n SimTreeMaker -c conda-forge python=3.11 slim=5.2 tskit=1.0.3 pyslim=1.1.1 numpy biopython matplotlib r-base r-ape -y
conda activate SimTreeMaker
```

**2. Tell SimTreeMaker where SLiM is.** Open `SLiMcancerCodes/slim_config.txt` and set the path, which you can find with `which slim` (macOS/Linux) or `where slim` (Windows):

```
SLIM_EXE = /path/to/slim
```

**3. Run:**

```
cd SLiMcancerCodes
python simtreemaker.py nomigration
```

SLiM may print a warning that "a neutral mutation type was defined and used"; this is expected and can be ignored. On Windows, keep the SimTreeMaker folder at a short path (e.g. `C:\SimTreeMaker` or `Documents\SimTreeMaker`), because Windows cannot create files whose full path is longer than 260 characters. A finished run ends with `Run complete.`

The results are in `SLIMcanceroutputs/TG/NoMigration/`:

```
NoMigration/
├── newick/TumorGrowth.nwk          true cell tree
├── pngTree/                        pictures of the tree
├── snv/TumorGrowth_snv.csv         true mutations of every cell
├── fasta/TumorGrowth_snv.fasta     sequences of the sampled cells (input for your method)
└── log.txt                         settings, seed and messages of every run
```

`SLIMcanceroutputs.zip` in this repository contains this example output for all modes, so you can look at it before installing anything.

---

## What you can simulate

Each feature is switched on with a word on the command line, and controlled by its own section of `config.txt`.

| Feature | Command | config.txt section |
|---|---|---|
| Tumor growth (single tumor) | `python simtreemaker.py nomigration` | Primary population, Mutation, Genome layout |
| Metastasis | `python simtreemaker.py migration` | Migration |
| Therapy (single tumor only) | `python simtreemaker.py nomigration treatment` | Treatment |
| Sequencing errors | add `noise` to any of the commands above | Sequencing noise |

Therapy is available only for a single tumor (`nomigration`). Migration runs are never treated.

### Tumor growth

A tumor grows from one founder cell (or more). Each generation, cells divide. Crowding keeps cells without drivers near a carrying capacity, while driver mutations make a cell fitter, so its descendants expand into large clones. Neutral (passenger) and driver mutations accumulate at the rates you set.

*Change:* tumor age, size, mutation rates, driver strength, genome size.

### Metastasis (migration)

Metastatic sites are founded by cells that migrate from the primary tumor or from another metastasis. They then keep receiving migrants at a rate you set. Each site grows and evolves on its own and is sampled separately, and all cells end up in one tree.

*Change:* the number of sites, when each is founded, its source, its size, and its migration rate.

### Therapy (treatment)

A single drug dose is given to a growing single tumor (no migration). Before the dose, resistance mutations arise by chance; they carry a small fitness cost. The dose kills most sensitive cells, and resistant cells survive. The tumor then regrows from the survivors. Resistant cells are marked in all outputs, and the numbers of sensitive and resistant cells in every generation are saved and plotted.

*Change:* when the dose is given, how strong it is, the cost of resistance, and the resistance mutation rate.

### Sequencing errors (noise)

The sequences are degraded as in single-cell DNA sequencing:
- missing data (`?`);
- false negatives (allelic dropout);
- false positives.

The amount of missing data can be the same for all cells or vary between cells and sites. Two versions are always written: a noise-free **control**, checked automatically against the true data, and the noisy **test**. Every introduced error is logged.

*Change:* the three error rates and how unevenly missing data is spread.

---

## Changing the settings (`config.txt`)

Everything about a simulation is set in `SLiMcancerCodes/config.txt`; no code needs to be edited.
1. Open it in any text editor (Notepad, TextEdit, …).
2. Change the value after `=` and save.
3. Run again.

The text after `#` on each line explains the setting. For example, to simulate an older tumor, change

```
simulationEndTick   = 200
```

to

```
simulationEndTick   = 300
```

Settings are checked before the simulation starts; a value that does not make sense (for example, overlapping genome regions) stops the run with a message saying what to fix. All settings are listed below. `SLiMcancerCodes/CONFIG_GUIDE.md` explains each one in more detail.

### Run

| Setting | Default | What it changes |
|---|---|---|
| `simulationEndTick` | 200 | Number of generations (tumor age). All populations are sampled at this generation. |
| `simSeed` | 12345 | Random seed. Same seed and settings = same tumor. Change it for a new replicate; leave empty for a random seed. |

### Primary population (p1)

| Setting | Default | What it changes |
|---|---|---|
| `initialPopSize` | 1 | Number of founder cells (1 = growing tumor; equal to `targetSize` = population of stable size). |
| `targetSize` | 500 | Carrying capacity for cells without drivers. Tumors with drivers grow beyond it. |
| `Pop1SampleCells` | 400 | Number of cells sampled. |

### Mutation

| Setting | Default | What it changes |
|---|---|---|
| `mutationRate` | 1e-8 | Neutral mutation rate (per base per generation). |
| `driverMutationRate` | 1e-8 | Driver mutation rate. 0 = no drivers. |
| `m2_effect` | 0.3 | Strength of selection on drivers (*s*). 0 = drivers have no effect (neutral evolution). |
| `m2_dominance` | 0.5 | Dominance of drivers (*h*). Each driver multiplies fitness by 1 + *hs* (1.15 by default). |

### Genome layout

| Setting | Default | What it changes |
|---|---|---|
| `m1Start`, `m1End` | 0 – 9,999,999 | Neutral region (10 Mb). |
| `m2Start`, `m2End` | 10,000,000 – 10,049,999 | Driver region (50 kb). |

Positions are counted from 0, ends included, and regions must not overlap. A longer region gives more mutations: doubling its length has the same effect as doubling the mutation rate.

### Migration (used with `migration`)

Each row is one metastatic site:

```
#        Origin  Time  IniPopSize  TargetSize  MigrationRate  SampleCells
Pop2   = Pop1    50    5           100         0.01           100
Pop3   = Pop1    100   5           80          0.01           100
Pop4   = Pop3    125   5           150         0.01           100
```

| Column | Meaning |
|---|---|
| Origin | Where the migrants come from (`Pop1` = primary tumor, or an earlier site). |
| Time | Generation the site is founded (before `simulationEndTick`). |
| IniPopSize | Number of founding cells. |
| TargetSize | Carrying capacity of the site. |
| MigrationRate | Chance, per cell and generation, that a cell of Origin migrates in afterwards. |
| SampleCells | Number of cells sampled from the site. |

- To add a site, add a row (`Pop5 = …`).
- To remove one, delete its row or put `#` at the start of it.

### Treatment (used with `nomigration treatment`; single tumor only)

| Setting | Default | What it changes |
|---|---|---|
| `treatment` | off | `on` = every no-migration run is treated (same as adding `treatment` to the command). |
| `treatmentTargetSize` | 2500 | Carrying capacity in treatment runs. It is larger than usual so that rare resistant cells are present when the drug is given. |
| `treatmentTick` | 150 | Generation of the dose. |
| `drugKill` | 0.7 | Chance that a sensitive cell dies at the dose. |
| `resistanceCost` | 0.05 | Fitness cost of resistance (0.05 = 5% less fit; 0 = no cost). |
| `resistanceMutationRate` | 2e-8 | Resistance mutation rate. |
| `m3Start`, `m3End` | 10,050,000 – 10,054,999 | Resistance region (5 kb). |

### Sequencing noise (used with `noise`)

| Setting | Default | What it changes |
|---|---|---|
| `noise` | off | `on` = every run adds noise (same as adding `noise` to the command). |
| `missingRate` | 0.3 | Average fraction of missing sites. |
| `cellSpread` | 0 | How much missing data varies between cells (0 = equal; closer to 1 = very unequal). |
| `siteSpread` | 0 | How much missing data varies between sites. |
| `fnRate` | 0.2 | False-negative (dropout) rate at mutated sites. |
| `fpRate` | 0.01 | False-positive rate at unmutated sites. |
| `noiseSeed` | 12345 | Seed for the noise. |

### Output

| Setting | Default | What it changes |
|---|---|---|
| `outputDir` | `../SLIMcanceroutputs` | Where results are written. |
| `treeOutputFile` | TumorGrowth | Name given to all output files. |
| `fastaSeed` | 12345 | Seed for the DNA bases in the FASTA files. |
| `writeFullFasta` | F | `T` = also write the full genome of every sampled cell (large). |

---

## Output files

Results go to `<outputDir>/TG/NoMigration/`, `Migration/` or `Treatment/`.

| File | Contents |
|---|---|
| `newick/<name>.nwk` | True cell tree of the sampled cells (Newick) |
| `pngTree/<name>_horizontal_labels.png`, `_no_labels.png` | Pictures of the tree |
| `snv/<name>_snv.csv` | One row per cell and mutation: population, type (m1 neutral, m2 driver, m3 resistance), position, generation it arose |
| `fasta/<name>_snv.fasta` | Sequences of the sampled cells at all mutated sites (one row per cell) |
| `fasta/<name>_snv_positions.tsv` | Genome position, mutation type and the two DNA bases of each column |
| `fasta/<name>_parent.fasta` | Genome of the founder cell |
| `fasta/<name>_full.fasta` | Full genome of every sampled cell (only if `writeFullFasta = T`) |
| `noise/control/`, `noise/test/` | Noise-free and noisy sequences, each with a Normal (unmutated) sequence, a log of every change and the missing rate of each cell |
| `trajectory/<name>_trajectory.csv`, `.png` | Treatment runs: total, sensitive and resistant cells in every generation |
| `<name>.trees` | Complete simulation history (tskit tree sequence) for advanced analyses |
| `log.txt` | Settings, seed and messages of every run |

Cells have the same name in every file, e.g. `gen0_ind12|pop1|m1+m2`. The fields are the cell ID, its population, and the types of mutation it carries.

---

## Reproducibility

- `simSeed` controls the simulation and which cells are sampled.
- `noiseSeed` controls the noise.
- `fastaSeed` controls the DNA bases.

With the same `config.txt` and command, the result is identical. The seed of every run is printed at the end and saved in `log.txt`. For replicates, change only `simSeed` (e.g. 1, 2, 3, …).

---

## How the code works

`simtreemaker.py` is the only program you run. It reads `config.txt`, chooses the SLiM model, runs SLiM, then hands the result to Python code that samples cells and writes the outputs. The table shows which code is responsible for each feature. Line numbers refer to the current version.

| Feature | What decides it | Where in the code |
|---|---|---|
| Tumor growth | Always | `TG.slim` (whole file) |
| Migration | `migration` on the command line | `TG.slim` lines 94–96 (marker) + `simtreemaker.py` `_POP_TEMPLATE` (line 97), `build_migration_slim_text()` (line 220) |
| Treatment | `treatment` on the command line, or `treatment = on` | `simtreemaker.py` line 612 selects `TG_Treatment.slim`; dose in `TG_Treatment.slim` lines 103–109 |
| Cell sampling and tree | Always | `slim_utilities.py` `convert_to_newick()` (line 550), `_merge_genome_copies()` (line 387) |
| FASTA | Always | `slim_utilities.py` `write_fasta_set()` (line 226) |
| Noise (gaps and errors) | `noise` on the command line, or `noise = on` | `simtreemaker.py` lines 801–814 call `add_noise.py` `add_noise()` (line 104) |

### Tumor growth: `TG.slim`

`TG.slim` is a plain, single-population model. The three key parts:

```
reproduction() {
    subpop.addCloned(individual);          // every cell makes one copy of itself
}

early() {
    n = p1.individualCount;
    p1.fitnessScaling = targetSize / n;    // crowding: survival = fitness x targetSize / n
    ////////////////////////////////
    //Add migration setting if any//       <- marker, lines 94-96 (does nothing on its own)
    ////////////////////////////////
}

s1 2 late() {
    sim.treeSeqOutput(treesPath);          // save the full ancestry at simulationEndTick
}
```

Drivers act through `initializeMutationType("m2", m2_dominance, "f", m2_effect)` (line 55): SLiM multiplies a cell's fitness by 1 + *hs* for each driver it carries.

### Migration: code inserted at the marker

`TG.slim` itself never migrates. With `migration`, `simtreemaker.py` does three things:
1. It reads `TG.slim`.
2. It replaces the 3-line marker with one block of code per `PopN` row in `config.txt` (`insert_migration_block()`, line 195). The block is written by `_POP_TEMPLATE` (line 97).
3. It runs the result from a temporary file (line 632), which is deleted after the run.

For `Pop2` the inserted block is:

```
if (sim.cycle >= Pop2Time) {
    source2 = sim.subpopulations[sim.subpopulations.id == Pop2Origin];
    if (sim.cycle == Pop2Time) {                              // founding
        sim.addSubpop("p2", 0);
        p2.takeMigrants(sample(source2.individuals, min(Pop2IniPopSize, source2.individualCount)));
    } else {                                                  // ongoing migration
        nOngoing = rbinom(1, source2.individualCount, Pop2MigrationRate);
        if (nOngoing > 0)
            p2.takeMigrants(sample(source2.individuals, nOngoing));
    }
    p2.fitnessScaling = Pop2TargetSize / p2.individualCount;  // crowding in p2
}
```

So **adding a `Pop5` row to `config.txt` adds one more block like this**; no code is edited by hand.

### Treatment: `TG_Treatment.slim`

With `treatment`, line 612 of `simtreemaker.py` runs `TG_Treatment.slim` instead of `TG.slim`, and lines 670–672 pass the treatment settings to SLiM. `TG_Treatment.slim` is `TG.slim` plus:

```
initializeMutationType("m3", 1.0, "f", -resistanceCost);    // line 64: resistance, with a fitness cost

if (sim.cycle == treatmentTick) {                           // lines 103-109: the single dose
    inds = p1.individuals;
    sensitive = inds[inds.countOfMutationsOfType(m3) == 0];  // cells without m3
    sensitive.fitnessScaling = 1.0 - drugKill;               // survive with probability 1 - drugKill
}
```

The `late()` block (lines 112–121) writes the numbers of sensitive and resistant cells in every generation to the trajectory file, and stops the run if no cell is left.

### Cell sampling and the true tree: `slim_utilities.py`

`convert_to_newick()` (line 550) opens the `.trees` file and samples `PopNSampleCells` cells per population (seeded by `simSeed`). It then calls `_merge_genome_copies()` (line 387), which joins each cell's two genome copies so that mutations on either copy are reported. Finally it writes the Newick tree and the SNV table.

### FASTA: `write_fasta_set()` in `slim_utilities.py` (line 226)

SLiM records where mutations are, not which DNA base they create, so the bases are drawn here (seeded by `fastaSeed`):

```
parent = rng.integers(0, 4, size=seq_len)                  # line 240: random founder genome (A/C/G/T)
allele_base = [ref_base] + [alt_pool[i % 3] ...]           # line 254: a different random base for each mutation
```

Each column of `_snv.fasta` is one mutated site: a cell shows the founder base if it lacks the mutation, and the new base if it has it. `_snv_positions.tsv` records the genome position and both bases of each column.

### Noise and gaps: `add_noise()` in `add_noise.py` (line 104)

`simtreemaker.py` calls `add_noise()` twice (lines 805–808): once with all rates 0 (control) and once with the `config.txt` rates (test). It then checks that the control equals the clean FASTA (line 814). Inside `add_noise()`:

```
cell_rate = beta_rates(rng, missing, cell_spread, n_cells)     # line 145: missing rate of each cell
site_rate = beta_rates(rng, missing, site_spread, n_sites)     #           and of each site
is_miss = rng.random(mat.shape) < p_miss                       # line 156: gaps
is_fn   = ~is_miss & mutated  & (flip < fn)                    # line 158: false negatives (mutant -> normal base)
is_fp   = ~is_miss & ~mutated & (flip < fp)                    # line 159: false positives (normal -> mutant base)
noisy[is_miss] = "?"                                           # line 164: gaps written as "?"
```

Gaps are decided first, so false negatives and false positives only affect sites that are observed. Every change is written to the `_log.tsv` file.

---

<details>
<summary><b>Model details</b></summary>

- **Growth.**
  - SLiM non-Wright–Fisher model of a clonal population of diploid cells.
  - Each generation, every cell produces one daughter. Each cell survives with probability min(1, fitness × `targetSize` / N), where N is the population size.
  - Drivers (m2) multiply fitness by 1 + *hs*; neutral mutations (m1) do not change fitness.
  - Mutations arise at per-base rates in the regions set in Genome layout.
  - The ancestry of all cells is recorded as a tree sequence.
- **Migration.**
  - The migration code is generated from `TG.slim` for each run, using the `PopN` rows, and runs from a temporary file.
  - Each site is created at its `Time` with `IniPopSize` migrants from its `Origin`. Afterwards, each Origin cell migrates in with probability `MigrationRate` each generation.
  - Each site has its own carrying capacity.
- **Treatment** (`TG_Treatment.slim`, single tumor only).
  - Resistance mutations (m3) have selection coefficient −`resistanceCost` and dominance 1. A cell with at least one m3 is resistant.
  - At `treatmentTick`, each sensitive cell dies with probability `drugKill`; resistant cells are unaffected. No further drug is given.
  - If no cell survives, the run reports `TUMOR ERADICATED` and writes no tree.
- **Sampling and sequences.**
  - `PopNSampleCells` cells are sampled from each population at `simulationEndTick`, and mutations on both genome copies are reported.
  - SLiM does not assign DNA bases. SimTreeMaker therefore draws a random founder genome and, for each mutated site, a random different base.
- **Noise.**
  - Missing sites are hidden first. False negatives (mutant → founder base) and false positives (founder → mutant base) are then applied to the observed sites.
  - With `cellSpread`/`siteSpread` > 0, each cell's or site's missing rate is drawn from a beta distribution with mean `missingRate`.
  - The control (all rates 0) must equal the clean alignment, or the run stops.

</details>

## Limitations

- Cells reproduce clonally. Only point mutations (SNVs) are modeled; copy-number changes and structural variants are not.
- There is no spatial structure within a population; metastatic sites are separate populations.
- Treatment is a single dose in a single population, with one resistance mechanism.
- Very large tumors (≫10⁴ cells) and long runs take more time and memory.

---

## Repository contents

| Path | Contents |
|---|---|
| `SLiMcancerCodes/simtreemaker.py` | Main program: run this |
| `SLiMcancerCodes/config.txt` | All settings |
| `SLiMcancerCodes/CONFIG_GUIDE.md` | Detailed explanation of every setting |
| `SLiMcancerCodes/slim_config.txt` | Location of SLiM on your computer |
| `SLiMcancerCodes/TG.slim`, `TG_Treatment.slim` | SLiM models (without and with treatment) |
| `SLiMcancerCodes/slim_utilities.py`, `add_noise.py` | Sampling, tree, sequence and noise code |
| `SLiMcancerCodes/plot_tree.R`, `plot_trajectory.R` | Plots (R + ape; matplotlib is used if R is missing) |
| `SLIMcanceroutputs.zip` | Example output of all modes |
| `GLOSSARY.md` | Definitions of terms |

`add_noise.py` can also add noise to any FASTA alignment that has a normal sequence: `python add_noise.py alignment.fasta 0.3 0.01 0.2 Normal --seed 1`. The numbers are the missing-data, false-positive and false-negative rates; see `python add_noise.py -h`.

Tested with SLiM 5.2, Python 3.11, tskit 1.0.3, R 4.5 and ape 5.8 on Windows.

---

## Citation

If you use SimTreeMaker, please cite:

Dawadi P, Miura S. SimTreeMaker: a flexible SLiM-based simulator for tumor evolution and cellular phylogeny benchmarking. *Manuscript in preparation.*

and the software it builds on:
- Haller BC, Messer PW (2023). SLiM 4: multispecies eco-evolutionary modeling. *The American Naturalist* 201(5):E127–E139.
- Paradis E, Schliep K (2019). ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. *Bioinformatics* 35(3):526–528.

## Contact

Questions, bugs and suggestions: please open an [issue](https://github.com/dprabin25/SimTreeMaker/issues), or contact Prabin Dawadi (pdawadi@go.olemiss.edu) or Sayaka Miura (smiura@olemiss.edu), University of Mississippi.
