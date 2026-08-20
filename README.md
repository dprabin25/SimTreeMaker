# SimTreeMaker

A pipeline for running SLiM cancer-evolution simulations and converting the resulting tree sequences into Newick phylogenies and PNG plots.

## Cancer simulation pipelines — overview

This project has four methods across three separate SLiM-based pipelines, each in its own self-contained folder. None of them share state — editing or running one never affects the others.

| Folder | Method | What it does |
|---|---|---|
| `TumorGrowth\` | 1 & 2 | Tumor growth, no migration (1) or with migration (2) |
| `Trees\` | 3 | Convert existing `.trees` files to Newick + SNV + PNG, no SLiM run |
| `CaseStudy\` | 4 | Patient-calibrated clonal hematopoiesis case study |

All three ultimately produce the same three deliverables from a SLiM `.trees` file: a **Newick tree** (`.nwk`), a **per-cell SNV CSV**, and **tree PNGs** (with and without tip labels). That conversion step is done by the same `slim_utilities.py` functions in every folder (each folder keeps its own copy, so the folders stay independent).

Unfamiliar term? See `GLOSSARY.md` — it defines every concept and every Eidos/SLiM code term (functions, properties, callbacks) actually used in each method's `.slim` file.

---

## Dependencies

- **SLiM** (the simulation engine) — download from [messerlab.org/slim](https://messerlab.org/slim/), then point `SLIM_EXE` in each folder's `slim_config.txt` to the executable. (`TumorGrowth\` and `CaseStudy\` each ship a `slim_config.example.txt` template — copy it to `slim_config.txt` in the same folder and edit the path. `Trees\` doesn't run SLiM at all, so it has no `slim_config.txt`.)
- **Python + R packages:**

  ```
  conda create -n SimTreeMaker -c conda-forge python=3.11 biopython=1.85 matplotlib pyslim=1.1.1 tskit=1.0.3 r-base r-ape -y
  conda activate SimTreeMaker
  ```

  R with the `ape` package is optional but preferred for tree plotting (`plot_tree.R`, run automatically via `Rscript` if it's on your PATH); if `Rscript` isn't found, the pipeline falls back to a built-in matplotlib renderer automatically — no action needed either way.

---

## Method 1 — Tumor growth, no migration

**Folder:** `TumorGrowth\`
**Model file:** `TG.slim`
**Runner:** `simtreemaker.py`

### At a glance

- **Working dir:** `TumorGrowth\`
- **What it does:** grows a tumor from a single founder cell to a target size, with spontaneous background and driver mutations, no migration.
- **Model file:** `TG.slim`
- **Requirements:** SLiM installed, its path known. Python 3. `Rscript` on PATH is optional (nicer tree PNGs); falls back to matplotlib automatically if missing.
- **Config:** `slim_config.txt` (SLiM executable path) + `config.txt` (model parameters).
- **Outputs:** `TumorGrowth\TG\NoMigration\` — `.trees`, `newick\`, `snv\`, `pngTree\`, `log.txt` (see full list below).

### The model

`TG.slim` is a nonWF (non-Wright-Fisher) SLiM model of tumor growth from a single founder cell:

- Generation 1 starts with `initialPopSize` founder cell(s) in population `p1`.
- Every living cell clones exactly one offspring per generation (`reproduction()`) — cells are never explicitly killed. Population size is instead regulated by density-dependent viability selection: each generation, `p1.fitnessScaling = targetSize / n`, where `n` is the current population size. This pulls the population toward `targetSize` without a hard cutoff.
- Background (passenger) mutations — type `m1`, neutral — arise at `mutationRate` per base across genome positions 0–49999.
- Driver mutations — type `m2`, selection coefficient `m2_effect`, dominance `m2_dominance` — arise spontaneously at `driverMutationRate` per base across positions 50001–99999. There is no guaranteed founder driver; every `m2` mutation, including the first, appears on its own via this background rate. Because drivers can arise anywhere in that range and at any time, this is a **multi-driver** design — multiple independent driver clones can appear and compete.
- The run stops at generation `simulationEndTick`, at which point the tree sequence is written to `treesPath`.

Required `-d` constants: `initialPopSize`, `targetSize`, `mutationRate`, `driverMutationRate`, `m2_dominance`, `m2_effect`, `simulationEndTick`, `treesPath`.

### Running it directly (no Python)

```
cd TumorGrowth
slim -d initialPopSize=1 -d targetSize=150 -d mutationRate=1e-8 -d driverMutationRate=1e-8 -d m2_dominance=0.5 -d m2_effect=0.6 -d simulationEndTick=200 -d treesPath="'TG.trees'" TG.slim
```

This only ever runs the plain no-migration model — `TG.slim` has no migration code in it at all (see Method 2).

### Running it through `simtreemaker.py` (recommended)

```
cd TumorGrowth
python simtreemaker.py
```

This does everything the terminal command above does, plus the Newick/SNV/PNG conversion afterward, by reading two config files instead of you typing `-d` flags by hand:

- **`slim_config.txt`** — one setting, `SLIM_EXE`, the path to your `slim` executable (on this machine: `C:\Users\newfaculty\anaconda3\envs\SimTreeMaker\Library\bin\slim.exe`, from the `SimTreeMaker` conda environment).
- **`config.txt`** — the actual model parameters: `mutationRate`, `driverMutationRate`, `m2_dominance`, `m2_effect`, `initialPopSize`, `targetSize`, `Pop1SampleCells` (how many `p1` cells get sampled for the tree/SNV output), `simulationEndTick`, and `treeOutputFile` (the base filename for outputs). It also holds the `Pop2`–`Pop4` migration settings and the `migrationEnabled` flag — those only matter for Method 2 (see below); with `migrationEnabled = F` (the default), this run ignores all of them.

You can override any single `config.txt` value for one run without editing the file, e.g.:

```
python simtreemaker.py targetSize=300
```

### Output

Everything lands in `TumorGrowth\TG\NoMigration\`:

- `TumorGrowth.trees` — the raw tree sequence
- `newick\TumorGrowth.nwk` — the lineage tree in Newick format
- `snv\TumorGrowth_snv.csv` — per-cell mutation table (cell ID, population, mutation type `m1`/`m2`, genomic position, generation the mutation arose)
- `pngTree\TumorGrowth_horizontal_labels.png` / `_no_labels.png` — rendered tree images (via R+`ape` if `Rscript` is on PATH, otherwise a matplotlib fallback)
- `log.txt` — one entry per run: the exact SLiM command used, SLiM's console output, and the resulting output paths (or the failure reason, if it failed)

---

## Method 2 — Tumor growth, with migration

**Folder:** `TumorGrowth\` (same folder as Method 1)
**Model file:** `TG.slim` (same file — see below)
**Runner:** `simtreemaker.py migration`

### At a glance

- **Working dir:** `TumorGrowth\` (same folder as Method 1)
- **What it does:** the same tumor growth model as Method 1, plus one or more additional populations seeded by migration from an existing population (metastasis-style spread).
- **Model file:** `TG.slim` (same file as Method 1 — no separate migration `.slim` file is ever saved; see below).
- **Requirements:** identical to Method 1 — same SLiM install, same optional `Rscript`.
- **Config:** the exact same `slim_config.txt` and `config.txt` as Method 1 — nothing separate. `config.txt`'s `migrationEnabled` flag and `PopN*` keys (ignored by Method 1) are what this method actually uses.
- **Inputs:** none — same as Method 1, fresh simulation every run.
- **Outputs:** `TumorGrowth\TG\Migration\` — same file set as Method 1, just in its own subfolder so the two methods' outputs never overwrite each other.

```
cd TumorGrowth
python simtreemaker.py migration
```

This is Method 1 plus one extra step, reading the exact same `config.txt` and `slim_config.txt` — nothing about Method 1 changes to support it.

**What's different:** `TG.slim` itself has no migration logic — just an inert 3-line comment marker inside its `early()` block:

```
////////////////////////////////
//Add migration setting if any//
////////////////////////////////
```

When you add `migration` on the command line, `simtreemaker.py` does the following (all inside `simtreemaker.py` itself — no separate generator file):

1. `detect_migration_pops(params)` scans `config.txt` for any `PopN*` keys (`PopN1Time`, `PopN2Time`, …) and collects whichever population numbers it finds — this isn't fixed to any particular count or set of numbers; add a new population by adding its `PopN*` block to `config.txt` and it's picked up automatically, no code changes needed.
2. `build_migration_slim_text(script_dir, pop_ids)` reads `TG.slim`'s text and calls `insert_migration_block()`, which replaces that marker with generated Eidos code — one block per population found in step 1 — that seeds each population from its configured source population at its configured generation, then draws ongoing migrants every generation after that.
3. `main()` writes that generated text to a temporary `.slim` file (`tempfile.NamedTemporaryFile(..., suffix=".slim", delete=False)`), runs SLiM against it, then removes it in a `finally` block right after — so the file exists only for the duration of that one run.

No separate migration `.slim` file is ever saved — everything needed to reproduce the migration model lives in `TG.slim` (the backbone), `simtreemaker.py` (the code that builds and runs the migration block), and `config.txt` (the population settings).

### What the substitution looks like

The marker, exactly as it sits in `TG.slim`:

```
    ////////////////////////////////
    //Add migration setting if any//
    ////////////////////////////////
```

After `python simtreemaker.py migration` builds the temporary script — shown here for one population, `p2` (`config.txt` defining `Pop2Time`, `Pop2Origin`, `Pop2IniPopSize`, `Pop2TargetSize`, `Pop2MigrationRate`). The three marker lines are kept as-is; the generated block is appended right after them:

```eidos
    ////////////////////////////////
    //Add migration setting if any//
    ////////////////////////////////
    if (sim.cycle >= Pop2Time) {
        source2 = sim.subpopulations[sim.subpopulations.id == Pop2Origin];
        if (sim.cycle == Pop2Time) {
            // Founding pull: seed p2 from its configured source population.
            sim.addSubpop("p2", 0);
            if (size(source2) > 0) {
                nMig = min(Pop2IniPopSize, source2.individualCount);
                migrants = sample(source2.individuals, nMig);
                p2.takeMigrants(migrants);
            }
        } else if (size(source2) > 0) {
            // Ongoing migration: draw this generation's migrant count from the source.
            nOngoing = rbinom(1, source2.individualCount, Pop2MigrationRate);
            if (nOngoing > 0)
                p2.takeMigrants(sample(source2.individuals, nOngoing));
        }
        n2 = p2.individualCount;
        if (n2 > 0)
            p2.fitnessScaling = Pop2TargetSize / n2;
    }
```

If `config.txt` defines more than one population (e.g. `Pop2` and `Pop3`), one of these blocks is appended per population found — `Pop2Time`/`Pop2Origin`/etc. for the first, `Pop3Time`/`Pop3Origin`/etc. for the second, back to back after the same marker. `Pop2Time`, `Pop2Origin`, and the rest are Eidos identifiers here — they arrive as `-d` constants on the SLiM command line, built from the matching `config.txt` values. This spliced text only ever exists as a temporary file for the duration of one run (see step 3 above) — never saved permanently.

### Migration population settings in `config.txt`

Each population is a number `N` (any number ≥ 2 — not limited to 2/3/4, and not necessarily consecutive) with five keys:

```
PopN Time            = ...   # generation this population is created and seeded
PopN IniPopSize      = ...   # migrants pulled in at PopN Time
PopN TargetSize      = ...   # carrying capacity thereafter
PopN MigrationRate   = ...   # per-source-individual probability of migrating in, each generation after PopN Time
PopN SampleCells     = ...   # cells sampled from this population for the Newick tree + SNV CSV
PopN Origin          = ...   # which existing population this one migrates from (e.g. "Pop1", "p1", or a bare integer)
```

(written without the space between `Pop` and `N` in the actual file, e.g. `Pop2Time`, `Pop5Origin`.) Override any of these the same way as any other config value, e.g. to move when population `2` appears:

```
python simtreemaker.py migration Pop2Time=30
```

**Output:** same structure as Method 1, just in `TumorGrowth\TG\Migration\` instead of `TG\NoMigration\` — `TumorGrowth.trees`, `newick\`, `snv\`, `pngTree\`, `log.txt`, same file-naming convention, same log format.

---

## Method 3 — Process existing trees (no SLiM run)

**Folder:** `Trees\`
**Runner:** `process_trees.py`

### At a glance

- **Working dir:** `Trees\`
- **What it does:** converts existing `.trees` files to Newick + SNV CSV + tree PNGs. No simulation runs here.
- **Model file:** none — no SLiM run happens in this method at all.
- **Requirements:** Python 3 only — no SLiM installation needed at all for this method.
- **Config:** none — sampling defaults (all populations, up to 400 cells, no fixed seed) are fixed directly in `slim_utilities.py`'s `load_tree_options()`.
- **Inputs:** **required** — one or more `.trees` files must already exist directly inside `Trees\` before running this (from any of the other methods, or from elsewhere). With none present, it just reports that and exits.
- **Outputs:** `Trees\newick\`, `Trees\snv\`, `Trees\pngTree\`, `Trees\log.txt` — one set per `.trees` file, named after that file.

For when you already have `.trees` files — from any of the methods above, or from elsewhere — and just want the Newick/SNV/PNG conversion without running SLiM again.

```
cd Trees
python process_trees.py
```

Drop any `.trees` file into this folder and running the command scans for every `*.trees` file and converts each one it hasn't already processed (checked by whether a matching `newick\<name>.nwk` already exists). Re-run any time — already-done files are skipped automatically.

- `python process_trees.py --force` — reprocess everything, even already-done files
- `python process_trees.py someFile.trees` — process just that one file

**Output:** `newick\`, `snv\`, `pngTree\`, `log.txt`, directly inside `Trees\` — same naming convention as the other methods, just keyed off each `.trees` file's own name instead of a run's parameters.

---

## Method 4 — Case study (clonal hematopoiesis)

**Folder:** `CaseStudy\`
**Model file:** `Simulation4_ClonalHematopoiesis.slim`
**Runner:** `casestudy.py`

### At a glance

- **Working dir:** `CaseStudy\`
- **What it does:** grows a single hematopoietic stem cell to 100,000 cells, with a driver mutation deliberately introduced at a chosen generation, modeling clonal hematopoiesis.
- **Model file:** `Simulation4_ClonalHematopoiesis.slim` (separate from `TG.slim` — not related to Methods 1/2).
- **Requirements:** SLiM installed, its path known. Python 3. Same optional `Rscript` for nicer PNGs as the other methods.
- **Config:** `slim_config.txt` only (SLiM executable path). No `config.txt` — `AGE`/`DRIVER_TICK`/`REP` have defaults built into `casestudy.py` itself and are overridden straight from the terminal.
- **Outputs:** `CaseStudy\CaseStudyOutputs\` — `.trees`, `newick\`, `snv\`, `pngTree\`, `log.txt` (see full list below).

A separate, patient-calibrated model (modeled after Mitchell et al. 2022, *Nature*) — a single ancestral hematopoietic stem cell expands to 100,000 cells, with a driver mutation introduced into a random newborn cell at a configurable generation, and stem-cell attrition after 35 generations approximating age-related decline.

Three `-d` constants, edited straight from the terminal — no `config.txt` for this one, since there's only one model file and three values:

- `AGE` — generation the simulation stops at (the patient's age at sampling; ≤ 500)
- `DRIVER_TICK` — generation the driver mutation is introduced
- `REP` — a replicate label used only to keep filenames from repeated runs apart (it does **not** seed SLiM's RNG — two runs with the same `AGE`/`DRIVER_TICK`/`REP` aren't guaranteed identical results)

```
cd CaseStudy
python casestudy.py                              # defaults: AGE=81 DRIVER_TICK=38 REP=1
python casestudy.py AGE=85 DRIVER_TICK=40 REP=2   # override any/all
```

`slim_config.txt` (same `SLIM_EXE` setting as `TumorGrowth\`) is still read, since SLiM still needs to run.

**Output:** `CaseStudyOutputs\`:

- `Sim4_CHIP<AGE>_rep<REP>.trees`
- `newick\Sim4_CHIP<AGE>_rep<REP>.nwk`
- `snv\Sim4_CHIP<AGE>_rep<REP>_snv.csv`
- `pngTree\Sim4_CHIP<AGE>_rep<REP>_horizontal_labels.png` / `_no_labels.png`
- `log.txt`

---

## Shared conventions across all four

- Every runner prints the equivalent raw terminal `slim` command before running it, so you can see exactly what's being executed.
- Every runner appends one entry per run to its own `log.txt` — the command used, SLiM's console output, and the resulting file paths (or the failure reason).
- Newick tip labels follow `gen0_ind<i>|pop<N>|<mutations>` (e.g. `no_mut`, `m1`, `m2`, or `m1+m2`).
- SNV CSVs have one row per (cell, mutation) pair: `cell_id`, `population`, `mutation_type` (`m1`/`m2`/`none`), `position`, `generations_before_sampling`, `origin_tick`.
- Tree PNGs render via R + `ape` (`plot_tree.R`) if `Rscript` is on your PATH, otherwise fall back to a matplotlib-based renderer — both produce a labeled and an unlabeled horizontal tree image.
- How many cells get sampled for the tree/SNV output: `TumorGrowth\config.txt`'s `Pop1SampleCells` (and `PopN SampleCells` for migration populations) overrides the pooled default per population. Everywhere else — the pooled fallback in `TumorGrowth\`, and all sampling in `CaseStudy\`/`Trees\` — uses the fixed defaults in `slim_utilities.py`'s `load_tree_options()` (all populations, up to 400 cells, no fixed seed); edit that function directly to change them.
