# SimTreeMaker

A pipeline for running SLiM cancer-evolution simulations and converting the resulting tree sequences into Newick phylogenies and PNG plots.

## Dependencies

**SLiM** (the simulation engine) — download from https://messerlab.org/slim/, then point `SLIM_EXE` in `slim_config.txt` to the executable.

**Python + R packages:**

```bash
conda create -n SimTreeMaker -c conda-forge python=3.11 biopython=1.85 matplotlib pyslim=1.1.1 tskit=1.0.3 r-base r-ape -y
conda activate SimTreeMaker
```

R with the `ape` package is optional but preferred for tree plotting (`plot_tree.R`, run automatically via `Rscript` if it's on your `PATH`); if `Rscript` isn't found, the pipeline falls back to a built-in matplotlib renderer automatically — no action needed either way.

## File Structure

```
SimTreeMaker/
├── simtreemaker.py                    # Main pipeline script — run this
├── slim_newick.py                     # Converts a .trees tree sequence to Newick (Snapshot sampler)
├── plot_tree.R                        # R/ape tree renderer, called automatically when Rscript is available
├── slim_config.txt                    # SLIM_EXE path — the one file you normally need to edit
├── requirements.txt                   # Pinned Python package versions
├── SimOptions/                        # Simulation parameter files — one .txt per model
│   ├── MutationSpread_NWF.txt / .md
│   ├── MutationSpread_WF.txt  / .md
│   ├── ClonalGrowth_NWF.txt   / .md
│   ├── ClonalGrowth_WF.txt    / .md
│   ├── Metastasis_NWF.txt     / .md
│   └── Metastasis_WF.txt      / .md
├── TreeOptions/
│   └── Snapshot.txt                   # Tree-sampling parameters — final-generation cells only
├── CaseStudy/
│   ├── CHIP2EnvN.slim                 # 3-gene neutral CHIP model (AGE, REP)
│   └── Simulation4_ClonalHematopoiesis.slim   # Patient-calibrated CHIP model (AGE, DRIVER_TICK, REP)
├── CaseStudyOutputs/                  # CaseStudy results (created automatically)
├── ReadyTrees/                        # Drop your own .trees/.tree files here
└── ReadyTreesOutputs/                 # ReadyTrees results (created automatically)
```

Each `.txt` file in `SimOptions/` has a matching `.md` file with the same parameters explained in more depth (SLiM mechanics + plain-language description of what the model simulates) — this README covers the same parameters more concisely, organized by model family.

## How to Run

```bash
python simtreemaker.py <ModelName> [tree]
python simtreemaker.py CaseStudy [tree] [KEY=VALUE ...]
python simtreemaker.py Tree
```

- **`<ModelName>`** — one of the six files in `SimOptions/` (without `.txt`): `MutationSpread_NWF`, `MutationSpread_WF`, `ClonalGrowth_NWF`, `ClonalGrowth_WF`, `Metastasis_NWF`, `Metastasis_WF`.
- **`tree`** — skip re-running SLiM and just rebuild the Newick tree/PNGs from an existing `.tree` file.
- **`KEY=VALUE`** *(CaseStudy only)* — overrides one of the SLiM script's `-d` constants (e.g. `AGE=100`, `REP=1`, `DRIVER_TICK=38`). See **CaseStudy** below — this is required if you want to change patient age, driver-mutation timing, or replicate number; the script will otherwise silently use its built-in defaults.

Examples:

```bash
python simtreemaker.py MutationSpread_NWF
python simtreemaker.py MutationSpread_NWF tree
python simtreemaker.py CaseStudy
python simtreemaker.py CaseStudy AGE=100 REP=1
python simtreemaker.py CaseStudy DRIVER_TICK=50 AGE=90 REP=3
python simtreemaker.py Tree
```

## SimOptions/*.txt — Simulation Parameter Reference

Each `.txt` file configures one cancer-evolution model as `key = value` pairs (one per line; lines starting with `#` are comments). `modelType` selects **WF** (Wright-Fisher: fixed population size, non-overlapping generations, all cells replaced each generation) or **nonWF** (non-Wright-Fisher: explicit birth/death, population size can fluctuate, cells have individual ages).

| Model file | Type | What it simulates |
| --- | --- | --- |
| `MutationSpread_NWF.txt` | nonWF | Small, age-structured population. Driver injected into a newborn; short-lived cells die after `maxAge`. |
| `MutationSpread_WF.txt` | WF | Fixed-size population. Driver seeded once, then manually copied into extra cells for a window of generations to force early spread. |
| `ClonalGrowth_NWF.txt` | nonWF | Clonal expansion with a carrying capacity and near-immortal cells (large `maxAge`); single early driver. |
| `ClonalGrowth_WF.txt` | WF | Population grows exponentially, splits into a second (metastatic) population, and carries two sequential driver mutations. |
| `Metastasis_NWF.txt` | nonWF | Two populations. Primary tumour (p1) grows first; a burst of cells migrates to found the metastatic site (p2) at a single tick. |
| `Metastasis_WF.txt` | WF | Two populations present from generation 1. Driver spreads in p1, then a continuous migration rate moves cells to p2 over a window. |

### Shared parameters (all six models)

| Parameter | Description |
| --- | --- |
| `ModelName` | Output folder name. No spaces — should match the filename. |
| `modelType` | `WF` or `nonWF` (see above). |
| `mutationRate` | Background mutation rate per base per generation. `0` = no random background mutations, driver only. |
| `m1_dominance` | Dominance of the neutral marker mutation `m1`. `0`=recessive, `0.5`=co-dominant, `1`=dominant. |
| `m1_effect` | Fitness effect of `m1`. Always `0.0` — `m1` never affects survival, it's a passenger/marker only. |
| `m2_dominance` | Dominance of the driver mutation `m2`. |
| `m2_effect` | Fitness advantage from the driver. E.g. `0.6` = cells carrying `m2` are 60% more fit. Higher = faster clonal sweep. |
| `m2_convertToSubstitution` | `F` = keep tracking `m2` after it fixes in the population (recommended). `T` = drop it from tracking once fixed. |
| `m1_proportion` | Fraction of genomic sites assigned to the neutral marker type. `0.001` = sparse placeholder; higher values (e.g. `0.4`) give a denser marker useful for visualizing clonal inheritance. |
| `chromosomeEnd` | Last base-pair position of the simulated chromosome (`99999` = 100,000 bp). |
| `recombinationRate` | Recombination rate per base. `0` = clonal inheritance, no crossover — standard for tumour models. |
| `mutationIntroTick` | Generation at which the driver (`m2`) is introduced. Must be less than `simulationEndTick`. |
| `mutationPosition` | Genomic position where `m2` is placed (`0` to `chromosomeEnd`). |
| `simulationEndTick` | Generation at which the simulation stops and the tree sequence is written. |
| `treeOutputFile` | Filename stem for this model's outputs (`.tree`, `.nwk`, PNGs). |

### WF-only parameters

| Parameter | Description |
| --- | --- |
| `targetSize` | Population size at generation 1, held fixed every generation thereafter. |
| `mutationSpreadStart` / `mutationSpreadEnd` | Generation window during which `m2` is manually copied into an extra random cell each tick, to force early spread before selection alone would achieve it. |

### nonWF-only parameters

| Parameter | Description |
| --- | --- |
| `initialPopSize` | Number of cells at generation 1 (`1` = single founder cell, typical for tumour-origin models). |
| `cloneCount` | Daughter cells produced per cell per generation. |
| `targetSize` | Carrying capacity — population is pushed toward this size via fitness scaling as it grows. |
| `maxAge` | Max generations a cell can live before it's forced to die. Low (`2`–`5`) = high turnover; high (`≥100`) = effectively immortal. |

### ClonalGrowth_WF — additional parameters

| Parameter | Description |
| --- | --- |
| `initialPopSize` | Starting size at generation 1, before exponential growth begins. |
| `growthAcc` | Growth accumulator starting value (`1.0` = growth begins at the first growth tick). |
| `growthRate` | Per-generation exponential growth multiplier for p1 (`1.05` = 5%/gen). |
| `growthStartTick` / `growthEndTick` | Window of exponential growth for p1. |
| `populationSplitTick` | Generation when a subset of p1 cells split off to found p2. Must be after `growthEndTick`. |
| `splitPopSize` | Number of cells seeded into p2 at the split (small = tight bottleneck). |
| `secondMutationTick` | Generation when a second driver mutation is seeded in p1. Must be after `mutationSpreadEnd`. |
| `secondMutationPosition` | Chromosome position for the second driver. |
| `secondSpreadStart` / `secondSpreadEnd` | Manual-spread window for the second driver. |
| `pop2GrowthStartTick` / `pop2GrowthEndTick` | Exponential-growth window for p2, at or after `populationSplitTick`. |

### Metastasis — migration parameters

| Parameter | Description |
| --- | --- |
| `initialPop2Size` | Starting size of the metastatic site p2. `0` = empty, populated only by migration (NWF standard). `>0` = pre-seeded from generation 1 (WF standard). |
| `migrationTick` *(NWF)* | Single generation at which a burst of cells moves from p1 to p2. |
| `migrantCount` *(NWF)* | Number of cells transferred at `migrationTick`. |
| `migrationStartTick` / `migrationEndTick` *(WF)* | Window during which the p1→p2 migration rate is active. Set equal for a single-generation pulse. |
| `migrationRate` *(WF)* | Fraction of p1 that migrates to p2 per generation while active (e.g. `0.01` = 1%/gen). |

## TreeOptions/Snapshot.txt — Tree Sampling Reference

Controls how cells are subsampled from a `.trees` file to build the Newick tree. Snapshot samples only cells alive at the final generation (like a biopsy) — good for comparing subclones and VAF-style analysis.

| Parameter | Description |
| --- | --- |
| `valid_pops` | Which populations to include. `all`, `1` (primary only), `2` (metastatic only), or `1,2` (both — use for metastasis models). |
| `snapshotSamples` | Max cells drawn from the final generation (all alive cells are used if fewer exist). |
| `seed` | Integer seed for reproducible sampling, or `none` for a different random sample each run. |

## CaseStudy — Predefined SLiM Scripts

`CaseStudy/` holds hand-written `.slim` scripts (not driven by the `SimOptions/` CSV/TXT system). Each script declares its own default `-d` constants in a first-line comment, e.g.:

```
//slim -d AGE=81 -d DRIVER_TICK=38 -d REP=1 Simulation4_ClonalHematopoiesis.slim
```

`python simtreemaker.py CaseStudy` runs **every** `.slim` file in the folder using those built-in defaults. **To change age, driver timing, or replicate, pass `KEY=VALUE` pairs on the command line** — they override the script's defaults for that run only (nothing is written back to the `.slim` file):

| Script | Overridable constants | What they control |
| --- | --- | --- |
| `CHIP2EnvN.slim` | `AGE` (default `100`), `REP` (default `1`) | `AGE` = generation the simulation stops at. `REP` = replicate label, used only in output filenames. |
| `Simulation4_ClonalHematopoiesis.slim` | `AGE` (default `81`), `DRIVER_TICK` (default `38`), `REP` (default `1`) | `AGE` = patient age / generation the simulation stops at. `DRIVER_TICK` = generation the driver mutation (`m2`, s=0.6) is introduced into a random newborn cell. `REP` = replicate label. Modeled after Mitchell et al. 2022 (*Nature*): default values represent an 81-year-old patient with the driver acquired at generation 38. |

```bash
python simtreemaker.py CaseStudy AGE=100 REP=1
python simtreemaker.py CaseStudy DRIVER_TICK=50 AGE=90 REP=3
```

Any constants you don't override keep the script's built-in default. Passing `KEY=VALUE` applies to *every* `.slim` script run that call — a script that doesn't declare that constant simply ignores it.

## Outputs

**Named-model runs** (`MutationSpread_NWF`, etc.) — output to `<ModelName>/Snapshot/`:

```
<ModelName>/Snapshot/
├── tree/      <treeOutputFile>.tree
├── newick/    <treeOutputFile>.nwk
├── pngTree/   <stem>_horizontal_labels.png / _no_labels.png
│              <stem>_vertical_labels.png   / _no_labels.png
└── log.txt    parameters, SLiM command, and SLiM output for this run
```

**CaseStudy** — output to `CaseStudyOutputs/<script>_<overrides>/Snapshot/`, e.g. `CaseStudyOutputs/Simulation4_ClonalHematopoiesis_AGE81_DRIVER_TICK38_REP1/Snapshot/`. The `<overrides>` suffix records every `-d` constant used for that run, so different `AGE`/`REP`/`DRIVER_TICK` combinations never overwrite each other:

```
CaseStudyOutputs/<script>_<overrides>/Snapshot/
├── newick/    <stem>.nwk
├── pngTree/   <stem>_horizontal_labels.png / _no_labels.png / _vertical_*.png
└── log.txt
```

The raw `.tree` file itself is written wherever the `.slim` script's `treeSeqOutput()` call points it (both current scripts write to the SimTreeMaker root or `CaseStudyOutputs/` directly) — the pipeline searches for it automatically after each SLiM run.

**Tree** (`ReadyTrees/`) — each dropped-in file is processed independently, output to `ReadyTreesOutputs/<stem>Output/`:

```
ReadyTreesOutputs/<stem>Output/
├── newick/    <stem>.nwk
├── pngTree/   <stem>_horizontal_labels.png / _no_labels.png / _vertical_*.png
└── log.txt
```

## Reference

If you use SimTreeMaker in your work, please cite:

Savannah L. Wilson¹, Prabin Dawadi¹*, Dikshya Niraula¹, Sayaka Miura¹*
¹Department of Biology, The University of Mississippi, University City, MS 38677, USA


## License

Copyright 2025, Authors and University of Mississippi

BSD 3-Clause "New" or "Revised" License

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
