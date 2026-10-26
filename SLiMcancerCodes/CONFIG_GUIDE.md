# TG (TumorGrowth) configuration guide

`config.txt` holds the settings; this guide explains them. Format is
`KEY = VALUE   # comment`. Any key can be overridden for one run on the command
line without editing the file:

```
python simtreemaker.py                      # p1 only (same as nomigration)
python simtreemaker.py nomigration          # p1 only
python simtreemaker.py migration            # p1 + every PopN row
python simtreemaker.py nomigration noise    # ... then add sequencing noise
python simtreemaker.py migration noise
python simtreemaker.py nomigration treatment  # p1 with one drug dose (or set treatment = on)
python simtreemaker.py migration Pop2Time=30 targetSize=150
```

Outputs go to `<outputDir>/TG/NoMigration/` or `<outputDir>/TG/Migration/` (noise in
their `noise/` subfolder); each run is recorded in that folder's `log.txt`.
`outputDir` (Output section of config.txt) defaults to `../SLIMcanceroutputs`,
relative to this folder; leave it blank to write `TG/` next to the scripts.

---

## Run

| Key | Meaning |
|---|---|
| `simSeed` | Random seed for the whole run: SLiM's simulation **and** which cells are sampled. Same seed + same settings = identical tree, SNVs and noise. Leave empty for a random seed; the seed actually used is printed at the end ("Seed") and saved in `log.txt`, so any run can be repeated. |
| `simulationEndTick` | Number of generations simulated. There is one clock for the whole run: every population is saved and sampled together at this generation, and every `PopN` Time counts on the same clock. Must be later than every PopN Time. |

Longer runs give more mutations, more driver clones, larger populations and
slower runs. 200 is a reasonable default; use ~50 for quick tests.

## Primary population (p1)

| Key | Meaning |
|---|---|
| `initialPopSize` | Founder cells at generation 1. |
| `targetSize` | Carrying capacity for **driver-free** cells. |
| `Pop1SampleCells` | Cells sampled from p1 at the end for the tree, SNV CSV and FASTA. If p1 is smaller, all its cells are used. |

**Why p1 grows past `targetSize`:** each generation every cell makes one clone,
then each cell survives with probability `fitness x targetSize / N`. Driver-free
cells have fitness 1, so they settle at `targetSize`. Each driver multiplies
fitness by 1.15, which also raises that clone's effective carrying capacity. A
run with ~9 drivers per cell can end at roughly 10x `targetSize`.

## Mutation

| Key | Meaning |
|---|---|
| `mutationRate` | Neutral (m1) mutation rate per base, per genome copy, per generation. 0 = off. |
| `driverMutationRate` | Driver (m2) mutation rate, same units. 0 = no drivers. |
| `m2_effect` | Driver selection coefficient *s*. |
| `m2_dominance` | Driver dominance *h*. |

Cells are diploid and a new driver hits one copy, so each driver gives fitness
`1 + h*s` = 1 + 0.5 x 0.3 = **1.15** with the defaults. Drivers multiply. (s = 0.3 was chosen so that tumors carry a realistic ~7 drivers per cell; s = 0.6 gave ~12.)

Rates are per base, so the number of mutations scales with region length. The
run warns when fewer than one mutation is expected in a region.

## Migration (only with `migration`)

Each row adds one subpopulation. Columns, in order:

| Column | Meaning |
|---|---|
| Origin | Population it draws migrants from: `Pop1`, `Pop2`, ... (any population that exists by then). |
| Time | Generation it is created, on the shared clock. |
| IniPopSize | Founding migrants taken from Origin at Time. |
| TargetSize | Its carrying capacity for driver-free cells. |
| MigrationRate | Per-cell probability, each later generation, that an Origin cell migrates in. |
| SampleCells | Cells sampled from it at the end (all, if fewer). |

```
Pop2   = Pop1    50    5    100   0.01   100
```

- **Add a population:** add a `Pop5 = ...` row. No code changes.
- **Override one value for a run:** use `PopN<Column>=value`, e.g. `Pop2Time=30`.
- A population created at Time *t* grows on its own for `simulationEndTick - t`
  generations; its cells still carry the mutations they inherited before
  migrating, so all populations join into one tree.

## Treatment (no-migration runs only)

A single drug dose given to a single population (no migration). Runs
`TG_Treatment.slim`: the normal growth model plus resistance mutations (m3).
Output goes to `<outputDir>/TG/Treatment/`.

Turning it on and off:

| How | Effect |
|---|---|
| `treatment = on` in config.txt | Every no-migration run includes the dose (`python simtreemaker.py`) |
| `treatment = off` in config.txt (default) | No drug |
| `treatment` on the command line | Dose for this run, whatever config.txt says |
| `notreatment` on the command line | No drug for this run, whatever config.txt says |

Migration runs never get the dose: with `treatment = on` they print a note and
run untreated, and `migration treatment` on the command line is an error.
`manuscript.py` always runs untreated.

| Key | Meaning |
|---|---|
| `treatmentTargetSize` | Carrying capacity used instead of `targetSize` in treatment runs only (default 2500, 5x the normal 500). A larger tumor lets resistance start at a more realistic, lower fraction while still containing resistant cells; runs take about 20 s instead of 4 s. |
| `treatmentTick` | The one generation the drug is given (between 2 and `simulationEndTick` - 1). |
| `drugKill` | Chance that a sensitive cell dies at the dose (default 0.7, about a RECIST partial response: ≥30% smaller diameter ≈ 65–70% smaller volume). Resistant cells are unaffected. |
| `resistanceCost` | Fitness cost of resistance at every generation: m3 has selection coefficient `-resistanceCost` with dominance 1, so 0.05 means resistant cells are 5% less fit. |
| `resistanceMutationRate` | m3 mutation rate per base per generation. |
| `m3Start`, `m3End` | Resistance region (same rules as the other regions; must not overlap them). |

How it behaves:
- A cell is resistant if it carries at least one m3. Resistance arises by
  chance throughout growth, so resistant cells may or may not exist when the
  dose is given; whether they do differs between seeds.
- Cells can double every generation, so a 90% kill is regrown in about 4
  generations; the lasting effect is on the resistant share, not on size.
- Afterwards resistance can fade (its cost, and competition with driver clones)
  or take over; the outcome varies between runs.
- With `drugKill = 1` and no resistant cell, the tumor is eradicated: the run
  stops, `log.txt` records `ERADICATED`, and no tree is written.

Outputs, in addition to the usual tree/SNV/FASTA files (m3 appears in tip
labels, the SNV CSV and the positions file):

| File | Contents |
|---|---|
| `trajectory/..._trajectory.csv` | Per generation: total, sensitive and resistant cells, and whether the dose was given. |
| `trajectory/..._trajectory.png` | Top: cells per generation; bottom: resistant share; dose marked. |

The end-of-run summary reports the size before and after the dose, how many
generations it took to regrow, and the resistant share before, right after,
and at the end.

## Genome layout

| Key | Meaning |
|---|---|
| `m1Start`, `m1End` | Neutral (m1) region. |
| `m2Start`, `m2End` | Driver (m2) region. |

Rules: whole numbers >= 0, 0-based, both ends inclusive (0-49999 is 50,000
bases), start <= end, regions must not overlap. They may be in either order
and gaps are allowed; positions outside both regions never mutate. Genome
length = largest end + 1. Invalid layouts stop the run with an error before
SLiM starts.

## Output

| Key | Meaning |
|---|---|
| `treeOutputFile` | Base name for `.trees`, `.nwk`, `_snv.csv`, FASTA and PNG files. |
| `fastaSeed` | Seed for the parent genome and the mutated bases. The model tracks *where* mutations are, not which base they change to, so bases are assigned afterwards from this seed. |
| `writeFullFasta` | `T` also writes `_full.fasta`: the full genome for every sampled cell (~100 KB per cell with the default layout). |

### FASTA files

| File | Length per cell | Contents |
|---|---|---|
| `_parent.fasta` | Full genome | The founder's starting genome. |
| `_snv.fasta` | One letter per SNV site | Each cell's base at only the mutated positions, so every sequence in a run has the same length (= number of SNV sites). |
| `_snv_positions.tsv` | - | Which genome position, mutation type, REF and ALT each `_snv.fasta` column stands for. |
| `_full.fasta` | Full genome | Only with `writeFullFasta = T`. |

---

## Sequencing noise (on / off)

Turns the clean `_snv.fasta` into what a single-cell DNA experiment would
observe. Runs as the last step of both modes when it is switched on:

| How | Effect |
|---|---|
| `noise = on` in config.txt | Every run adds noise |
| `noise = off` in config.txt (default) | No noise; the rate settings below are ignored |
| `noise` on the command line | Noise for this run, whatever config.txt says |
| `nonoise` on the command line | No noise for this run, whatever config.txt says |

`manuscript.py` always simulates with `nonoise` and applies its own noise
sweeps afterwards.

| Key | Applies to | Effect |
|---|---|---|
| `missingRate` | every site of every cell | **Average** chance a site becomes `?` (no coverage). |
| `cellSpread` | cells | 0 = every cell has `missingRate`. Higher = some cells much better or worse (each cell's rate is drawn from a Beta distribution with mean `missingRate`). Must be < 1. |
| `siteSpread` | sites | Same idea per site: some positions drop out far more often. Must be < 1. |
| `fpRate` | observed sites where the cell has the normal base | Replaced by that site's ALT base (amplification/sequencing error). |
| `fnRate` | observed sites where the cell is mutated | Replaced by the normal base (allelic dropout). |
| `noiseSeed` | - | Same seed + same run = identical noise. Leave empty for different noise each run. |

FP and FN apply only to sites that are not missing. Typical published ranges:
missing 0.1-0.5, FP 0.00001-0.01, FN 0.1-0.3. Real cells vary a lot in
quality, so a `cellSpread` of ~0.3-0.6 is more realistic than 0.

The SNV FASTA has no normal cell, so one named `Normal` is built from the REF
column of `_snv_positions.tsv` (the parent's base at each SNV site) and
written last, unchanged.

Outputs in `<outputDir>/TG/<mode>/noise/` (cleared at the start of every run, since old
noise would describe the previous run's cells), in two folders made
automatically:

| Folder | Contents |
|---|---|
| `control/` | Same step with every rate 0: the clean cells + `Normal`, file name `..._Miss_0_FP_0_FN_0`. The run checks it matches the clean SNV FASTA exactly and stops with an error if not. |
| `test/` | The same cells with the `config.txt` rates applied. |

Run any downstream tool (tree inference, clustering) on both: control is the
best that tool can do with these SNVs, so the gap between control and test is
the effect of the noise alone. Each folder holds:

| File | Contents |
|---|---|
| `..._snv_Miss_<m>_FP_<fp>_FN_<fn>.fasta` | Noisy cells + `Normal`. |
| `..._log.tsv` | Every changed cell/site: MISS, FP or FN, true and observed base. |
| `..._cells.tsv` | Each cell's target and observed missing rate. |

### Running it on any FASTA

`add_noise.py` also works on its own, e.g. on an older run or a FASTA from
elsewhere (output goes to a `noise/` folder next to the input):

```
python add_noise.py <fasta> <missing> <FP> <FN> <normal name> [--seed N] [--cell-spread S] [--site-spread S] [--out-dir DIR]
python add_noise.py Test.fasta 0.5 0.5 0.5 Nor
```
