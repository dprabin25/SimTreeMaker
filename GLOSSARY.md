# Glossary

Terms used across `README.md`, organized the same way as that document — Method 1 defines the full shared vocabulary; each later method only adds the terms it introduces on top of that, rather than repeating everything. Each method section is split into general concepts/files, and the actual Eidos/SLiM code terms (functions, properties, callbacks) used literally inside that method's `.slim` file.

---

## Method 1 — Tumor growth, no migration

### Concepts & project files

**SLiM** — the simulation program (installed as `slim` / `slim.exe`) that runs `.slim` model files and produces `.trees` output.

**Eidos** — the scripting language `.slim` files are written in (similar in spirit to R or Python, but SLiM-specific).

**`.slim` file** — a text file containing an Eidos script defining a simulation model (e.g. `TG.slim`). This is what SLiM actually runs.

**nonWF (nonWF model)** — "non Wright-Fisher": a SLiM model type where you control births, deaths, and population size yourself, instead of SLiM enforcing fixed, non-overlapping generations automatically. `TG.slim` is a nonWF model.

**`-d` flag** — a command-line flag (e.g. `-d initialPopSize=1`) that passes a constant value into a `.slim` script from the terminal, without editing the script itself.

**Subpopulation (`p1`, `p2`, ...)** — a named group of cells within the simulation. `p1` is always the founding population.

**Founder cell** — the single starting cell (or cells) a simulation begins from at generation 1, before any growth or mutation has happened.

**Carrying capacity / `targetSize`** — the population size density-dependent selection pulls toward (not a hard limit — actual size fluctuates around it).

**Density-dependent selection** — a regulation method where fitness/survival depends on how crowded the population currently is, rather than killing individuals at a fixed age or count.

**Mutation type (`m1` / `m2`)** — a category of mutation defined in a `.slim` file. Here, `m1` = neutral/background mutation, `m2` = driver mutation (confers a fitness advantage).

**Driver mutation** — a mutation that gives the cell carrying it a selective/fitness advantage over other cells, letting its lineage (clone) expand faster.

**Passenger / background mutation** — a neutral mutation that doesn't affect fitness; it just accumulates over time and is useful for tracking lineage history.

**Selection coefficient (`s`)** — how much fitness advantage (or disadvantage) a mutation confers. Written as `m2_effect` in `config.txt` (e.g. `0.6` = 60% more fit).

**Dominance** — how much of a mutation's fitness effect shows up when only one copy is present (vs. two). Written as `m2_dominance` in `config.txt`.

**Mutation rate** — the probability of a new mutation occurring per base pair per generation (e.g. `mutationRate = 1e-8`).

**Genomic element / genomic position** — a defined stretch of positions in the simulated genome (e.g. positions 0–49999) that a particular mutation type is allowed to occur in.

**Recombination rate** — the probability that a genome breaks and recombines between two given positions during reproduction. Set to `0` in every model in this project (no recombination — each genome passes down as one intact block).

**Tree sequence** — the internal data structure SLiM's `treeSeqOutput()` saves to a `.trees` file: the full ancestry (who's descended from whom) of every cell in the simulation, plus their mutations.

**`.trees` file** — the raw binary output SLiM writes via `treeSeqOutput()`. Contains the full tree sequence with mutation metadata. Not human-readable directly — needs conversion.

**tskit** — the Python library used to read/process `.trees` files (tree sequence toolkit). Used inside `slim_utilities.py`, not something you run directly.

**Newick format / `.nwk` file** — a standard, simple text format for representing a tree's branching structure. This project's readable tree output.

**Tip label** — the text label attached to each leaf (sampled cell) in a Newick tree, formatted here as `gen0_ind<i>|pop<N>|<mutations>`.

**Branch length** — the distance (in generations, here) between a node and its parent in a tree — represents how much time/evolution separates them.

**Root / leaf / node** — tree vocabulary: the root is the tree's starting point (common ancestor); leaves (or tips) are the sampled cells at the end; nodes are any branching point in between.

**SNV (single nucleotide variant)** — a single-position mutation/variant in the genome.

**SNV CSV** — a spreadsheet-style table (one row per cell/mutation pair) listing `cell_id`, `population`, `mutation_type`, `position`, `generations_before_sampling`, `origin_tick`.

**Tree PNG** — a rendered image of the Newick tree, with (`labels`) or without (`no_labels`) tip text, via R+`ape` or a matplotlib fallback.

**Snapshot sampling** — the sampling strategy used throughout this project: only cells alive at the very last generation of the run are eligible to be sampled for the tree/SNV output, like a biopsy at the end.

**`snapshotSamples`** — the max number of cells to sample (pooled across populations) for the tree/SNV output. Fixed at `400` directly in `slim_utilities.py`'s `load_tree_options()` — not read from a config file.

**`valid_pops`** — which population(s) to sample from ("all", or a specific list of population numbers). Fixed at `"all"` directly in `slim_utilities.py`'s `load_tree_options()` — not read from a config file.

**`stem`** — the base filename (without extension) used to name a run's set of outputs consistently — e.g. `TumorGrowth` — so the `.trees`, `.nwk`, `_snv.csv`, and `.png` files from the same run are easy to match up.

**`TG.slim`** — the Method 1/2 model file, the "backbone" tumor growth model, with no migration code in it.

**`simtreemaker.py`** — the Python runner for Method 1 and Method 2.

**`config.txt`** (TumorGrowth) — holds the model parameters for Method 1/2 (mutation rates, population sizes, migration settings, etc.).

**`slim_config.txt`** — holds one setting, `SLIM_EXE`: the file path to your installed `slim` executable.

**`log.txt`** — a running record appended to by every method's runner: one entry per run, recording the exact SLiM command used, SLiM's console output, and the resulting output paths (or the failure reason).

### Eidos/SLiM code terms used in `TG.slim`

**`initialize()` callback** — the block that runs once, before the simulation starts, to set up the model's static configuration (mutation types, genomic structure, model type, etc.).

**`early()` callback** — a block of code that runs at the start of every tick, before reproduction happens. Used here for density-dependent fitness scaling.

**`late()` callback** — a block of code that runs at the end of every tick, after reproduction and fitness evaluation. Used here to check whether the simulation should stop and write output.

**`reproduction()` callback** — a block of code that runs once per tick per individual, to generate that individual's offspring.

**Tick range prefix (e.g. `1:1000000`)** — the number(s) written before a callback block, restricting which ticks it runs during. `1:1000000` means "every tick from 1 to 1,000,000" — a generous upper bound so the block keeps running until the model's own stopping condition (`simulationEndTick`) is reached; it isn't the model's real intended length.

**`c()`** — Eidos's function for combining values into a vector (list), e.g. `c(mutationRate, driverMutationRate)`.

**`catn()`** — prints a line of text to the console, with a newline at the end. Used here for status messages like "SLiM done."

**`initializeSLiMModelType()`** — declares whether the model is `"WF"` (Wright-Fisher, SLiM's default) or `"nonWF"`. Must be called first inside `initialize()`.

**`initializeTreeSeq()`** — turns on tree-sequence recording, so the simulation's full ancestry can later be written out with `treeSeqOutput()`.

**`initializeMutationRate()`** — sets how often new mutations occur — either one flat rate for the whole genome, or, as in `TG.slim`, different rates for different regions via paired vectors.

**`initializeMutationType()`** — defines a named category of mutation (like `m1` or `m2`): its dominance coefficient and its fitness-effect distribution. `"f"` here means a fixed/constant effect size, not drawn from a random distribution.

**`convertToSubstitution`** — a mutation-type property, set to `F` here so that even a mutation that eventually spreads to every individual (fixes) is still kept in the record individually, rather than being collapsed into a simpler "substitution" entry.

**`initializeGenomicElementType()`** — defines a named type of genome region (like `g1` or `g2`) and which mutation type(s) can occur in it.

**`initializeGenomicElement()`** — places a genomic element type at a specific range of positions in the simulated genome (e.g. `g1` covering positions 0–49999).

**`initializeRecombinationRate()`** — sets the recombination rate (see "Recombination rate" above); `0` in every model here.

**`sim.addSubpop()`** — creates a new subpopulation with a given ID and starting size.

**`sim.cycle`** — the current tick/generation number, as the simulation runs.

**`sim.treeSeqOutput()`** — writes the full recorded tree sequence out to a `.trees` file at the given path.

**`sim.simulationFinished()`** — tells SLiM to stop running after the current tick completes.

**`subpop.addCloned()`** — creates one cloned (asexual, genetically identical) offspring from a given parent individual.

**`individualCount`** — a subpopulation property giving its current number of living individuals.

**`fitnessScaling`** — a per-individual or per-subpopulation property SLiM uses to scale survival/reproduction likelihood; set here each generation to `targetSize / currentPopulationSize`.

---

## Method 2 — Tumor growth, with migration

Everything above still applies — Method 2 is Method 1 plus these additional terms only.

### Concepts & config

**Migration (in this project)** — cells moving from one subpopulation into another during the simulation (e.g. from `p1` into `p2`), simulating metastasis or spatial spread.

**`migrationEnabled`** — a `config.txt` value (`T` or `F`) that `simtreemaker.py` reads to decide whether to build and run the migration variant of `TG.slim`. Never passed to SLiM itself — Python-only.

**`PopN` (`PopNTime`, `PopNOrigin`, etc.)** — the set of five `config.txt` keys describing one migration-seeded population `N`: its creation generation, starting size, capacity, migration rate, sample count, and source population. `N` can be any number ≥ 2, not fixed to 2/3/4.

**Marker comment** — the inert 3-line comment inside `TG.slim`'s `early()` block (`//Add migration setting if any//`) that `simtreemaker.py` finds and replaces with generated migration code when building the migration variant.

**Temporary migration script** — the short-lived `.slim` file `simtreemaker.py` writes the generated migration code to before running SLiM against it, then deletes immediately afterward. Nothing named `TG_Migration.slim` is saved permanently.

### Eidos/SLiM code terms used in the generated migration code

(These appear only in the temporary migration `.slim` file `simtreemaker.py` generates — not in `TG.slim` itself.)

**`sim.subpopulations`** — a list of every subpopulation that currently exists in the simulation.

**`id` (subpopulation property)** — a subpopulation's numeric ID (e.g. `1` for `p1`), used here to look up a specific source population by number.

**`size()`** — returns how many elements are in a vector/list — e.g. how many individuals are in a filtered group.

**`min()`** — returns the smaller of two (or more) values — used here to cap how many migrants are pulled in at once.

**`sample()`** — randomly draws a given number of elements from a vector without replacement — e.g. picking specific random migrant individuals.

**`individuals` (subpopulation property)** — the list of every living individual currently in a subpopulation.

**`takeMigrants()`** — moves a given set of individuals into a subpopulation from wherever they currently are.

**`rbinom()`** — draws a random number from a binomial distribution — used here to decide how many individuals migrate in a given generation, from a population size and a per-individual probability.

---

## Method 3 — Process existing trees (no SLiM run)

Reuses the file-format and output terms from Method 1 (`.trees`, Newick, SNV CSV, Tree PNG, `stem`, `log.txt`). No `.slim` file runs in this method, so no new Eidos/SLiM code terms — the only new term is the tool itself:

**`process_trees.py`** — the Python runner for Method 3: batch-converts any `.trees` files already sitting in the `Trees\` folder into Newick/SNV/PNG, without running SLiM at all.

---

## Method 4 — Case study (clonal hematopoiesis)

Reuses the general Eidos/SLiM terms from Method 1 (`initialize()`, `early()`/`late()`, `reproduction()`, `initializeSLiMModelType()`, `initializeTreeSeq()`, `initializeMutationRate()`, `initializeMutationType()`, `initializeGenomicElementType()`, `initializeGenomicElement()`, `initializeRecombinationRate()`, `sim.addSubpop()`, `sim.cycle`, `sim.treeSeqOutput()`, `sim.simulationFinished()`, `subpop.addCloned()`, `individualCount`, `fitnessScaling`, `catn()`). New terms specific to this case study:

### Concepts & project files

**Clonal hematopoiesis** — a real biological phenomenon where a single blood stem cell's descendants (a "clone") come to make up an unusually large share of a person's blood cells, usually due to a driver mutation acquired with age. This is the condition this case study's model represents (modeled after Mitchell et al. 2022, *Nature*).

**Clone / subclone** — a group of cells all descended from one common ancestor cell that acquired a particular mutation — they share that mutation and everything else that ancestor had.

**`Simulation4_ClonalHematopoiesis.slim`** — the Method 4 model file (not related to `TG.slim`).

**`casestudy.py`** — the Python runner for Method 4.

**`AGE`** — the generation the simulation stops at, representing the patient's age at sampling.

**`DRIVER_TICK`** — the generation a driver mutation is introduced into a random newborn cell.

**`REP`** — a replicate label used only to keep output filenames from repeated runs apart. Does not seed SLiM's randomness, so it doesn't guarantee reproducible results by itself.

### Eidos/SLiM code terms used in `Simulation4_ClonalHematopoiesis.slim`

**`age` (individual property)** — how many ticks/generations an individual has been alive. Used here to remove cells 35 generations or older, approximating stem-cell attrition.

**`for` loop (`for (x in y) {...}`)** — Eidos's loop syntax for running the same code once per element in a list — used here to check every individual's age each generation.

**Boolean vector indexing (e.g. `inds[inds.age == 0]`)** — filtering a list down to just the elements where a condition is true — here, picking out only the newborn cells (age 0) as candidates for the driver mutation.

**`haplosomes` (individual property) / haplosome** — an individual's copies of its genome; indexing into it (`ind.haplosomes[hapIndex]`) picks one specific copy to place a mutation on.

**`addNewDrawnMutation()`** — manually places a specific mutation (of a given type) at a specific genome position on a given haplosome — used here to deliberately introduce the driver mutation at `DRIVER_TICK`, rather than waiting for it to arise on its own via the background mutation rate.

**Range literal (e.g. `0:1`)** — shorthand for a vector of consecutive integers — here, `0:1` is the vector `[0, 1]`, used to randomly pick one of an individual's two haplosome copies.
