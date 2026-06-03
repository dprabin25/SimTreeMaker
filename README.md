# SimTreeMaker

A pipeline for running SLiM cancer evolution simulations and converting tree sequences into Newick phylogenies and visualizations.

## Dependencies

### Anaconda

Install Anaconda: https://www.anaconda.com/distribution/

Then create the environment:

```bash
conda create -n SimTreeMaker -c conda-forge python=3.11 biopython=1.85 matplotlib pyslim=1.1.1 tskit=1.0.3 -y
conda activate SimTreeMaker
```

## File Structure

```
SimTreeMaker/
├── simtreemaker.py       # Main pipeline script
├── slim_newick.py        # Converts SLiM tree sequences to Newick format
├── slim_config.txt       # SLiM executable path — only edit SLIM_EXE if path changes
├── Options/
│   ├── MutationSpread.csv
│   ├── ClonalGrowth.csv
│   └── Metastasis.csv
├── CaseStudy/
│   ├── CHIP2EnvN.slim    # Forward-time SLiM simulation script
│   └── Slim2tree3.py     # Post-processing / tree conversion for case study
└── ReadyTrees/           # Drop .trees files here for direct processing
```

## Core Files

| File                | Description                                       |
| ------------------- | ------------------------------------------------- |
| `simtreemaker.py`   | Main pipeline script                              |
| `slim_newick.py`    | Converts SLiM tree sequences to Newick format     |
| `slim_config.txt`   | SLiM executable path (edit `SLIM_EXE` if needed) |

## Input Files (User-Editable)

| Location                     | Description                                    |
| ---------------------------- | ---------------------------------------------- |
| `Options/MutationSpread.csv` | Mutation dynamics simulation parameters        |
| `Options/ClonalGrowth.csv`   | Clonal expansion simulation parameters         |
| `Options/Metastasis.csv`     | Metastasis simulation parameters               |
| `CaseStudy/CHIP2EnvN.slim`   | Predefined forward-time SLiM simulation        |
| `CaseStudy/Slim2tree3.py`    | Post-processing and tree conversion            |
| `ReadyTrees/*.trees`         | Pre-existing tree sequence files               |

## User Workflow Options

### Option 1: Parameter-driven simulation

Use the CSV files in `Options/` to define simulation behavior, then run using the CSV filename (without `.csv`):

```bash
python simtreemaker.py MutationSpread
python simtreemaker.py ClonalGrowth
python simtreemaker.py Metastasis
```

---

### Option 2: Case study execution

Directly runs SLiM scripts in `CaseStudy/`. Uses `Slim2tree3.py` for tree conversion when present, otherwise falls back to `slim_newick.py`.

```bash
python simtreemaker.py CaseStudy
```

---

### Option 3: Tree sequence processing

Place `.trees` files in `ReadyTrees/` and run:

```bash
python simtreemaker.py Tree
```

---

## Outputs

For each model row processed, outputs are organized under `<ModelName>/<stem>/`:

```
<ModelName>/<stem>/
├── scripts/   <stem>.slim                         ← generated SLiM script
├── tree/      <stem>.tree                         ← raw SLiM tree sequence
├── newick/    <stem>.nwk                          ← Newick phylogeny
└── pngTree/   <stem>_horizontal_labels.png
               <stem>_horizontal_no_labels.png
               <stem>_vertical_labels.png
               <stem>_vertical_no_labels.png
```

For **CaseStudy**, outputs go to `CaseStudyOutputs/CaseStudyTrees/`.  
For **ReadyTrees**, outputs go to `ReadyTreesOutputs/<stem>Output/`.

## Command Reference
<img width="647" height="116" alt="image" src="https://github.com/user-attachments/assets/c1928fc5-5cd4-4161-ae69-400227c41019" />
Make sure you are within working dir of SimTreeMaker. You can then type the following command for your prefered action. 

| Command                                   | Action                                          |
| ----------------------------------------- | ----------------------------------------------- |
| `python simtreemaker.py MutationSpread`   | Runs mutation spread simulation                 |
| `python simtreemaker.py ClonalGrowth`     | Runs clonal growth simulation                   |
| `python simtreemaker.py Metastasis`       | Runs metastasis simulation                      |
| `python simtreemaker.py CaseStudy`        | Runs predefined SLiM case study scripts         |
| `python simtreemaker.py Tree`             | Processes `.trees` files from `ReadyTrees/`     |

## Reference

If you use SimTreeMaker in your work, please cite:

> Prabin Dawadi,...............(2026) SLiM cancer evolution simulations and converting tree sequences into Newick phylogenies and visualizations. Under Review.

## License

Copyright 2025, Authors and University of Mississippi

BSD 3-Clause "New" or "Revised" License

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

