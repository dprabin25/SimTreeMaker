# SimTreeMaker

## Dependencies
### Anaconda
- Open Anaconda terminal and then create conda environment for SimTreeMaker.

Please install Anaconda : https://www.anaconda.com/distribution/

To create the environment, include to install correct Python version and R.

```bash
conda create -n SimTreeMaker -c conda-forge python=3.11 biopython=1.85 matplotlib pyslim=1.1.1 tskit=1.0.3 -y

conda activate SimTreeMaker
```

## Files required
SimTreeMaker/
├── simtreemaker.py        # Main pipeline script
├── slim_newick.py         # SLiM tree-sequence → Newick conversion
├── slim_config.txt        # SLiM executable + working directory config
│
├── Options/               # Simulation parameters (edit here)
│   ├── MutationSpread.csv
│   ├── ClonalGrowth.csv
│   └── Metastasis.csv
│
├── CaseStudy/
│   ├── *.slim             # SLiM simulation scripts
│   └── *.py               # Post-processing scripts for SLiM outputs
│
└── ReadyTrees/
    └── *.trees            # SLiM tree-sequence input files
