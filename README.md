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
Inside "SimTreeMaker" dir there are these files. 

# Core Files
| File              | Type   | Description                                   |
| ----------------- | ------ | --------------------------------------------- |
| `simtreemaker.py` | Script | Main pipeline script                          |
| `slim_newick.py`  | Script | Converts SLiM tree-sequences to Newick format |
| `slim_config.txt` | Config | SLiM executable path + working directory      |


# Input Files
(User-defined / Editable)
| Location                     | Type          | Description                                |
| ---------------------------- | ------------- | ------------------------------------------ |
| `Options/`                   | Folder        | Simulation parameter files                 |
| `Options/MutationSpread.csv` | CSV           | Mutation dynamics parameters               |
| `Options/ClonalGrowth.csv`   | CSV           | Clonal expansion parameters                |
| `Options/Metastasis.csv`     | CSV           | Metastasis simulation parameters           |
| `CaseStudy/*.slim`           | SLiM script   | Forward-time simulation models             |
| `CaseStudy/*.py`             | Python script | Post-processing / analysis of SLiM outputs |
| `ReadyTrees/*.trees`         | Tree file     | SLiM tree-sequence input files             |

# User Workflow Options
- Option 1: Parameter-driven simulation

Use files in Options/

Run:

Clonal growth simulation

Mutation spread simulation

Metastasis simulation


- Option 2: Case study execution

Directly run SLiM scripts in CaseStudy/

Includes predefined simulation + analysis workflows


- Option 3: Tree sequence processing

Place .trees files in ReadyTrees/

Run pipeline to generate:
Newick (.nwk) trees
PNG visualizations
