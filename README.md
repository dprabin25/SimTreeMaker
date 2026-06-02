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
  simtreemaker.py
  slim_newick.py
  slim_config.txt       ## edit SLIM_EXE path 
  Options/              ## Contains paramters for simulation
    MutationSpread.csv  ### Add/Edit paramters for Mutation spread
    ClonalGrowth.csv    ### Add/Edit paramters for ClonalGrowth
    Metastasis.csv      ### Add/Edit paramters for Metastasis
  CaseStudy/            ## Case study
   *.slim               #### Case study slim codes
   *.py                 #### 
  ReadyTrees/           ## drop .trees files here if you want to generate output for 

