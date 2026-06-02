# SimTreeMaker

## Dependencies
1. Anaconda
- Open Anaconda terminal and then create conda environment for SimTreeMaker.

Please install Anaconda : https://www.anaconda.com/distribution/

To create the environment, include to install correct Python version and R.

```bash
conda create -n Sample_Tree python=3.12.2 r-base=4.3.0 pandas -y

conda activate Sample_Tree
```

#### Python packages

Installing pandas, numpy and graphviz packages
```
conda install -c conda-forge pandas=2.2.2 numpy=1.26.4
conda install -c conda-forge graphviz python-graphviz
```

Installing opanai package
```
conda install conda-forge::openai==2.30.0   
```
