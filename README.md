# Anneal-Path

Simulated Annealing-based branched novel biochemical pathway prediction tool

## Quick Installation of Environment

> Create a conda environment for RDkit, `conda create -c rdkit -n rdkit3 rdkit python=3.7`.

> For using RDKit in Jupyter Notebooks, `conda install -n rdkit3 nb_conda_kernels`.

> Activating the conda environment, `conda activate rdkit3` and install other required packages using `pip` or `conda`,
```
pip install pymysql
pip install tqdm
```

## Change Path to R in AP_output.py 

Change Path to R in AP_output.py line301 for generating graphical output of the program, 

e.g., `R_path="/bin/R"` or `"C:\\R\\R-4.1.1\\bin\\Rscript.exe"`

