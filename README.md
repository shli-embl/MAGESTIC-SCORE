# MAGESTIC-SCORE
## A genome-wide annotation of difficult-to-edit regions based on MAGESTIC and SCORE

### Install conda environment
Install conda environment using CONDA (replace {myenv_name} with the name of environment to create):
```
git clone https://github.com/shli-embl/MAGESTIC-SCORE.git
cd MAGESTIC-SCORE/
conda env create -f envs/environment.yaml -p envs/{myenv_name}
```
Or using MAMBA (recommended):
```
mamba env create -f envs/environment.yaml -p envs/{myenv_name}
```

Then, activate the environment before running the snakemake pipelines:
```
conda activate envs/{myenv_name}
```

### Snakemake pipelines
#### pipeline 1: Single clonal whole-genome sequencing data analysis and editing outcome (pre-)calling
#### pipeline 2: Target site sequence and genomic feature annotations
#### pipeline 3: Machine learning model construction and tuning
#### pipeline 4: Bulk whole-genome sequencing analysis for structral variant characterizations
