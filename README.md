# MAGESTIC-SCORE
## A genome-wide annotation of difficult-to-edit regions based on MAGESTIC and SCORE

### 1. Install conda environment
a. Install conda environment using CONDA (replace {myenv_name} with the name of environment to create):
```
git clone https://github.com/shli-embl/MAGESTIC-SCORE.git
cd MAGESTIC-SCORE/
conda env create -f envs/environment.yaml -p envs/{myenv_name}
```
Or using MAMBA (recommended):
```
git clone https://github.com/shli-embl/MAGESTIC-SCORE.git
cd MAGESTIC-SCORE/
mamba env create -f envs/environment.yaml -p envs/{myenv_name}
```

b. Then, activate the environment before running the snakemake pipelines:
```
conda activate envs/{myenv_name}
```

### 2. Snakemake pipelines
#### pipeline 1: Single clonal whole-genome sequencing data analysis and editing outcome (pre-)calling
The single WGS analysis pipeline takes in .fastq read files from multiple Cas9 edited strains isolated from MAGESTIC pool, so that each single genome contains a barcode linked to a pre-designed gRNA-donor-DNA pair to introduce a variant on the target locus. The pipeline contains read processing and on-target genotyping (both for SNV/indel detection via GATK4 and SV detection via a 3-in-1 pre-selection approach). The SV pre-scan includes GATK4 (to detect unintended SNV/indels proximate to the target), SvABA and CNV detection. We also extended the outcome analysis to a list of off-target loci predicted by edit-distance from the on-target sequence. 
1a. 
#### pipeline 2: Target site sequence and genomic feature annotations
#### pipeline 3: Machine learning model construction and tuning
#### pipeline 4: Bulk whole-genome sequencing analysis for structral variant characterizations
