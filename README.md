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
#### Pipeline 1: Single clonal whole-genome sequencing data analysis and editing outcome (pre-)calling
The single clonal WGS analysis pipeline takes in .fastq read files from multiple Cas9 edited strains isolated from cell pool created by MAGESTIC, so that each single genome contains a barcode linked to a pre-designed gRNA-donor-DNA pair to introduce a variant on the target locus. The pipeline contains read processing and on-target genotyping (both for SNV/indel detection via GATK4 and SV detection via a 3-in-1 pre-selection approach). The SV pre-scan includes GATK4 (to detect unintended SNV/indels proximate to the target), SvABA and CNV detection. We also extended the outcome analysis to a list of off-target loci predicted by edit-distance from the on-target sequence. 
#### 1a. Modify the configuration files
```
cd 1_outcome_precalling
ls inputs/
```
file 1. gRNA_list.txt (sample id and gRNA sequence): 
```
MAGESTIC_REDI_1	CTCCTCAGTTTGTCCACGGT
MAGESTIC_REDI_10	AGAGGAAGTCTCAACGGCAG
MAGESTIC_REDI_100	ATACTGGCCACGTTTGACAA
MAGESTIC_REDI_1000	TTGTGATTTTATTGATTCTG
MAGESTIC_REDI_1001	GGTAACAAAGTCACGGCTCC
...
```
file 2. input.yaml (required by snakemake for locating file paths
```
---
samplesheet: 'inputs/samples.txt'
grna: 'inputs/gRNA_list.txt'
genome: 'inputs/fasta/yeast_reference.fa'
tmp_dir: './tmpdir'
output_dir: './output'
offtar: 'inputs/off_target.txt'
```
file 3. samples.txt (sample information, read file paths, with header)
```
sample	read1fq	read2fq	chr	strand	PAM_coordinate	variant_coordinate	template_start_coordinat
e	template_end_coordinate	reference_allele	mutation_allele	template_sequence	target_m
utations	synthetic_errors
MAGESTIC_REDI_1	/g/steinmetz/incoming/solexa/2018-10-24-HCW3HAFXY/HCW3HAFXY_CST85_18s005308-1-1_Szu-Tu_l
ane1012AA1_1_sequence.txt.gz	/g/steinmetz/incoming/solexa/2018-10-24-HCW3HAFXY/HCW3HAFXY_CST85_18s005
308-1-1_Szu-Tu_lane1012AA1_2_sequence.txt.gz	chr8	-	104423	104423	104371	104479	CCCA	
ACCC	TGGCTGCTATCGTTGAAATTATCGACCAAAAGAAGGTATGTTGAACCTAAAAACCCCCGTGGACAAACTGAGGAGGAAATTGTAAGGAAGAGAAAG
TCCCCGTATGTTC	chr8_104423_C_to_A,chr8_104426_A_to_C	
...
```
Note: template_sequence and synthetic_errors were obtained from REDI sequencing of the barcode locus prior to the WGS. <br><br>
file 4. off_target.txt (list of off-target location predicted from edit distance <=3)
```
sample	guide	chr	PAM_coordinate	match	strand	mismatch
MAGESTIC_REDI_10	AGAGGAAGTCTCAACGGCAGNNN	chr4	1503943	AGAGGAAGTCTCAACaGCAGAGG	+	1
MAGESTIC_REDI_10	AGAGGAAGTCTCAACGGCAGNNN	chr4	1503961	AGAGGAAGTCTCAACaGCAGAGG	+	1
MAGESTIC_REDI_10	AGAGGAAGTCTCAACGGCAGNNN	chr4	1503889	AGAGGAAGTCTCAgCaGCAGAGG	+	2
MAGESTIC_REDI_10	AGAGGAAGTCTCAACGGCAGNNN	chr4	1503925	AGAaGAAGTCTCAACaGCAGAGG	+	2
MAGESTIC_REDI_10	AGAGGAAGTCTCAACGGCAGNNN	chr4	1503979	AGAGGAAGTCTCAgCaGCAGAGG	+	2
MAGESTIC_REDI_10	AGAGGAAGTCTCAACGGCAGNNN	chr4	1503853	AGAaGAAaTCTCAACaGCAGAGG	+	3
...
```
<br>
#### 1b. Run pipeline
We use a HPC cluster for computing the task, for which the configuration of job submission and control are saved in profile/config.yaml. 
```
cd 1_outcome_precalling
snakemake --profile profile
```


#### Pipeline 2: Target site sequence and genomic feature annotations
#### Pipeline 3: Machine learning model construction and tuning
#### Pipeline 4: Bulk whole-genome sequencing analysis for structral variant characterizations
