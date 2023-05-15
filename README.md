# MAGESTIC-SCORE
## A genome-wide annotation of difficult-to-edit regions based on MAGESTIC and SCORE

### Install conda environment
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

### Snakemake pipelines
### Pipeline 1: Single clonal whole-genome sequencing data analysis and editing outcome (pre-)calling
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
file 2. input.yaml (required by snakemake for locating file paths):
```
---
samplesheet: 'inputs/samples.txt'
grna: 'inputs/gRNA_list.txt'
genome: 'inputs/fasta/yeast_reference.fa'
tmp_dir: './tmpdir'
output_dir: './output'
offtar: 'inputs/off_target.txt'
```
file 3. samples.txt (sample information, read file paths, with header):
```
sample	read1fq	read2fq	chr	strand	PAM_coordinate	variant_coordinate	template_start_coordinate	template_end_coordinate	reference_allele	mutation_allele	template_sequence	target_mutations	synthetic_errors
MAGESTIC_REDI_1	/g/steinmetz/incoming/solexa/2018-10-24-HCW3HAFXY/HCW3HAFXY_CST85_18s005308-1-1_Szu-Tu_lane1012AA1_1_sequence.txt.gz	/g/steinmetz/incoming/solexa/2018-10-24-HCW3HAFXY/HCW3HAFXY_CST85_18s005308-1-1_Szu-Tu_lane1012AA1_2_sequence.txt.gz	chr8	-	104423	104423	104371	104479	CCCA	ACCC	TGGCTGCTATCGTTGAAATTATCGACCAAAAGAAGGTATGTTGAACCTAAAAACCCCCGTGGACAAACTGAGGAGGAAATTGTAAGGAAGAGAAAGTCCCCGTATGTTC	chr8_104423_C_to_A,chr8_104426_A_to_C	
...
```
Note: template_sequence and synthetic_errors were obtained from REDI sequencing of the barcode locus prior to the WGS. <br><br>
file 4. off_target.txt (list of off-target location predicted from edit distance <= 3):
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
#### 1b. Run pipeline
We use a HPC cluster for computing the task, for which the configuration of job submission and control are saved in profile/config.yaml. 
```
cd 1_outcome_precalling
snakemake --profile profile
```


### Pipeline 2: Target site sequence and genomic feature annotations
For a given list of PAM sites (and gRNA sequences), our pipeline 2 outputs a full list of annotations that are used for constructing the predictive model used in Pipeline 3. These annotations includes target sequence features, chromatin accessibility, histone modification, as well as the customized sequence repetitiveness metrics (and more details about feature list can be found in our publication). 
#### 1a. Prepare annotation files
Run the script to download pre-processed annotation files for the SCORE model (yeast chromatin, histone modification data). Please find the information about how the processed data were generated from 2_target_list_annotation/readmeplease.txt. 
```
cd 2_target_list_annotation/
sh download_annotation.sh
```

#### 1b. Modify the configuration files
```
ls inputs/
```
file 1. input.yaml (required by snakemake for locating file paths):
```
---
targetsheet: 'inputs/targets.txt'
genome: 'inputs/fasta/yeast_reference.fa'
tmp_dir: './tmpdir'
output_dir: './output'
atac_seq_bed: 'inputs/data/lin_all.merge.ins.bedgraph'
centromere_coords: 'inputs/data/yeast_CEN.coords'
transcription_forward_strand: 'inputs/data/WT.polyA.forward.depth'
transcription_reverse_strand: 'inputs/data/WT.polyA.reverse.depth'
hist_mod_list: 'inputs/data/histone_modification_binned_bed_files.txt'
```
file 2. targets.txt (list of samples with corresponding gRNA targets):
```
sample	chr	strand	coord	gRNA_sequence
MAGESTIC_REDI_1	chr8	-	104423	CTCCTCAGTTTGTCCACGGT
MAGESTIC_REDI_10	chr4	+	1503892	AGAGGAAGTCTCAACGGCAG
MAGESTIC_REDI_100	chr11	+	220649	ATACTGGCCACGTTTGACAA
MAGESTIC_REDI_1000	chr10	+	92631	TTGTGATTTTATTGATTCTG
MAGESTIC_REDI_1001	chr6	+	20768	GGTAACAAAGTCACGGCTCC
```

#### 1c. Run pipeline
We use a HPC cluster for computing the task, for which the configuration of job submission and control are saved in profile/config.yaml. 
```
cd 2_target_list_annotation/
snakemake --profile profile
```

### Pipeline 3: Machine learning model construction and tuning
We use the R caret package to perform model searching across a number of hyperparameters: 1) ML algorithms (gradient-boosted trees, random forest, L1- or L2-regularized regressions; 2) Hyperparamters within each ML algorithm (e.g. tree number ...); 3) Feature subsets (target sequence, chromatin context and genomic repetitiveness); 4) Data imbalance treatment (over- or down-sampling, with SMOTE or ADASYN). Repeated k-fold cross validation was used to assess the model performance. 

### Pipeline 4: Bulk whole-genome sequencing analysis for structral variant characterizations
#### 1a. Modify the configuration files
file 1. input.yaml (required by snakemake for locating file paths):
```
---
observed_outcome: 'inputs/edit_outcome.txt'
feature_table: 'inputs/combined_features.txt'
tmp_dir: './tmpdir'
output_dir: './output'
```
file 2. edit_outcome.txt (measured editing oucomes for samples encoded as 0/1 matrix, where 1 represents a matched unwanted outcome):
```
New_id	Editing_outcome_NE	Editing_outcome_SV	Editing_outcome_LD	Editing_outcome_NRecT
MAGESTIC_REDI_1	0	0	0	0
MAGESTIC_REDI_10	NA	NA	NA	NA
MAGESTIC_REDI_100	NA	NA	NA	NA
MAGESTIC_REDI_1000	1	0	0	0
MAGESTIC_REDI_1001	1	0	0	0
MAGESTIC_REDI_1002	0	0	0	0
MAGESTIC_REDI_1003	0	0	0	0
```
file 3. combined_features.txt (feature table obtained from Pipeline 2, with additional columns in the end to assign editing system and variant type as correction factors)
```
id	ATACseq_ins_freq	bp_to_telomere	bp_to_centromere	norm_distance_to_telomere	H2AK5ac	H2AS129ph	H3K4ac	H3K4me	H3K4me2	H3K4me3	H3K9ac	H3K14ac	H3K18ac	H3K23ac	H3K27ac	H3K36meH3K36me2	H3K36me3	H3K56ac	H3K79me	H3K79me3	H3S10ph	H4K5ac	H4K8ac	H4K12ac	H4K15acH4K20me	H4R3me	H4R3me2s	HTZ1	IGR100	IGR1000	IGR10000	IGR50	IGR500	IGR5000	ILR100	ILR1000	ILR10000	ILR50	ILR500	ILR5000	m0	m1	m2	m3	T_score	target_sequence_30mer	Transcription_Cas9_bind_template_strand	Transcription_Cas9_bind_non_template_strand	TRL2500x2	TRL250x2	TRL25x2	TRL5000x2	TRL500x2	TRL50x2	Editing_system_1	Editing_system_2	Editing_system_3	variant_type_SNV	variant_type_INDEL	variant_type_MNV
MAGESTIC_REDI_1	3.06	104423	1163	0.989	0.571911931	-1.413944272	1.82410521	-1.913294733	-1.06959095	2.134772687	1.647827909	1.381763726	2.527676605	1.081172126	0.557087315	-0.340625001	-0.460996809	-0.281504667	1.408581249	0.038505213	-0.193029944	0.632719882	0.000316058	-0.412611799	-0.306177367	-0.719997155	-1.018213765	-0.464537889	-1.073701456	-0.6606927	0.005	0.135	0.048	0	0.013	0.086	-0.00110.0166	0.0122	-0.0114	0.0314	0.0089	1	0	0	0	3.5	ATTTCCTCCTCAGTTTGTCCACGGTGGGTT	-5.644	9.1	17	10	0	21	17	6	1	0	0	0	0	1
```
