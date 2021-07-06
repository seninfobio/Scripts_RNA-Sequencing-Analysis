# Scripts_RNA Sequencing Analysis

* The major task's are installation of tools in the respective environment operating system.
* The top recommendation is to install CONDA package to get easy accessiblity of the all the tools for RNA sequencing data.

#Installation


# download latest conda installer

$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

# run the installer
```
$ bash Miniconda3-latest-Linux-x86_64.sh

conda update --yes conda


./configure --prefix=/NABIC/HOME/senthil83/programs/mummer-4.0.0rc1
make
make install


conda install -c anaconda gcc_linux-64

conda install -c anaconda gxx_linux-64

```

# Step 1: Quality Control of RAW reads 

```

#--Make a script

vi run_fastqc.sh


#!/bin/bash

set -e

#--Quality control of rna seq data

fastqc \
--threads 16 \
-o /NABIC/HOME/senthil83/analysis/002_quality_control_rna_data \
/NABIC/HOME/senthil83/datafiles/001_sesamums_rna.seq_reads/*.gz 


#--Run

$bash run_fastqc.sh &> log &

exit to save the job in background

```


# Step 2: Trinity Installation

```
#-- Create a conda environment

conda create --name trinity_env

#--Activate the created environment

source activate trinity_env


#---Install trinity

conda install -c bioconda trinity

#--Check if trinity is working well

Trinity



#--exit from trinity environment

source deactivate trinity_env
```

# Step 3: Making script for running Trinity Program

#--Tutorial website  

https://biohpc.cornell.edu/lab/doc/trinity_workshop_part1.pdf

https://angus.readthedocs.io/en/2017/assembly-trinity.html

https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity


```
#--Run trinity for S. al datasets

cd /NABIC/HOME/senthil83/analysis/004_denovo_assembly_s_al_rna

#--Create a bash file 

vi run_trinity.sh

#!/bin/bash

set -e

source activate trinity_env

#--Create link file for input file

cd /NABIC/HOME/senthil83/analysis/004_denovo_assembly_s_al_rna

ln -fs /NABIC/HOME/senthil83/datafiles/001_sesamums_rna.seq_reads/S_ala_*.fastq.gz .

#--Run trinity 

Trinity \
--seqType fq \
--left S_ala_seed_1.fastq.gz,S_ala_stem_1.fastq.gz,S_ala_leaf_1.fastq.gz,S_ala_flower_1.fastq.gz,S_ala_capsule_1.fastq.gz \
--right S_ala_seed_2.fastq.gz,S_ala_stem_2.fastq.gz,S_ala_leaf_2.fastq.gz,S_ala_flower_2.fastq.gz,S_ala_capsule_2.fastq.gz \
--SS_lib_type RF \
--max_memory 500G \
--CPU 30 \
--output /NABIC/HOME/senthil83/analysis/004_denovo_assembly_s_alatum_rna/s_alatum_trinity_out

source deactivate trinity_env

#--END


$/usr/bin/time -o out.ram.time.txt -v bash run_trinity.sh &> log &


$jobs

#--Running should appear. If not, resolve the problem that occured. Please check the log file.

$exit

```

# Step 4: Installation of Mapping tools

```
#create conda environment

conda create --name bowtie_env

#--Activate the created environment

source activate bowtie_env


#---Install bowtie

conda install -c bioconda bowtie

#--Check if bowtie is working well

bo



#--exit from bowtie environment

source deactivate bowtie_env

```
# Assembler stat tool

```
cdconda install -c bioconda assembly-stats
conda install -c bioconda/label/cf201901 assembly-stats


#Run assembly stat

$assembly-stats file.fasta another_file.fastq

```
# Quality Assessment Tool for Genome Assemblies (BUSCO_secondary assembler)

```
# create specific environment for busco
conda create --name busco_env

conda install -c bioconda busco

conda install -c bioconda/label/broken busco

onda install -c bioconda/label/cf201901 busco

CondaSystemExit: Exiting.
```

# Quality Assessment Tool for Genome Assemblies

```
conda install -c bioconda quast
conda install -c bioconda/label/cf201901 quast

#-- Create a conda environment

conda create --name quast_env

#--Activate the created environment

source activate quast_env


#---Install quast

conda install -c bioconda quast

#--Check if quast is working well

quast



#--exit from quast environment

source deactivate quast_env

```
# CD_HIT (secondary assembler)_Installation

$conda create --name cd_hit_env

source activate cd_hit_env

$conda install -c bioconda cd-hit


https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide#CDHITEST

cd-hit-est -i est_human -o est_human95 -c 0.95 -n 10 -d 0 -M 16000 -T 8

https://www.researchgate.net/post/How_to_use_CD-hit_to_remove_the_duplicates

*You need to change the clustering threshold (-c) - at the moment you are clustering at 90% identity 
- if you want to remove all exact duplicates you need to use -c 1.0   

https://www.biostars.org/p/418872/

https://www.biostars.org/p/108519/


*Question: CD-HIT Clustering evaluation

https://www.biostars.org/p/425716/


cat out.clstr | awk ' /Cluster/ { no+=1;}; !/Cluster/ { id=substr($3, 2, length($3)-4); printf("%s\t%s\n", no, id) } '

justification of cd hit usage

https://www.biobam.com/cd-hit-clustering-omicsbox/?cn-reloaded=1

```
#--Running

cd /NABIC/HOME/senthil83/analysis/015_seed_transcriptomics/03.cd_hit_alatum_seed

source activate cd_hit_env

cd-hit-est -i Trinity.fasta -o est_alatum_seed -c 0.95 -n 10 -d 0 -M 16000 -T 8 &> log &

cd /NABIC/HOME/senthil83/analysis/015_seed_transcriptomics/06.cd_hit_angolense_seed

source activate cd_hit_env
cd-hit-est -i Trinity.fasta -o est_angolense_seed -c 0.95 -n 10 -d 0 -M 16000 -T 8 &> log &

```

# Step 5.1: Differential Gene Expression (DEG) analysis

* Install RSEM in Trinity_envs
Source activate trinity_env

* install
conda install -c biobuilds rsem


* vi run_rsem.sh

```

#!/bin/bash

set -e

#--Go to:QW the working directory


cd /NABIC/HOME/senthil83/analysis/005_RSEM_sesamum_alatum


#--Link the datasets
ln -fs /NABIC/HOME/senthil83/datafiles/001_sesamums_rna.seq_reads/S_ala_* .
ln -fs /NABIC/HOME/senthil83/analysis/004_denovo_assembly_s_alatum_rna/s_alatum_trinity_out/Trinity.fasta .

source activate trinity_env

#--Align and abundance estimation


for reads in *_1.fastq.gz

do

base=$(basename $reads _1.fastq.gz)

mkdir -p output_dir_${base}

align_and_estimate_abundance.pl --thread_count  96 \
--output_prefix ${base}_out_RSEM \
--output_dir output_dir_${base} \
--transcripts Trinity.fasta \
--seqType fq \
--left ${base}_1.fastq.gz \
--right ${base}_2.fastq.gz \
--est_method RSEM \
--aln_method bowtie \
--trinity_mode \
--prep_reference 

done

source deactivate trinity_env


#--Bye


$ /usr/bin/time -o out.time.ram.txt -v bash run_rsem.sh &> log.rsem &

```
# Step 5.2: Using RSEM.isoforms.results from different profiles 

# Data preparation

#---Create symbolic links for RSEM.isoforms.results from diverse organs

```

cd /NABIC/HOME/senthil83/analysis/006_edgeR_alatum

ln -fs /NABIC/HOME/senthil83/analysis/005_RSEM_sesamum_alatum/output_dir_S_ala_capsule/S_ala_capsule_out_RSEM.isoforms.results .
ln -fs /NABIC/HOME/senthil83/analysis/005_RSEM_sesamum_alatum/output_dir_S_ala_flower/S_ala_flower_out_RSEM.isoforms.results .
ln -fs /NABIC/HOME/senthil83/analysis/005_RSEM_sesamum_alatum/output_dir_S_ala_leaf/S_ala_leaf_out_RSEM.isoforms.results .
ln -fs /NABIC/HOME/senthil83/analysis/005_RSEM_sesamum_alatum/output_dir_S_ala_seed/S_ala_seed_out_RSEM.isoforms.results .
ln -fs /NABIC/HOME/senthil83/analysis/005_RSEM_sesamum_alatum/output_dir_S_ala_stem/S_ala_stem_out_RSEM.isoforms.results .

#---Merge them into one matrix (solution here : https://github.com/deweylab/RSEM#de)

source activate trinity_env


rsem-generate-data-matrix S_ala_capsule_out_RSEM.isoforms.results \
S_ala_flower_out_RSEM.isoforms.results \
S_ala_leaf_out_RSEM.isoforms.results \
S_ala_seed_out_RSEM.isoforms.results \
S_ala_stem_out_RSEM.isoforms.results > output_isoforms_alatum.matrix


#---Run DEG with edgeR | No replicate

run_DE_analysis.pl --matrix output_isoforms_alatum.matrix --method edgeR  --dispersion 0.2 --output /NABIC/HOME/senthil83/analysis/006_edgeR_alatum 

```

# Step 5.3: To visualize in Heat map using R Program

```

# Installation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")

# Libraries
library(edgeR)

#---General descriptive analysis

# Data importation
count_data = read.csv("iso_counts.csv", h = T, sep = ",", row.names = 1)
group = factor(c("cap","flo","lea","see","ste" ))

# Create the list file
dge <- DGEList(counts=count_data, group=group)

# Calculating library sizes from column totals.
dim(dge)

# Normalize by total count using exact test with a bcv = 0.2
bcv = 0.2
et <- exactTest(dge, dispersion=bcv^2)
et

# Multidimention scaling plot

plotMDS(dge)

# plor BCV
plotBCV(dge)


# MD
plotMD(dge)


# PLOT SMEAR

plotSmear(dge)


# heatmap
library(RColorBrewer)
logcpm <- cpm(dge, log=TRUE)
heatmap(logcpm)
heatmap(logcpm, col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))

```

