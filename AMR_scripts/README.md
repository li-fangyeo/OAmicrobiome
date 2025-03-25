# Antimicrobial Resistance Analysis

## Overview
This folder contains scripts involved in the antimicrobial resistance (AMR) analysis using metagenomics sequencing data. The workflow involves aligning sequencing reads to resistance gene databases, analyzing and visualizing patterns.

### 1. Data Processing  
- Input: `fastq.gz` files containing sequencing reads.  
- Reads were aligned to the **ResFinder** and **ResFinderFG** databases separately using **KMA**.  
- Alignment outputs were concatenated separately for **ResFinder** and **ResFinderFG**.  
- **ResFinder** results were used for the primary analysis, while **ResFinderFG** results were used in a complementary analysis.  

### 2. Data Analysis using R
- ResFinder result were used to generate:
  - **Distribution Box Plot**: Visualizing the distribution of antimicrobial resistance genes RPKM and Shannon index.
  - **Heatmap of Gene RPKM**: Displaying gene abundance levels related to specific drug classes.
  - **Venn Diagram**: Identifying shared and unshared resistance genes.
- ResFinderFG result were used to generate:
  - **Stacked Bar Plot**: Visualize the proportion of gene related to specific drug class in the total gene pool.

## Tools
- **KMA (Version 1.4.15)**: For sequence alignment to AMR databases.
- **R (Version 4.4.0) & RStudio**: For data visualization and statistical analysis.
- **ResFinder (Version 4.0) & ResFinderFG (Version 2.0)**: Databases used for resistance gene identification.

## Running the Pipeline  

### 0. Installation & Setup
KMA was installed with following command:
```bash
# Create and activate a conda environment
conda create --name kma
conda activate kma

# Clone and compile KMA
cd ~/project_directory
git clone https://bitbucket.org/genomicepidemiology/kma.git
cd kma && make
```

Resfinder database was downloaded and installed with following command:
```bash
# Setting up Resfinder script and database
cd ~/project_directory
git clone https://git@bitbucket.org/genomicepidemiology/resfinder.git
cd resfinder
git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git db_resfinder

# KMA indexing
cd ~/project_directory/resfinder/db_resfinder
mkdir kma_index
cd ~/project_directory/kma
./kma index -i ~/project_directory/resfinder/db_resfinder/all.fsa -o ~/project_directory/resfinder/db_resfinder/kma_index/resfinder

# Create an output directory for batch run
cd ~/project_directory/
mkdir kma_out
```

ResfinderFG database was downloaded and installed with following command:
```bash
# Download database
cd ~/project_directory
mkdir resfinderFG
cd resfinderFG
wget https://github.com/RemiGSC/ResFinder_FG_Construction/blob/main/output/RFG_db/ResFinder_FG.fasta

# Index ResFinderFG database
mkdir resfinderFG/kma_index
./kma index -i ~/project_directory/resfinderFG/ResFinder_FG.fasta -o ~/project_directory/resfinderFG/kma_index/resfinderFG
```

### 1. Sequence Alignment with KMA  
Running KMA:
```bash
./kma -ipe <input_read_1.fastq.gz> <input_read_2.fastq.gz> \
  -o <output_directory/sample_name> \
  -t_db <index_directory/index_prefix> \
  -1t1 -nc -na -nf -ef -ID 80 -ml 60
```
Example, running against Resfinder:
```bash
./kma -ipe sample1_r1.fastq.gz sample1_r2.fastq.gz \
  -o ~/project_directory/kma_out/sample1 \
  -t_db ~/project_directory/resfinder/db_resfinder/kma_index/resfinder \
  -1t1 -nc -na -nf -ef -ID 80 -ml 60
```

**Note:**  
- KMA is not added to the system path, so it should be run by specifying the path.  

**Batch Processing:** 
- A bash script (batch_kma.sh) is provided in the repository to automate batch runs.  

### 2. Data Analysis in R  
- Run the provided R script to generate plots.
- ResFinder result:
  - **00_concatenate.R**: Concatenate individual outputs
  - **01_diversity.R**: Generate RPKM and Shannon index distribution boxplots.
  - **02_heatmap.R**: Generate heatmap showing abundance level related to specific drug class.
  - **03_venn.R**: Generate Venn diagram to visualized the number of shared and unshared gene.
- ResFinderFG result:
  - **barPlot_resfinderFG.R**: Concatenate the individual outputs and a generate stacked bar plot to show proportion related to specific drug class.

## Citation
If you use this workflow, please cite the relevant tools and databases.
