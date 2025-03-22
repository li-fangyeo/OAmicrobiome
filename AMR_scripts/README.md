# Antimicrobial Resistance Analysis

## Overview
This folder consists of script involved in the antimicrobial resistance (AMR) analysis using metagenomics sequencing data. The workflow involves aligning sequencing reads to resistance gene databases, analyzing, and visualizing patterns.

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

### 1. Sequence Alignment with KMA  
**Batch Processing:** A bash script is provided in the repository to automate batch runs.  

**Note:**  
- The command for indexing will be updated soon.  
- KMA is not yet added to the system path, so it should be run by specifying the path.  

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
