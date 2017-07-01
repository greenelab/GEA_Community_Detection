#!/bin/bash

# Example of GEA community detection simulation and application to HGSC data
set -o errexit

### -------- Parameters
iterations = 100
num_paths_min = 2
num_paths_max = 8
percent_min = .3
percent_max = 1
addit_min = .1
addit_max = 1
file_name = 'all_iterations.csv'

### -------- Simulation pipeline

# Prepare ontology and variable files
python --Scripts/Ontology_Prep.py

# Run simulations
python --Scripts/Enrichment_Main.py data_generation $iterations $num_paths_min $num_paths_max $percent_min $percent_max $addit_min $addit_max $file_name

# Visualize results
Rscript Scripts/gea_paper_figures.R

### -------- HGSC pipeline

# Convert from Gene Symbols to Entrezid
Rscript --convert_gene_id.R

# Run HGSC analysis
python --hgsc_analysis.py
