#!/bin/bash

# Example of GEA community detection simulation and application to HGSC data
set -o errexit

### -------- Simulation pipeline

# Prepare KEGG ontology
python Scripts/ontology_prep.py c2.cp.kegg.v5.1.entrez.gmt.txt KEGG

# Run simulations
python Scripts/enrichment_main.py 100 2 8 .3 1 .1 1 all_iterations_data.csv 

# Visualize results
Rscript Scripts/gea_paper_figures.R

### -------- HGSC pipeline

# Prepare PID ontology 
python Scripts/ontology_prep.py PID.Entrez.DB.txt PID

# Convert from Gene Symbols to Entrezid
Rscript Scripts/convert_gene_id.R

# Run HGSC analysis
python Scripts/hgcs_analysis.py
