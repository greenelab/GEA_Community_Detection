# GEA_Community_Detection

## Summary

This repository performs gene enrichment analysis using either the KEGG,
Reactome, or PID ontologies. The experiment is set up to contain both a control
and experimental arm where the control arm is enrichment of a gene list of m
pathways using only n% of the genes in each pathway with a% additional random
genes from the ontology. This gene list is then subjected to enrichment analysis
and the relevant enriched pathways are determined. The overall power (or more
precisely, the true positives) and average false positives are returned. The
experimental condition is just like the control except that community detection
is performed before enrichment analysis. In particular, one can select
fastgreedy, walktrap, infomap, or multilevel as the possible clustering method.
Again, power and average false positive are returned. 

## Data Download

Data for this project can be downloaded using the script and URL text file
located in the Data Download Script folder. This is required before running the
pipeline.

```bash
bash Data/data_files.sh
```

To reproduce all analyses including simulations and HGSC applications:

```bash
bash Scripts/gea_pipeline.sh
```

## Contact

*About the code: [Lia Harrington](lia.x.harrington.gr@dartmouth.edu)

*About the project or collaboration:
[Jennifer Doherty](jennifer.a.doherty@dartmouth.edu) or
[Casey Greene](csgreene@mail.med.upenn.edu).

