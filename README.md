# GEA_Community_Detection

## Summary

This repository performs gene enrichment analysis using either the KEGG, Reactome, or PID ontologies. The experiment is set up to contain both a control and experimental arm where the control arm is enrichment of a gene list of m pathways using only n% of the genes in each pathway with a% additional random genes from the ontology. This gene list is then subjected to enrichment analysis and the relevant enriched pathways are determined. The overall power (or more precisely, the true positives) and average false positives are returned. The experimental condition is just like the control except that community detection is performed before enrichment analysis. In particular, one can select fastgreedy, walktrap, infomap, or multilevel as the possible clustering method. Again, power and average false positive are returned. 

## Data Download

Data for this project can be downloaded using the script and URL text file located in the Data Download Script folder. Please be sure to have the URL text file and bash script both in the same location. Once the bash script is run, all data files will be downloaded into a folder called Data in the directory where you run the script.

## Contact

*About the code: Lia Harrington (lia.x.harrington.gr@dartmouth.edu)

*About the project or collaboration: Jennifer Doherty, PhD (jennifer.a.doherty@dartmouth.edu) or Casey Greene at csgreene@mail.med.upenn.edu.
