# Lia Harrington 2016 - Community Detection
# convert_gene_id.R
#
# Usage:
# Run in command line:
#
#       Rscript convert_gene_id.R
#
# Output:
# Converts Gene Symbols to Entrezid and saves as tsv file. 
#
# Notes: First all genes that have a one-one mapping are immediately mapped to their
# entrezid using the "SYMBOL" option in the select command. Of those genes that have a
# one-many mapping, the "ALIAS" option is used to coerce a one-one mapping. 

library(org.Hs.eg.db)

genes_df <- read.csv("./Data/hgsc_cluster_genes.tsv", sep = "\t", row.names = 1)

# list of all genes in the hgsc file
all_genes <- row.names(genes_df)

# performs mapping of all genes symbols to entrezid
mapping <- select(org.Hs.eg.db, all_genes, c("ENTREZID", "GENENAME"), "SYMBOL")

# any genes that do not initally map to an entrezid, thus the genes that have an 
# "NA" for entrezid s
genesNA <- mapping[is.na(mapping$ENTREZID),]$SYMBOL

# perform second stage mapping of the genes that did not initially map with "SYMBOL"
# option
mappingNA <- select(org.Hs.eg.db, genesNA, c("ENTREZID", "GENENAME"), "ALIAS")

# creates an association between the gene name and entrezid for those genes that did 
# not map using "SYMBOL" option
alias_lst <- mappingNA$ENTREZID
names(alias_lst) <- mappingNA$ALIAS

# iterates through the SYMBOL feature of the mapping object and for those genes that
# did not map through first pass, it replaces the "NA" in enterezid with the mapped
# alias enterezid 
for (s in 1:length(mapping$SYMBOL)){
  symbol = mapping$SYMBOL[s]
  if (is.element(symbol, names(alias_lst))){
  mapping$ENTREZID[s] <- alias_lst[symbol]
  }
}

# checks that mapping of all enterezid complete 
mapping [is.na(mapping$ENTREZID),]$ENTREZID

# sets the rownames of orignal genes_df to the new complete entrezid mapping 
row.names(genes_df) <- mapping$ENTREZID

# writes genes_df to tsv file
write.table(genes_df, "./Data/entrezid_hgsc.txt", sep = "\t")
