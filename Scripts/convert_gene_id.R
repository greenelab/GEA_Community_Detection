library(org.Hs.eg.db)

genes_df = read.csv('hgsc_cluster_genes.tsv', sep = '\t', row.names = 1)

all_genes = row.names(genes_df)

mapping = select(org.Hs.eg.db, all_genes, c("ENTREZID", "GENENAME"), "SYMBOL")

genesNA = mapping[is.na(mapping$ENTREZID),]$SYMBOL

mappingNA = select(org.Hs.eg.db, genesNA, c("ENTREZID", "GENENAME"), "ALIAS")

alias_lst = mappingNA$ENTREZID
names(alias_lst) = mappingNA$ALIAS

for (s in 1:length(mapping$SYMBOL)){
  symbol = mapping$SYMBOL[s]
  if (is.element(symbol, names(alias_lst))){
  mapping$ENTREZID[s] = alias_lst[symbol]
  }
}

mapping [is.na(mapping$ENTREZID),]$ENTREZID

row.names(genes_df) = mapping$ENTREZID

write.table(genes_df, 'entrezid_hgsc.txt', sep = '\t')

