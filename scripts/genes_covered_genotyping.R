
gene_names <- lapply(genes$Gene.s., function(str) return(unique(strsplit(as.character(str), ',')[[1]])))
gene_names_unique <- unique(unlist(gene_names))


asd_genes_included <- top109_ASDgenes$GENE_NAME[top109_ASDgenes$GENE_NAME %in% gene_names_unique]

asd_genes_excluded <- top109_ASDgenes$GENE_NAME[!(top109_ASDgenes$GENE_NAME %in% gene_names_unique)]

#covers 103 of the 109 genes

markers_per_gene <- sapply(asd_genes_included, function(gene_name) return(sum(unlist(gene_names) == gene_name)))
                           
mean(markers_per_gene)