# Load limma tables
tb0 <- read.csv("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_aging.txt", sep="")

# Get cosmic genes
all_cosmic = read.csv("/home/ubuntu/cosmic_genes.csv")
oncogene = all_cosmic$Gene.Symbol[which(all_cosmic$Role.in.Cancer == "oncogene")]
TSG = all_cosmic$Gene.Symbol[which(all_cosmic$Role.in.Cancer == "TSG")]

# Get escape genes
escape_genes <- read.delim("/home/ubuntu/escape_genes.txt", header=FALSE)
colnames(escape_genes) <- escape_genes[1,]
escape_genes <- escape_genes[-1,]
head(escape_genes)

escape_gene_names = escape_genes$HUGO_gene_id

# get cosmic and escape genes that significantly change with age
aging_gene = tb0$gene_name[which(tb0$P.Value<0.05)]
aging_oncogene = intersect(oncogene, aging_gene)
aging_tsg = intersect(TSG, aging_gene)
aging_escape = intersect(escape_gene_names, aging_gene)


