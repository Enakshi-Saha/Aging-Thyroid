# Sex agnostic aging genes
tb <- read.csv("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_aging.txt", sep="")

increasing_gene.GTEx = tb$gene_name[which(tb$P.Value<0.05 & tb$t>0)]
decreasing_gene.GTEx = tb$gene_name[which(tb$P.Value<0.05 & tb$t<0)]

GSE165724_tissue <- read.csv("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GSE165724_normalVSnormalAdjacent.txt",sep="")

increasing_gene = GSE165724_tissue$gene_name[which(GSE165724_tissue$P.Value<0.05 & GSE165724_tissue$t>0)]
decreasing_gene = GSE165724_tissue$gene_name[which(GSE165724_tissue$P.Value<0.05 & GSE165724_tissue$t<0)]

write.table(increasing_gene.GTEx, file="/home/esaha/Aging_thyroid/results/agingGenes_CLUEreg/GTEx_agingGenes_increasing.txt", row.names = F, col.names = F)
write.table(decreasing_gene.GTEx, file="/home/esaha/Aging_thyroid/results/agingGenes_CLUEreg/GTEx_agingGenes_decreasing.txt", row.names = F, col.names = F)
write.table(increasing_gene, file="/home/esaha/Aging_thyroid/results/agingGenes_CLUEreg/GSE165724_agingGenes_increasing.txt", row.names = F, col.names = F)
write.table(decreasing_gene, file="/home/esaha/Aging_thyroid/results/agingGenes_CLUEreg/GSE165724_agingGenes_decreasing.txt", row.names = F, col.names = F)

###################################################################################
# Sex-specific aging genes
# Load limma tables
tb_male <- read.csv("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_aging_male.txt", sep="")
tb_female <- read.csv("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_aging_female.txt", sep="")

male_increasing = tb_male$gene_name[which((tb_male$logFC >= 1) & (tb_male$P.Value < 0.05))]
male_decreasing = tb_male$gene_name[which((tb_male$logFC <= -1) & (tb_male$P.Value < 0.05))]

female_increasing = tb_female$gene_name[which((tb_female$logFC >= 1) & (tb_female$P.Value < 0.05))]
female_decreasing = tb_female$gene_name[which((tb_female$logFC <= -1) & (tb_female$P.Value < 0.05))]

write.table(male_increasing, file="/home/esaha/Aging_thyroid/results/agingGenes_CLUEreg/GTEx_agingGenes_male_increasing.txt", row.names = F, col.names = F)
write.table(male_decreasing, file="/home/esaha/Aging_thyroid/results/agingGenes_CLUEreg/GTEx_agingGenes_male_decreasing.txt", row.names = F, col.names = F)
write.table(female_increasing, file="/home/esaha/Aging_thyroid/results/agingGenes_CLUEreg/GTEx_agingGenes_female_increasing.txt", row.names = F, col.names = F)
write.table(female_decreasing, file="/home/esaha/Aging_thyroid/results/agingGenes_CLUEreg/GTEx_agingGenes_female_decreasing.txt", row.names = F, col.names = F)


##########################
# Find overlapping genes (same direction between sexes):
intersect(male_increasing, female_increasing)
intersect(male_decreasing, female_decreasing)

# Find conflicting genes (different direction between sexes):
intersect(male_increasing, female_decreasing)
intersect(male_decreasing, female_increasing)

###########################
# Check on which chromosomes the aging-associated genes are on
male_chr_increasing = tb_male$chr[which((tb_male$logFC >= 1) & (tb_male$P.Value < 0.05))]
male_chr_decreasing = tb_male$chr[which((tb_male$logFC <= -1) & (tb_male$P.Value < 0.05))]

female_chr_increasing = tb_female$chr[which((tb_female$logFC >= 1) & (tb_female$P.Value < 0.05))]
female_chr_decreasing = tb_female$chr[which((tb_female$logFC <= -1) & (tb_female$P.Value < 0.05))]

