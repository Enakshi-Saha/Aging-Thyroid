######### Heatmap for disease effect #########
disease = "hashimoto"

library(ggplot2)
# Load GSEA tables
GTEx_age = get(load("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_Thyroid_Aging_SexDiff_ageEffect_withoutSexInteraction_GSEA.RData"))
GTEx_disease_male = get(load(paste("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_Thyroid_diseaseEffect_", disease, "_male_GSEA.RData", sep = "")))
GTEx_disease_female = get(load(paste("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_Thyroid_diseaseEffect_", disease, "_female_GSEA.RData", sep = "")))

tab = data.frame(GTEx_age$NES)
rownames(tab) = GTEx_age$pathway
tab$male = GTEx_disease_male$NES[match(rownames(tab), GTEx_disease_male$pathway)]
tab$female = GTEx_disease_female$NES[match(rownames(tab), GTEx_disease_female$pathway)]
colnames(tab) = c("age effect", "disease effect (male)", "disease effect (female)")
head(tab)

tab_pval = data.frame(GTEx_age$padj)
rownames(tab_pval) = GTEx_age$pathway
tab_pval$male = GTEx_disease_male$padj[match(rownames(tab_pval), GTEx_disease_male$pathway)]
tab_pval$female = GTEx_disease_female$padj[match(rownames(tab_pval), GTEx_disease_female$pathway)]
colnames(tab_pval) = c("age effect", "disease effect (male)", "disease effect (female)")
head(tab_pval)

# Heatmap with significant pathways
sig_pathways = GTEx_age$pathway[which(GTEx_age$padj <0.05)]

tab_subset = tab[match(sig_pathways, rownames(tab)),]
rownames(tab_subset) = stringr::str_to_title(lapply(strsplit(rownames(tab_subset), split = "_"), function(x){paste(x[-1], collapse = " ")}))
head(tab_subset)

# A matrix of characters (stars) to indicate which results are significant
lab = matrix(NA, nrow(tab_pval), ncol(tab_pval))
lab[which(tab_pval < 0.05)] = "*"
lab[which(tab_pval >= 0.05)] = ""
lab_subset = lab[match(sig_pathways, rownames(tab_pval)),]

# Make heatmap comparing all groups

library(ggplot2)
library(gplots)
library(RColorBrewer)


mycol2 <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
pdf(paste("/home/esaha/Aging_thyroid/results/plots/heatmap/GTEx_thyroid_diseaseEffect_", disease, "_bySex_heatmap.pdf", sep = ""),h=25,w=16)
heatmap.2(as.matrix(tab_subset),density.info="none",trace="none",col=mycol2,symbreaks=T,symkey=T, 
          cexRow=1, cexCol=1, srtCol = 45, keysize=0.5, mar=c(20,50), key.title=NULL, key.xlab="NES", cellnote=lab_subset, notecol="black", notecex=5)
dev.off()
