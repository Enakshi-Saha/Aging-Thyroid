

# Import validated pathway names for each disease
hashimoto_info = get(load("/home/esaha/Aging_thyroid/results/disease_specific_validated_pathways_list/hashimoto_validated_pathways.RData."))
PTC_info = get(load("/home/esaha/Aging_thyroid/results/disease_specific_validated_pathways_list/PTC_pathways.RData."))
ATC_info = get(load("/home/esaha/Aging_thyroid/results/disease_specific_validated_pathways_list/ATC_pathways.RData."))
FTC_info = get(load("/home/esaha/Aging_thyroid/results/disease_specific_validated_pathways_list/FTC_pathways.RData."))
follicular_adenoma_info = get(load("/home/esaha/Aging_thyroid/results/disease_specific_validated_pathways_list/follicular adenoma_validated_pathways.RData."))

# Import relevant datasets for NES values
hashimoto_gsea = get(load("~/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_Thyroid_diseaseEffect_hashimoto_GSEA.RData"))
PTC_gsea = get(load("~/Aging_thyroid/results/limma_GSEA_xcell_tables/GSE27155/sexAgnosticGSEA/GSE27155_diseaseEffect_Papillary Thyroid Carcinoma_GSEA.RData"))
ATC_gsea = get(load("~/Aging_thyroid/results/limma_GSEA_xcell_tables/GSE33630/sexAgnosticGSEA/GSE33630_diseaseEffect_anaplastic thyroid carcinoma (ATC)_GSEA.RData"))
FTC_gsea = get(load("~/Aging_thyroid/results/limma_GSEA_xcell_tables/GSE27155/sexAgnosticGSEA/GSE27155_diseaseEffect_Follicular Thyroid Carcinoma_GSEA.RData"))
follicular_adenoma_gsea = get(load("~/Aging_thyroid/results/limma_GSEA_xcell_tables/GSE29315/sexAgnosticGSEA/GSE29315_diseaseEffect_follicular adenoma_GSEA.RData"))

# Get pathways validated for at least 1 disease
validated_pathways = unique(c(hashimoto_info$validated_pathways, PTC_info$validated_pathways, 
ATC_info$validated_pathways, FTC_info$validated_pathways, follicular_adenoma_info$validated_pathways))

# Table for NES
h = hashimoto_gsea$NES[match(validated_pathways,hashimoto_gsea$pathway)]
ptc = PTC_gsea$NES[match(validated_pathways,PTC_gsea$pathway)]
atc = ATC_gsea$NES[match(validated_pathways,ATC_gsea$pathway)]
ftc = FTC_gsea$NES[match(validated_pathways,FTC_gsea$pathway)]
fa = follicular_adenoma_gsea$NES[match(validated_pathways,follicular_adenoma_gsea$pathway)]

tab = data.frame(cbind(h,ptc,atc,ftc,fa))
colnames(tab) = c("HT", "PTC", "ATC", "FTC", "FA")
rownames(tab) = stringr::str_to_title(lapply(strsplit(validated_pathways, split = "_"), function(x){paste(x[-1], collapse = " ")}))
head(tab)

# A matrix of characters (stars) to indicate which results are significant
lab = data.frame(matrix(NA, nrow(tab), ncol(tab)))
rownames(lab) = validated_pathways
colnames(lab) = c("HT", "PTC", "ATC", "FTC", "FA")
lab$HT[which(rownames(lab) %in% hashimoto_info$validated_pathways)] = "*"
lab$PTC[which(rownames(lab) %in% PTC_info$validated_pathways)] = "*"
lab$ATC[which(rownames(lab) %in% ATC_info$validated_pathways)] = "*"
lab$FTC[which(rownames(lab) %in% FTC_info$validated_pathways)] = "*"
lab$FA[which(rownames(lab) %in% follicular_adenoma_info$validated_pathways)] = "*"

rownames(lab) = stringr::str_to_title(lapply(strsplit(rownames(lab), split = "_"), function(x){paste(x[-1], collapse = " ")}))
head(lab)

# Make heatmap comparing all groups

library(ggplot2)
library(gplots)
library(RColorBrewer)

mycol2 <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
pdf("/home/esaha/Aging_thyroid/results/plots/heatmap/diseaseEffect_validatedPathways_heatmap.pdf",h=25,w=20)
heatmap.2(as.matrix(tab),density.info="none",trace="none",col=mycol2,symbreaks=T,symkey=T, 
          cexRow=1.4, cexCol=2, srtCol = 0, keysize=0.5, mar=c(20,30), key.title=NULL, key.xlab="NES", cellnote=lab, notecol="black", notecex=5)
dev.off()


