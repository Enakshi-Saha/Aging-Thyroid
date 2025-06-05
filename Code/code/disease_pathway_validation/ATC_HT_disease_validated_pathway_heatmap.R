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

# Get pathways validated for at least 1 disease between HD and ATC
validated_pathways = unique(c(hashimoto_info$validated_pathways, ATC_info$validated_pathways))

################ Heatmap with pathways arranged in KEGG Hierarchical groups ######
# Load KEGG pathways
pathways <- gmtPathways("/home/ubuntu/c2.cp.kegg.v7.1.symbols.gmt")

# Load KEGG pathway hierarchy
pathway_hierarchy = get(load("/home/esaha/Aging_thyroid/data/KEGG_pathway_hierarchy_byFunction.RData"))

pathway_superclass = c("metabolism",
                       "genetic information processing",
                       "cell signaling, cell growth and death",
                       "cell communication",
                       "immune system and immune diseases",
                       "endocrine, circulatory, nervous system and related diseases",
                       "metabolic diseases",
                       "infectious diseases")

pathway_ordering = c(2,4,14,15,17,20,23,21, 24,25, 34, 26, 27, 29, 33, 35, 37, 38)
pathway_order = unlist(pathway_hierarchy[pathway_ordering])
pathway_order = unname(pathway_order[which(pathway_order %in% validated_pathways)])

# Table for NES
h = hashimoto_gsea$NES[match(pathway_order,hashimoto_gsea$pathway)]
atc = ATC_gsea$NES[match(pathway_order,ATC_gsea$pathway)]

tab = data.frame(cbind(h,atc))
colnames(tab) = c("HT", "ATC")
rownames(tab) = stringr::str_to_title(lapply(strsplit(pathway_order, split = "_"), function(x){paste(x[-1], collapse = " ")}))
head(tab)

# A matrix of characters (stars) to indicate which results are significant
lab = data.frame(matrix(NA, nrow(tab), ncol(tab)))
rownames(lab) = pathway_order
colnames(lab) = c("HT", "ATC")
lab$HT[which(rownames(lab) %in% hashimoto_info$validated_pathways)] = "*"
lab$ATC[which(rownames(lab) %in% ATC_info$validated_pathways)] = "*"

rownames(lab) = stringr::str_to_title(lapply(strsplit(rownames(lab), split = "_"), function(x){paste(x[-1], collapse = " ")}))
head(lab)

# Make heatmap comparing all groups

library(ggplot2)
library(gplots)
library(RColorBrewer)

mycol1 <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
pdf("/home/esaha/Aging_thyroid/results/plots/heatmap/ATC_HT_diseaseEffect_validatedPathways_heatmap.pdf",h=25,w=20)
heatmap.2(as.matrix(tab),density.info="none",trace="none",col=mycol1,symbreaks=T,symkey=T, Rowv = F,
          cexRow=2, cexCol=2, srtCol = 0, keysize=0.5, mar=c(5,90), key.title=NULL, key.xlab="NES", cellnote=lab, notecol="black", notecex=5)
dev.off()

###############################
# Order pathways by age direction
# Load GTEx GSEA age tables 
GTEx = get(load("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/GTEx_Thyroid_aging_GSEA.RData"))
paths = stringr::str_to_title(lapply(strsplit(GTEx$pathway, split = "_"), function(x){paste(x[-1], collapse = " ")}))
pos_paths = paths[which(GTEx$NES>0)]
neg_paths = paths[which(GTEx$NES<=0)]
relevant_pos_paths = intersect(pos_paths, rownames(tab))
relevant_neg_paths = intersect(neg_paths, rownames(tab))
tab0 = tab[match(c(relevant_pos_paths,relevant_neg_paths),rownames(tab)),]

# Different color pallette to distinguish disease effect from age effect
mycol2 <- colorRampPalette(rev(brewer.pal(11,"PiYG")))(50)
pdf("/home/esaha/Aging_thyroid/results/plots/heatmap/ATC_HT_diseaseEffect_validatedPathways_heatmap_differentColor_ageOrdered.pdf",h=25,w=20)
heatmap.2(as.matrix(tab0),density.info="none",trace="none",col=mycol2,symbreaks=T,symkey=T, Rowv = F, Colv=F,
          cexRow=2, cexCol=2, srtCol = 0, keysize=0.5, mar=c(5,90), key.title=NULL, key.xlab="NES", cellnote=lab, notecol="black", notecex=5)
dev.off()

# Heatmap with pathways significant for both age and disease
paths = stringr::str_to_title(lapply(strsplit(GTEx$pathway, split = "_"), function(x){paste(x[-1], collapse = " ")}))
sig_pos_paths = paths[which(GTEx$NES>0 & GTEx$padj < 0.05)]
sig_neg_paths = paths[which(GTEx$NES<=0 & GTEx$padj < 0.05)]
sig_relevant_pos_paths = intersect(sig_pos_paths, rownames(tab))
sig_relevant_neg_paths = intersect(sig_neg_paths, rownames(tab))
tab1 = tab[match(c(sig_relevant_pos_paths,sig_relevant_neg_paths),rownames(tab)),]
lab1 = lab[match(c(sig_relevant_pos_paths,sig_relevant_neg_paths),rownames(lab)),]

# Different color pallette to distinguish disease effect from age effect
mycol2 <- colorRampPalette(rev(brewer.pal(11,"PiYG")))(50)
pdf("/home/esaha/Aging_thyroid/results/plots/heatmap/ATC_HT_diseaseEffect_validatedPathways_heatmap_differentColor_ageOrdered_ageSignificant.pdf",h=25,w=20)
heatmap.2(as.matrix(tab1),density.info="none",trace="none",col=mycol2,symbreaks=T,symkey=T, Rowv = F, Colv=F,
          cexRow=3, cexCol=3, srtCol = 0, keysize=0.5, mar=c(5,90), key.title=NULL, key.xlab="NES", cellnote=lab1, notecol="black", notecex=5)
dev.off()

##########################################
# Order pathways by sex direction
# Load GTEx GSEA age tables

GTEx = get(load("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/GTEx_Thyroid_sexDifference_GSEA.RData"))
paths = stringr::str_to_title(lapply(strsplit(GTEx$pathway, split = "_"), function(x){paste(x[-1], collapse = " ")}))
pos_paths = paths[which(GTEx$NES>0)]
neg_paths = paths[which(GTEx$NES<=0)]
relevant_pos_paths = intersect(pos_paths, rownames(tab))
relevant_neg_paths = intersect(neg_paths, rownames(tab))
tab0 = tab[match(c(relevant_pos_paths,relevant_neg_paths),rownames(tab)),]

# Different color pallette to distinguish disease effect from age effect
mycol2 <- colorRampPalette(rev(brewer.pal(11,"PiYG")))(50)
pdf("/home/esaha/Aging_thyroid/results/plots/heatmap/ATC_HT_diseaseEffect_validatedPathways_heatmap_differentColor_sexOrdered.pdf",h=25,w=20)
heatmap.2(as.matrix(tab0),density.info="none",trace="none",col=mycol2,symbreaks=T,symkey=T, Rowv = F, Colv=F,
          cexRow=2, cexCol=2, srtCol = 0, keysize=0.5, mar=c(5,90), key.title=NULL, key.xlab="NES", cellnote=lab, notecol="black", notecex=5)
dev.off()

# Heatmap with pathways significant for both sex and disease
paths = stringr::str_to_title(lapply(strsplit(GTEx$pathway, split = "_"), function(x){paste(x[-1], collapse = " ")}))
sig_pos_paths = paths[which(GTEx$NES>0 & GTEx$padj < 0.05)]
sig_neg_paths = paths[which(GTEx$NES<=0 & GTEx$padj < 0.05)]
sig_relevant_pos_paths = intersect(sig_pos_paths, rownames(tab))
sig_relevant_neg_paths = intersect(sig_neg_paths, rownames(tab))
tab1 = tab[match(c(sig_relevant_pos_paths,sig_relevant_neg_paths),rownames(tab)),]
lab1 = lab[match(c(sig_relevant_pos_paths,sig_relevant_neg_paths),rownames(lab)),]

# Different color pallette to distinguish disease effect from age effect
mycol2 <- colorRampPalette(rev(brewer.pal(11,"PiYG")))(50)
pdf("/home/esaha/Aging_thyroid/results/plots/heatmap/ATC_HT_diseaseEffect_validatedPathways_heatmap_differentColor_sexOrdered_sexSignificant.pdf",h=25,w=20)
heatmap.2(as.matrix(tab1),density.info="none",trace="none",col=mycol2,symbreaks=T,symkey=T, Rowv = F, Colv=F,
          cexRow=3, cexCol=3, srtCol = 0, keysize=0.5, mar=c(5,90), key.title=NULL, key.xlab="NES", cellnote=lab1, notecol="black", notecex=5)
dev.off()





