library(ggplot2)
library(fgsea)

# Sex-specific Heatmaps
gender = "F"
# Load KEGG pathways
pathways <- gmtPathways("/home/ubuntu/c2.cp.kegg.v7.1.symbols.gmt")
path_ordering = names(pathways)

# Set datasetID to be one of these three: GSE29315, GSE33630, GSE27155
datasetID = "GSE29315"
directory=paste("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/", datasetID, "/", sep = "")
setwd(directory)

# Load GSEA tables 
files_ls = list.files(path=".", pattern=paste(datasetID,gender,sep="_"), all.files=TRUE,full.names=TRUE)
files_ls = files_ls[grep("RData",files_ls)]
tab = data.frame(matrix(0,nrow=length(path_ordering),ncol=length(files_ls)))
rownames(tab) = path_ordering
tab_pval = tab
sig_pathways = c()

for (i in 1:length(files_ls)){
  f = files_ls[i]
  disease = unlist(strsplit(f,split="_"))[4]
  fgsea = get(load(f))
  NES = fgsea$NES[match(path_ordering, fgsea$pathway)]
  p = fgsea$padj[match(path_ordering, fgsea$pathway)]
  tab[,i] = NES
  colnames(tab)[i] = disease
  tab_pval[,i] = p
  colnames(tab_pval)[i] = disease
  sig_pathways = c(sig_pathways, rownames(tab_pval)[which(p<0.05)])
}
sig_pathways = unique(sig_pathways)

# A matrix of characters (stars) to indicate which results are significant
lab = matrix(NA, nrow(tab_pval), ncol(tab_pval))
lab[which(tab_pval < 0.05)] = "*"
lab[which(tab_pval >= 0.05)] = ""
lab_subset = lab[match(sig_pathways, rownames(tab_pval)),]
tab_subset = tab[match(sig_pathways, rownames(tab_pval)),]

rownames(tab_subset) = stringr::str_to_title(lapply(strsplit(rownames(tab_subset), split = "_"), function(x){paste(x[-1], collapse = " ")}))
head(tab_subset)

# Make heatmap comparing all groups

library(ggplot2)
library(gplots)
library(RColorBrewer)


mycol2 <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
pdf(paste("/home/esaha/Aging_thyroid/results/plots/heatmap/", datasetID, "_", gender, "_diseaseEffect_heatmap.pdf", sep = ""),h=25,w=16)
heatmap.2(as.matrix(tab_subset),density.info="none",trace="none",col=mycol2,symbreaks=T,symkey=T, 
          cexRow=1, cexCol=1, srtCol = 45, keysize=0.5, mar=c(20,50), key.title=NULL, key.xlab="NES", cellnote=lab_subset, notecol="black", notecex=5)
dev.off()


