library(ggplot2)
library(fgsea)

#tb <- read.csv("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_aging.txt", sep="")

#increasing_gene.GTEx = tb$gene_name[which(tb$P.Value<0.05 & tb$t>0)]
#decreasing_gene.GTEx = tb$gene_name[which(tb$P.Value<0.05 & tb$t<0)]

GSE165724_tissue <- read.csv("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GSE165724_normalVSnormalAdjacent.txt",sep="")
# GSE165724_tissue <- read.csv("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GSE165724_normalVSnormalAdjacent_maleOnlyModel.txt",sep="")

GSE165724_age <- read.csv("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GSE165724_normalThyroid_aging.txt",sep="")
# GSE165724_age <- read.csv("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GSE165724_normalThyroid_aging_maleOnlyModel.txt",sep="")

increasing_gene = GSE165724_age$gene_name[which(GSE165724_age$P.Value<0.05 & GSE165724_age$t>0)]
decreasing_gene = GSE165724_age$gene_name[which(GSE165724_age$P.Value<0.05 & GSE165724_age$t<0)]

maintitle = "Aging gene targeting in GSE165724: normal vs adjacent normal"
score_increasing = GSE165724_tissue$t[which(rownames(GSE165724_tissue) %in% increasing_gene)]
score_decreasing = GSE165724_tissue$t[which(rownames(GSE165724_tissue) %in% decreasing_gene)]
score = c(score_increasing, score_decreasing)
stype = c(rep("increasing", length(score_increasing)), rep("decreasing", length(score_decreasing)))
fulldata = data.frame(score)
fulldata$aging_trend = factor(stype, levels = c("increasing", "decreasing"))
ggplot(fulldata, aes(x = aging_trend, y = score, fill = aging_trend)) + geom_boxplot() + theme_bw() + ggtitle(maintitle) + xlab("") + ylab("t-statistic (tumor adjacent normal - normal)") + geom_hline(yintercept=0, linetype="dashed",color = "red")  

ggplot(fulldata, aes(x = aging_trend, y = score, fill = aging_trend)) + geom_half_violin(side = "l") + geom_half_boxplot(side = "r")  + theme_bw() + ggtitle(maintitle) + xlab("") + ylab("t-statistic (tumor adjacent normal - normal)") + geom_hline(yintercept=0, linetype="dashed",color = "red") +
  scale_fill_manual(values=c("red", "blue"))
