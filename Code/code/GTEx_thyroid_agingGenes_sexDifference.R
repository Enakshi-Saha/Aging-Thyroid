library(ggplot2)
library(fgsea)
library(gghalves)

tb <- read.csv("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_aging.txt", sep="")

increasing_gene.GTEx = tb$gene_name[which(tb$P.Value<0.05 & tb$t>0)]
decreasing_gene.GTEx = tb$gene_name[which(tb$P.Value<0.05 & tb$t<0)]

GTEx_sex <- read.csv("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_sexDiff_female.txt",sep="")

maintitle = "Aging gene targeting in GTEx: Female vs Male"
score_increasing = GTEx_sex$t[which(GTEx_sex$gene_name %in% increasing_gene.GTEx)]
score_decreasing = GTEx_sex$t[which(GTEx_sex$gene_name %in% decreasing_gene.GTEx)]
score = c(score_increasing, score_decreasing)
stype = c(rep("increasing", length(score_increasing)), rep("decreasing", length(score_decreasing)))
fulldata = data.frame(score)
fulldata$aging_trend = factor(stype, levels = c("increasing", "decreasing"))
ggplot(fulldata, aes(x = aging_trend, y = score, fill = aging_trend)) + geom_boxplot() + theme_bw() + ggtitle(maintitle) + xlab("") + ylab("t-statistic (Female - Male)") + geom_hline(yintercept=0, linetype="dashed",color = "red")  

ggplot(fulldata, aes(x = aging_trend, y = score, fill = aging_trend)) + geom_half_violin(side = "l") + geom_half_boxplot(side = "r")  + theme_bw() + ggtitle(maintitle) + xlab("") + ylab("t-statistic (Female - Male)") + geom_hline(yintercept=0, linetype="dashed",color = "red") +
  scale_fill_manual(values=c("red", "blue"))
