library(ggplot2)
library(fgsea)

tb <- read.csv("/home/esaha/Aging_thyroid/results/limma_GSEA_tables/limma_GTEx_aging.txt", sep="")

tb_smoker <- read.csv("/home/esaha/Aging_thyroid/results/limma_GSEA_tables/limma_GTEx_aging_smoker.txt", sep="")
tb_nonsmoker <- read.csv("/home/esaha/Aging_thyroid/results/limma_GSEA_tables/limma_GTEx_aging_nonsmoker.txt", sep="")

increasing_gene = rownames(tb)[which(tb$P.Value<0.05 & tb$t>0)]
decreasing_gene = rownames(tb)[which(tb$P.Value<0.05 & tb$t<0)]

maintitle = "Gene score in GTEx (genes with increasing targeting with age)"
score_smoker = tb_smoker$t[which(rownames(tb_smoker) %in% increasing_gene)]
score_nonsmoker = tb_nonsmoker$t[which(rownames(tb_nonsmoker) %in% increasing_gene)]
score = c(score_smoker, score_nonsmoker)
stype = c(rep("smoker", length(score_smoker)), rep("nonsmoker", length(score_nonsmoker)))
fulldata = data.frame(score)
fulldata$sex = factor(stype, levels = c("smoker", "nonsmoker"))
ggplot(fulldata, aes(x = sex, y = score, fill = sex)) + geom_boxplot() + theme_bw() + ggtitle(maintitle) + xlab("") + ylab("t-statistic (age coefficient)") + geom_hline(yintercept=0, linetype="dashed",color = "red")  

wilcox.test(score_smoker, score_nonsmoker, alternative = "less")  

maintitle = "Gene score in GTEx (genes with decreasing targeting with age)"
score_smoker = tb_smoker$t[which(rownames(tb_smoker) %in% decreasing_gene)]
score_nonsmoker = tb_nonsmoker$t[which(rownames(tb_nonsmoker) %in% decreasing_gene)]
score = c(score_smoker, score_nonsmoker)
stype = c(rep("smoker", length(score_smoker)), rep("nonsmoker", length(score_nonsmoker)))
fulldata = data.frame(score)
fulldata$sex = factor(stype, levels = c("smoker", "nonsmoker"))
ggplot(fulldata, aes(x = sex, y = score, fill = sex)) + geom_boxplot() + theme_bw() + ggtitle(maintitle) + xlab("") + ylab("t-statistic (age coefficient)") + geom_hline(yintercept=0, linetype="dashed",color = "red")  

wilcox.test(score_smoker, score_nonsmoker, alternative = "greater")  

