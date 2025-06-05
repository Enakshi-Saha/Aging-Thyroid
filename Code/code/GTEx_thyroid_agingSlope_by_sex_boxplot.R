library(ggplot2)
library(gghalves)
library(fgsea)

tb <- read.csv("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_aging.txt", sep="")

tb_male <- read.csv("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_aging_male.txt", sep="")
tb_female <- read.csv("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_aging_female.txt", sep="")

increasing_gene = rownames(tb)[which(tb$P.Value<0.05 & tb$t>0)]
decreasing_gene = rownames(tb)[which(tb$P.Value<0.05 & tb$t<0)]

maintitle = "Gene score in GTEx (genes with increasing targeting with age)"
score_female = tb_female$t[which(rownames(tb_female) %in% increasing_gene)]
score_male = tb_male$t[which(rownames(tb_male) %in% increasing_gene)]
score = c(score_female, score_male)
stype = c(rep("Female", length(score_female)), rep("Male", length(score_male)))
fulldata = data.frame(score)
fulldata$sex = factor(stype, levels = c("Female", "Male"))
ggplot(fulldata, aes(x = sex, y = score, fill = sex)) + geom_boxplot() + theme_bw() + ggtitle(maintitle) + xlab("") + ylab("t-statistic (age coefficient)") + geom_hline(yintercept=0, linetype="dashed",color = "red")  
ggplot(fulldata, aes(x = sex, y = score, fill = sex)) + geom_half_violin(side = "l") + geom_half_boxplot(side = "r") + theme_bw() + ggtitle(maintitle) + xlab("") + ylab("t-statistic (age coefficient)") + geom_hline(yintercept=0, linetype="dashed",color = "red")  

wilcox.test(score_female, score_male, alternative = "less")  
ks.test(score_female, score_male, alternative = "greater")  

maintitle = "Gene score in GTEx (genes with decreasing targeting with age)"
score_female = tb_female$t[which(rownames(tb_female) %in% decreasing_gene)]
score_male = tb_male$t[which(rownames(tb_male) %in% decreasing_gene)]
score = c(score_female, score_male)
stype = c(rep("Female", length(score_female)), rep("Male", length(score_male)))
fulldata = data.frame(score)
fulldata$sex = factor(stype, levels = c("Female", "Male"))
ggplot(fulldata, aes(x = sex, y = score, fill = sex)) + geom_boxplot() + theme_bw() + ggtitle(maintitle) + xlab("") + ylab("t-statistic (age coefficient)") + geom_hline(yintercept=0, linetype="dashed",color = "red")  
ggplot(fulldata, aes(x = sex, y = score, fill = sex)) + geom_half_violin(side = "l") + geom_half_boxplot(side = "r") + theme_bw() + ggtitle(maintitle) + xlab("") + ylab("t-statistic (age coefficient)") + geom_hline(yintercept=0, linetype="dashed",color = "red")  

wilcox.test(score_female, score_male, alternative = "greater")  
ks.test(score_female, score_male, alternative = "less")  

