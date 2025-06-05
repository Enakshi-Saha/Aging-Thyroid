library(ggplot2)
library(gghalves)

tb_male <- read.csv("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_aging_male.txt", sep="")
tb_female <- read.csv("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_aging_female.txt", sep="")

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


maintitle = "GTEx: sex difference in the effect of aging on gene targeting"
score_tsg_male = tb_male$t[which(tb_male$gene_name %in% TSG)]
score_tsg_female = tb_female$t[which(tb_female$gene_name %in% TSG)]
score_onco_male = tb_male$t[which(tb_male$gene_name %in% oncogene)]
score_onco_female = tb_female$t[which(tb_female$gene_name %in% oncogene)]
score_escape_male = tb_male$t[which(tb_male$gene_name %in% escape_gene_names)]
score_escape_female = tb_female$t[which(tb_female$gene_name %in% escape_gene_names)]
score_others_male = tb_male$t[-which(tb_male$gene_name %in% c(all_cosmic$Gene.Symbol,escape_gene_names))]
score_others_female = tb_female$t[-which(tb_female$gene_name %in% c(all_cosmic$Gene.Symbol,escape_gene_names))]

# Plot with other genes (25363 genes)
score = c(score_tsg_male - score_tsg_female, score_onco_male -score_onco_female, score_escape_male-score_escape_female, score_others_male-score_others_female)
stype = c(rep("TSG", length(score_tsg_male)), rep("Oncogene", length(score_onco_male)), rep("Escape gene", length(score_escape_male)), rep("Other gene", length(score_others_male)))
stype = factor(stype, levels = c("Other gene", "Oncogene", "TSG", "Escape gene"))
fulldata = data.frame(score)
fulldata$gene_type = stype
ggplot(fulldata, aes(x = gene_type, y = score)) + geom_boxplot() + 
  theme_bw() + ggtitle(maintitle) + xlab("") + ylab("t-statistic (age coefficient)") + 
  geom_hline(yintercept=0, linetype="dashed", color = "red") 

# Plot without other genes
score = c(score_tsg_male, score_onco_male, score_escape_male, score_tsg_female, score_onco_female, score_escape_female)
stype = rep(c(rep("TSG", length(score_tsg_male)), rep("Oncogene", length(score_onco_male)), rep("Escape gene", length(score_escape_male))),2)
stype = factor(stype, levels = c("Oncogene", "TSG", "Escape gene"))
ngenes = length(score)/2
sex = c(rep("male", ngenes), rep("female", ngenes))
fulldata = data.frame(score)
fulldata$gene_type = stype
fulldata$sex = sex
ggplot(fulldata, aes(x = gene_type, y = score, fill = sex)) + geom_boxplot() + 
  theme_bw() + ggtitle(maintitle) + xlab("") + ylab("t-statistic (age coefficient)") + 
  geom_hline(yintercept=0, linetype="dashed", color = "red")

ggplot(fulldata, aes(x = gene_type, y = score, fill = sex)) + geom_half_violin(side = "l") + geom_half_boxplot(side = "r") +
  theme_bw() + ggtitle(maintitle) + xlab("") + ylab("t-statistic (age coefficient)") + 
  geom_hline(yintercept=0, linetype="dashed", color = "red")

# Plot with only oncogenes and TSGs
score = c(score_tsg_male, score_onco_male, score_tsg_female, score_onco_female)
stype = rep(c(rep("TSG", length(score_tsg_male)), rep("Oncogene", length(score_onco_male))),2)
stype = factor(stype, levels = c("Oncogene", "TSG"))
ngenes = length(score)/2
sex = c(rep("male", ngenes), rep("female", ngenes))
fulldata = data.frame(score)
fulldata$gene_type = stype
fulldata$sex = sex

ggplot(fulldata, aes(x = gene_type, y = score, fill = sex)) + geom_half_violin(side = "l") + geom_half_boxplot(side = "r") +
  theme_bw() + ggtitle(maintitle) + xlab("") + ylab("t-statistic (age coefficient)") + 
  geom_hline(yintercept=0, linetype="dashed", color = "red") + geom_jitter(size = 0.5)


#################################
wilcox.test(score_tsg_male)
wilcox.test(score_onco_male)
wilcox.test(score_escape_male)

wilcox.test(score_tsg_female, alternative = "less")
wilcox.test(score_onco_female, alternative = "less")
wilcox.test(score_escape_female, alternative = "less")

wilcox.test(score_tsg_male, score_others_male)
wilcox.test(score_onco_male, score_others_male)
wilcox.test(score_escape_male, score_others_male)

wilcox.test(score_tsg_female, score_others_female)
wilcox.test(score_onco_female, score_others_female)
wilcox.test(score_escape_female, score_others_female)

wilcox.test(score_tsg_male-score_tsg_female, score_others_male-score_others_female, alternative = "greater")
wilcox.test(score_onco_male-score_onco_female, score_others_male-score_others_female, alternative = "greater")
wilcox.test(score_escape_male-score_escape_female, score_others_male-score_others_female, alternative = "greater")
