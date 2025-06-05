library(ggplot2)
library(data.table)

male_xcell <- read.csv("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_thyroid_HashimotoExcluded_male_immuneInfiltration_aging.csv")
female_xcell <- read.csv("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_thyroid_HashimotoExcluded_female_immuneInfiltration_aging.csv")

# sig_cells = unique(c(male_xcell$X[which(male_xcell$p.value<0.05)],female_xcell$X[which(female_xcell$p.value<0.05)]))

age_coefficient = c(male_xcell$age.coefficient, female_xcell$age.coefficient)
sex = c(rep("male", nrow(male_xcell)), rep("female", nrow(female_xcell)))
sex = factor(sex, levels = c("female", "male"))
cell_type = c(male_xcell$X, female_xcell$X)
fulldata = data.frame(cell_type, sex, age_coefficient)
# Select cells that significantly change with age in either GTEx or TCGA
# fulldata = fulldata[which(fulldata$cell_type %in% sig_cells),]
head(fulldata)

ggplot(fulldata, aes(x=age_coefficient, y=cell_type, fill=sex)) + 
  geom_bar(position="dodge", stat="identity") + 
  labs(title="Change in Cell Composition with Age in GTEx",y="cell type", x = "t-statistic of age coefficient") +
  geom_vline(xintercept=qnorm(0.025), linetype="dashed",color = "red") +
  geom_vline(xintercept=qnorm(0.975), linetype="dashed",color = "red")

