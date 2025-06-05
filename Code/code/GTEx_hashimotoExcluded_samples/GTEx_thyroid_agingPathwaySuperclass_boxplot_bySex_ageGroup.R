# Sex effect on aging pathways
# to check if females "look older" after adjusting for age, smoking etc.
set.seed(1)
#if (!requireNamespace("BiocManager", quietly = TRUE))   
#  install.packages("BiocManager",repos = "http://cran.us.r-project.org")  
#BiocManager::install("fgsea")
#BiocManager::install("limma")
#BiocManager::install("Biobase")
#install.packages("ggplot2")
#install.packages("igraph")
library(data.table)
library(ggplot2)
library(dplyr)
library(fgsea)

# Load KEGG pathways
pathways <- gmtPathways("/home/ubuntu/c2.cp.kegg.v7.1.symbols.gmt")

# Load KEGG pathway hierarchy
pathway_hierarchy = get(load("/home/esaha/Aging_thyroid/data/KEGG_pathway_hierarchy_byFunction.RData"))

# Get GTEx aging pathways
# Load GSEA tables 
GTEx_gsea = get(load("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/GTEx_Thyroid_aging_GSEA.RData"))
sig_pathways = GTEx_gsea$pathway[which(GTEx_gsea$padj <0.05)]

# Get portal data and use only gtex portal samples
gtex.portal_thyroid <- read.csv("/home/ubuntu/gtex.portal_thyroid.txt", sep="")

##############################
bonobo = data.frame(fread("/home/esaha/BONOBO/network/BonoboPanda_GTEx_thyroid.txt"))
genes = bonobo$V1
bonobo = bonobo[,match(colnames(gtex.portal_thyroid), colnames(bonobo))]
bonobo_indices =  unlist(lapply(strsplit(colnames(bonobo), split = ".", fixed = T), function(x){paste(paste(x, collapse = "-"), 1, sep = ".")}))

# Get gene data
GTEx_thyroid = get(load("/home/esaha/BONOBO/data/GTEx_thyroid/GTEx_thyroid_RSE.RData"))
fdata.GTEx = as.data.frame(GTEx_thyroid@rowRanges)

# Discard Y genes
### Get GenecIDs for Y genes 
Y_genes.GTEx = fdata.GTEx$gene_id[which(fdata.GTEx$seqnames == "chrY")]
Y_genes.GTEx = gsub("\\..*", "", Y_genes.GTEx)

bonobo = bonobo[-which(genes %in% Y_genes.GTEx),]
genes = genes[-which(genes %in% Y_genes.GTEx)]

geneNames = fdata.GTEx$gene_name[match(genes,gsub("\\..*", "",fdata.GTEx$gene_id))]

# Get phenotypic data
phenotypes <- read.csv("~/BONOBO/data/GTEx_thyroid/thyroid_phenotypes.txt", sep="")
GTEx_phenotypes = as.data.frame(cbind(phenotypes$SEX, phenotypes$AGE, phenotypes$RACE, 
                                      phenotypes$BMI, phenotypes$TRISCH, 
                                      phenotypes$gtex.smrin, phenotypes$gtex.smnabtcht,
                                      phenotypes$MHSMKSTS))
rownames(GTEx_phenotypes) = rownames(phenotypes)
colnames(GTEx_phenotypes) = c("gender", "age", "race", "bmi", "ischemic_time", "rna_degrad", "batch", "smoking")
GTEx_phenotypes =  GTEx_phenotypes[match(bonobo_indices, rownames(GTEx_phenotypes)),]

# Define the covariates: gender, age, race etc
gender <- GTEx_phenotypes$gender
gender[which(gender == 1)] = "MALE"
gender[which(gender == 2)] = "FEMALE"
gender = factor(gender, levels = c("MALE", "FEMALE"))

age <- as.numeric(as.character(GTEx_phenotypes$age))
age[which(is.na(age))] <- mean(age,na.rm=TRUE)

# Categorize race 1, 4, 99 as single class "others", 2: "black", 3: "white"
race <- as.numeric(as.character(GTEx_phenotypes$race))
race[which(race != 2 & race != 3)] = "others"
race[which(race == 2)] = "black or african american"
race[which(race == 3)] = "white"
race <- as.factor(race)

bmi <- as.numeric(as.character(GTEx_phenotypes$bmi))
bmi[which(is.na(bmi))] <- mean(bmi,na.rm=TRUE)

ischemic_time = GTEx_phenotypes$ischemic_time
hour_minute = strsplit(ischemic_time, split=" ")
ischemic_timeH = sapply(hour_minute, function(x){as.numeric(as.character(x[1]))+as.numeric(as.character(x[1]))/60})

rna_degrad = as.numeric(as.character(GTEx_phenotypes$rna_degrad))
rna_degrad[which(is.na(rna_degrad))] <- mean(rna_degrad,na.rm=TRUE)

batch_effect = as.factor(GTEx_phenotypes$batch)

smoking_status = GTEx_phenotypes$smoking
smoking_status[which(is.na(smoking_status))] = "Unknown"
smoking_status = factor(smoking_status, levels = c("No", "Yes", "Unknown"))

###### Indegree matrix ######
indegree = bonobo
rownames(indegree) = genes


fulldata = data.frame(age)
fulldata$smoking_status = as.character(smoking_status)
fulldata$sex = factor(gender, levels = c("FEMALE", "MALE"))
fulldata$race = race

fulldata$age_group = fulldata$age
fulldata$age_group[which(fulldata$age <=52)] = "age<=52"
# fulldata$age_group[which(fulldata$age >40 & fulldata$age <=50)] = "age:40-50"
# fulldata$age_group[which(fulldata$age >50 & fulldata$age <=60)] = "age:50-60"
# fulldata$age_group[which(fulldata$age >60)] = "age>60"
fulldata$age_group[-which(fulldata$age <=52)] = "age>52"
fulldata$age_group = factor(fulldata$age_group, levels = c("age<=52", "age>52"))
head(fulldata)

############################################
# Bigger superclass of pathway trajectory
# Pathway groups:
# 1. Metabolism (carb, lipid): hierarchy c(2,4)
# 2. Genetic Information processing (transcription, translation, replication & repair): hierarchy c(14,15,17)
# 3. Cell signaling, cell growth and death (P53, Signaling molecules and interactions): hierarchy c(20,23)
# 4. Cell communication (transport, communication): hierarchy c(21, 24)
# 5. Immune system and immune diseases: hierarchy c(25, 34)
# 6. Endocrine system, circulatory system, Nervous system and related diseases (cancers: leukemia, neurodegenerative diseases): hierarchy c(26, 27, 29, 33, 35)
# 7. Metabolic diseases: hierarchy 37
# 8. Infectious diseases: hierarchy 38

pathway_superclass = c("metabolism",
                       "genetic information processing",
                       "cell signaling, cell growth and death",
                       "cell communication",
                       "immune system and immune diseases",
                       "endocrine, circulatory, nervous system",
                       "metabolic diseases",
                       "infectious diseases")

pathway_superclass_index = c("c(2,4)", "c(14,15,17)", "c(20,23)", "c(21, 24)", "c(25, 34)", "c(26, 27, 29, 33, 35)", "37", "38")
path_superclass = cbind(pathway_superclass, pathway_superclass_index)

pdf("/home/esaha/Aging_thyroid/results/plots/boxplot/GTEx_thyroid_hashimotoExcluded/KEGG_agingPathway_8Superclass_boxplot_bySex_ageGroup.pdf", 
    width=10,height=5)
for (i in 1:length(pathway_superclass)){
  pathway_class = pathway_superclass[i]
  
  j = eval(parse(text = path_superclass[i,2]))
  pathList = unlist(pathway_hierarchy[j])
  
  # Keep only pathways that are associated to aging in GTEx
  pathList = intersect(pathList,sig_pathways)
  
  if (length(pathList) != 0){
    pathgenes = unlist(pathways[which(names(pathways) %in% pathList)])
    
    indegree_pathgenes = apply(indegree[which(geneNames %in% pathgenes),], MARGIN=2, mean)
    fulldata$indegree = indegree_pathgenes
    
    # Age trajectory by sex
    pathway_name_uppercase = toupper(pathway_class)
    # boxplot without top & bottom 5% outliers
    p = ggplot(fulldata, aes(x = sex , y = indegree, fill = age_group)) +
      scale_fill_manual(values=c("darkgreen", "lightgreen")) +
      geom_boxplot() + labs(title=paste(pathway_name_uppercase),
                                                                   y="pathway indegree") +
      scale_y_continuous(limits = quantile(fulldata$indegree, c(0.05, 0.95))) +
      theme(text=element_text(size=22))
    print(p)
  }
}
dev.off()



