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
GTEx_gsea = get(load("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_Thyroid_Aging_SexDiff_ageEffect_withoutSexInteraction_GSEA.RData"))
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

age_cuts = quantile(age, probs = seq(0, 1, 0.05))
age_cuts[1] = age_cuts[1]-1
age_cuts[length(age_cuts)] = age_cuts[length(age_cuts)] + 1
mean_interval = function(x){
  y = rep(0, (length(x)-1))
  for (i in 2:length(x)){
    y[i-1] = (x[i-1] + x[i])/2
  }
  return(y)
}
fulldata$age_group = cut(fulldata$age,
                         breaks=age_cuts,
                         labels=mean_interval(age_cuts))
fulldata$age_group = as.numeric(as.character(fulldata$age_group))
head(fulldata)

pdf("/home/esaha/Aging_thyroid/results/plots/ageTrajectory/KEGG_agingPathwaySuperclass_trajectory.pdf")
for (i in 1:length(pathway_hierarchy)){
  pathway_class = names(pathway_hierarchy)[i]
  
  pathList = pathway_hierarchy[[i]]
  
  # Keep only pathways that are associated to aging in GTEx
  pathList = intersect(pathList,sig_pathways)
  
  if (length(pathList) != 0){
    pathgenes = unlist(pathways[which(names(pathways) %in% pathList)])
    
    indegree_pathgenes = apply(indegree[which(geneNames %in% pathgenes),], MARGIN=2, mean)
    fulldata$indegree = indegree_pathgenes
    
    # Age trajectory by sex
    pathway_name_lowercase = tolower(pathway_class)
    # Summarize data for each age
    fulldata_aggregated = fulldata %>% group_by(age_group, sex)  %>%
      summarise(mean_pathScore = mean(indegree))
    head(fulldata_aggregated)
    
    fulldata_aggregated$sex = factor(fulldata_aggregated$sex, levels = c("FEMALE", "MALE"))
    
    p = ggplot(fulldata_aggregated, aes(x = age_group , y = mean_pathScore, color = sex)) +
      geom_point(size = 0.5) + geom_smooth(method = "lm") + labs(title=paste("GTEx: targeting of pathways associated to\n",pathway_name_lowercase, "with age"),
                                                                 y="gene score", x = "age")
    print(p)
  }
}
dev.off()


# Check which classes aging pathways belong to
lapply(pathway_hierarchy, function(x){intersect(sig_pathways,x)})

############################################
# Bigger superclass of pathway trajectory
# Pathway groups:
# 1. Metabolism (carb, energy, lipid, diabetes): hierarchy 2:13
# 2. Genetic & Environmental Information processing (Signal transduction, signaling molecules): hierarchy 14:20
# 3. Cell growth-death: hierarchy 22:23
# 4. Cell communication (transport, communication): hierarchy c(21, 24)
# 5. Immune system and immune diseases: hierarchy (25, 34)
# 6. Cancers including thyroid cancer: hierarchy 33
# 7. Nervous system and neurodegenerative diseases: hierarchy (29, 35)
# 8. Endocrine system and circulatory system: hierarchy 26:27
# 9. Metabolic and circulatory diseases: hierarchy c(36, 37)
# 10. Infectious diseases: hierarchy 38

pathway_superclass = c("metabolism",
                       "genetic & environmental Information Processing",
                       "cell growth and death",
                       "cell communication",
                       "immune system and immune diseases",
                       "cancers including thyroid cancer",
                       "nervous system and neurodegenerative diseases",
                       "endocrine system and cisculatory system",
                       "metabolic and circulatory diseases",
                       "infectious diseases")

pathway_superclass_index = c("2:13", "14:20", "22:23", "c(21,24)", "c(25,34)", "33", "c(29,35)", "26:27", "36:37", "38")
path_superclass = cbind(pathway_superclass, pathway_superclass_index)

pdf("/home/esaha/Aging_thyroid/results/plots/ageTrajectory/KEGG_agingPathway_10Superclass_trajectory.pdf")
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
    pathway_name_lowercase = tolower(pathway_class)
    # Summarize data for each age
    fulldata_aggregated = fulldata %>% group_by(age_group, sex)  %>%
      summarise(mean_pathScore = mean(indegree))
    head(fulldata_aggregated)
    
    fulldata_aggregated$sex = factor(fulldata_aggregated$sex, levels = c("FEMALE", "MALE"))
    
    p = ggplot(fulldata_aggregated, aes(x = age_group , y = mean_pathScore, color = sex)) +
      geom_point(size = 0.5) + geom_smooth(method = "auto") + labs(title=paste("GTEx: targeting of pathways with age: pathways associated to\n",pathway_name_lowercase, "with age"),
                                                                 y="gene score", x = "age")
    print(p)
  }
}
dev.off()



