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

tb <- read.csv("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_aging.txt", sep="")

increasing_gene = rownames(tb)[which(tb$P.Value<0.05 & tb$t>0)]
decreasing_gene = rownames(tb)[which(tb$P.Value<0.05 & tb$t<0)]

# composite aging score with all aging genes
composite_score = prcomp(t(indegree[which(rownames(indegree) %in% c(increasing_gene,decreasing_gene)),]))$x[,1]

indegree_increasing = apply(indegree[which(rownames(indegree) %in% increasing_gene),], MARGIN=2, mean)
indegree_decreasing = apply(indegree[which(rownames(indegree) %in% decreasing_gene),], MARGIN=2, mean)

fulldata = data.frame(age)
fulldata$smoking_status = as.character(smoking_status)
fulldata$sex = gender
fulldata$race = race

fulldata$indegree_increasing = indegree_increasing
fulldata$indegree_decreasing = indegree_decreasing
fulldata$composite_score = composite_score

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

############################
# Age trajectory by sex

# Summarize data for each age
fulldata_aggregated = fulldata %>% group_by(age_group, sex)  %>%
  summarise(mean_pathScore = mean(indegree_increasing))
head(fulldata_aggregated)

fulldata_aggregated$sex = factor(fulldata_aggregated$sex, levels = c("FEMALE", "MALE"))

ggplot(fulldata_aggregated, aes(x = age_group , y = mean_pathScore, color = sex)) +
  geom_point(size = 0.5) + geom_smooth(method = "lm") + labs(title="GTEx: Genes increasingly targeted with age ",
                                                             y="gene score", x = "age")

# Summarize data for each age
fulldata_aggregated = fulldata %>% group_by(age_group, sex)  %>%
  summarise(mean_pathScore = mean(indegree_decreasing))
head(fulldata_aggregated)

fulldata_aggregated$sex = factor(fulldata_aggregated$sex, levels = c("FEMALE", "MALE"))

ggplot(fulldata_aggregated, aes(x = age_group , y = mean_pathScore, color = sex)) +
  geom_point(size = 0.5) + geom_smooth(method = "lm") + labs(title="GTEx: Genes decreasingly targeted with age ",
                                                             y="gene score", x = "age")


# Summarize data for each age
fulldata_aggregated = fulldata %>% group_by(age_group, sex)  %>%
  summarise(mean_pathScore = mean(composite_score))
head(fulldata_aggregated)

ggplot(fulldata_aggregated, aes(x = age_group , y = mean_pathScore, color = sex)) +
  geom_point(size = 0.5) + geom_smooth(method = "auto") + labs(title="GTEx: change in aging_score with age ",
                                                               y="gene score", x = "age")
