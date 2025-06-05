set.seed(1)
#if (!requireNamespace("BiocManager", quietly = TRUE))   
#  install.packages("BiocManager",repos = "http://cran.us.r-project.org")  
#BiocManager::install("fgsea")
#BiocManager::install("limma")
#BiocManager::install("Biobase")
#install.packages("ggplot2")
#install.packages("igraph")
library(dplyr)
library(limma)
library(Biobase)
library(ggplot2)
library(igraph)
library(data.table)
library(fgsea)
library(readr)
library(biomaRt)
library(edgeR)

# Get disease sample IDs
GTEx_thyroid_diseaseList <- readxl::read_excel("/home/ubuntu/GTEx_thyroid_diseaseList.xlsx")
hashimoto_samples = GTEx_thyroid_diseaseList$`Subject ID`[grep("hashimoto", GTEx_thyroid_diseaseList$`Pathology Categories`)]
goiter_samples = GTEx_thyroid_diseaseList$`Subject ID`[grep("goiter", GTEx_thyroid_diseaseList$`Pathology Categories`)]
adenoma_samples = GTEx_thyroid_diseaseList$`Subject ID`[grep("adenoma", GTEx_thyroid_diseaseList$`Pathology Categories`)]

# Get portal data and use only gtex portal samples
gtex.portal_thyroid <- read.csv("/home/ubuntu/gtex.portal_thyroid.txt", sep="")

##############################
bonobo = data.frame(fread("/home/esaha/Aging_thyroid/data/network/BonoboPanda_GTEx_thyroid_significantGene_subnetwork.txt"))
genes = bonobo$gene
TFs = bonobo$TF
rownames(bonobo) = paste(TFs, genes, sep = "_")
bonobo = bonobo[,match(colnames(gtex.portal_thyroid), colnames(bonobo))]
bonobo_indices =  unlist(lapply(strsplit(colnames(bonobo), split = ".", fixed = T), function(x){paste(paste(x, collapse = "-"), 1, sep = ".")}))

# Get gene data
GTEx_thyroid = get(load("/home/esaha/BONOBO/data/GTEx_thyroid/GTEx_thyroid_RSE.RData"))
fdata.GTEx = as.data.frame(GTEx_thyroid@rowRanges)

# Extract ESR and AR subnetworks
sub_net = bonobo[which(TFs %in% c("ESR1","ESR2","AR")),]
TF_name = TFs[which(TFs %in% c("ESR1","ESR2","AR"))]

# Get phenotypic data
phenotypes <- read.csv("~/BONOBO/data/GTEx_thyroid/thyroid_phenotypes.txt", sep="")
GTEx_phenotypes = as.data.frame(cbind(phenotypes$SUBJID,phenotypes$SEX, phenotypes$AGE, phenotypes$RACE, 
                                      phenotypes$BMI, phenotypes$TRISCH, 
                                      phenotypes$gtex.smrin, phenotypes$gtex.smnabtcht,
                                      phenotypes$MHSMKSTS))
rownames(GTEx_phenotypes) = rownames(phenotypes)
colnames(GTEx_phenotypes) = c("subject_ID","gender", "age", "race", "bmi", "ischemic_time", "rna_degrad", "batch", "smoking")
GTEx_phenotypes =  GTEx_phenotypes[match(bonobo_indices, rownames(GTEx_phenotypes)),]

# Define disease status: normal, hashimoto, goiter
disease = rep("normal", nrow(GTEx_phenotypes))
disease[which(GTEx_phenotypes$subject_ID %in% hashimoto_samples)] = "hashimoto"
disease[which(GTEx_phenotypes$subject_ID %in% goiter_samples)] = "goiter"
disease[which(GTEx_phenotypes$subject_ID %in% adenoma_samples)] = "adenoma"

disease = factor(disease, levels = c("normal", "hashimoto", "goiter", "adenoma"))

# Define the covariates: gender, age, race etc
gender <- GTEx_phenotypes$gender
gender[which(gender == 1)] = "MALE"
gender[which(gender == 2)] = "FEMALE"
gender = factor(gender, levels = c("FEMALE", "MALE"))

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

###### Fit limma model ######
#### Design Matrix
design = model.matrix(~ age + gender + race + bmi + ischemic_timeH + rna_degrad + batch_effect +
                        smoking_status + disease)

fit <- lmFit(bonobo, design)
fit <- eBayes(fit)

# Save table for disease effect
tb = topTable(fit,coef="diseasehashimoto",number=Inf)
tb$gene_id = lapply(strsplit(rownames(tb), split = "_"), function(x){x[2]})
tb$TF = lapply(strsplit(rownames(tb), split = "_"), function(x){x[1]})

gene_id = gsub("\\..*", "", fdata.GTEx$gene_id)
tb$gene_name = fdata.GTEx$gene_name[match(tb$gene_id,gene_id)]

head(tb)

# Load KEGG pathways
pathways <- gmtPathways("/home/ubuntu/c2.cp.kegg.v7.1.symbols.gmt")
# pathways <- gmtPathways("/home/ubuntu/c5.go.bp.v2022.1.Hs.symbols.gmt")
# pathways <- gmtPathways("/home/ubuntu/c2.cp.reactome.v2022.1.Hs.symbols.gmt")

# Number of significant genes
ESR1_genes = tb$gene_name[which(tb$adj.P.Val<0.05 & tb$TF == "ESR1")]
ESR2_genes = tb$gene_name[which(tb$adj.P.Val<0.05 & tb$TF == "ESR2")]
AR_genes = tb$gene_name[which(tb$adj.P.Val<0.05 & tb$TF == "AR")]

# Number of significant aging-related genes that are also disease-relevant
age_table = fread("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/limma_GTEx_aging.txt")
age_table = data.frame(age_table)
aging_genes = age_table$gene_name[which(age_table$P.Value<0.05)]

length(intersect(aging_genes, ESR1_genes))
length(intersect(aging_genes, ESR2_genes))
length(intersect(aging_genes, AR_genes))

# Create network object with only significant edges and connected to aging genes
tb0 = tb[which(tb$TF %in% c("ESR1", "ESR2", "AR")),]
tb1 = tb0[which(tb0$adj.P.Val<0.05 & tb0$gene_name %in% aging_genes),]
diffNet0 = data.frame(tb0$t)
diffNet1 = data.frame(tb1$t)
colnames(diffNet0) = "force"
diffNet0$from = tb0$TF
diffNet0$to = tb0$gene_name
diffNet0$arrows = "to" 
head(diffNet0)

# Select only top 20% edges
diffNet = diffNet0[1:floor(nrow(tb0)*0.2),]

# Plot graph
library(visNetwork)

# Network Plot
nodes <- data.frame(id = unique(as.vector(unlist(diffNet[,c(2,3)]))) , 
                    label=unique(as.vector(unlist(diffNet[,c(2,3)]))))
nodes$group <- ifelse(nodes$id %in% diffNet$from, "TF", "gene")

# Color male edges blue and female edges red
edge_col = rep("purple", nrow(diffNet))
edge_col[which(sign(diffNet[,1]) == -1)] = "green"

diffNet$color = edge_col

net <- visNetwork(nodes, diffNet, width = "100%", 
                  main = "Most Differential Edges in HT") 
net <- visGroups(net, groupname = "TF", shape = "square",
                 color = list(background = "teal", border="black"))
net <- visGroups(net, groupname = "gene", shape = "dot",       
                 color = list(background = "gold", border="black"))
net
visLegend(net, main="Legend", position="right", ncol=1) 

# Save plot as PNG
html_name <- tempfile(fileext = ".html")
visSave(net, html_name)

library(webshot); #webshot::install_phantomjs() #in case phantomjs was not installed 

webshot(html_name, zoom = 2, file = "/home/esaha/Aging_thyroid/results/plots/network_plot/HT_ESR_AR_subnetwork.png")

##############################
# Color genes by KEGG pathway
# Load aging fgsea results
aging_fgsea = get(load("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/GTEx_Thyroid_sexDifference_GSEA.RData"))
sig_paths_age = aging_fgsea$pathway[which(aging_fgsea$padj<0.05)]
sig_genes = c(ESR1_genes, ESR2_genes, AR_genes)
# Pathways associated to differential genes
pathlist = lapply(pathways, function(x){intersect(x,sig_genes)})
pathlength = lapply(pathlist, length)

# pathways that are also significant in aging
# pathlist = pathlist[which(names(pathlist) %in% sig_paths_age)]

# Select only top 10% edges
diffNet = diffNet0[1:floor(nrow(diffNet0)*0.10),]

# Plot graph
library(visNetwork)

# Network Plot
nodes <- data.frame(id = unique(as.vector(unlist(diffNet[,c(2,3)]))) , 
                    label=unique(as.vector(unlist(diffNet[,c(2,3)]))))
nodes$group <- ifelse(nodes$id %in% diffNet$from, "TF", "gene")

node_genes = nodes$id[which(nodes$group == "gene")]
pathtoplot = node_genes
for (i in 1:length(node_genes)){
  g = node_genes[i]
  y = unlist(lapply(pathlist, function(x){g %in% x}))
  if (sum(as.numeric(y)) > 0){
  pathtoplot[i] = names(which.max(pathlength[y]))
  }else {pathtoplot[i] == ""}
}
pathtoplot = stringr::str_to_title(lapply(strsplit(pathtoplot, split = "_"), function(x){paste(x[-1], collapse = " ")}))
nodes$group[which(nodes$group == "gene")] = pathtoplot

# Color male edges blue and female edges red
edge_col = rep("purple", nrow(diffNet))
edge_col[which(sign(diffNet[,1]) == -1)] = "green"

diffNet$color = edge_col

# Remove genes that are not in any pathway
nodes0 = nodes
nodes = nodes0[-which(nodes$group == ""),]

# Adjust node size for better visibility
nodes$size <- rep(25, nrow(nodes))  # Increase node size
nodes$font.size <- rep(25, nrow(nodes))  # Make labels more readable

# Make gene groups square and different color by pathway
gene_groups = unique(nodes$group)
# Generate distinct colors for each group
set.seed(123)  # For reproducibility
library(RColorBrewer)
# Generate color blind friendly colors using the "Dark2" palette
group_colors <- colorRampPalette(brewer.pal(n=9, name = "Set1"))(length(gene_groups))

net <- visNetwork(nodes, diffNet, width = "100%", 
                  main = "Most Differential Edges in HT") 

for (i in 1:length(gene_groups)) {
  net <- visGroups(net, groupname = gene_groups[i], 
                   shape = "square", 
                   color = list(background = group_colors[i], border = "black"))
}
net <- visGroups(net, groupname = "TF", shape = "triangle",
                 color = list(background = "blue", border="black"))

# Adjust physics settings to spread nodes out
net <- visPhysics(net, solver = "forceAtlas2Based", 
                  forceAtlas2Based = list(gravity = -40)) %>% 
  visLegend(stepY=10) # Negative gravity spreads nodes

# net <- visLegend(net, main="Legend", position="right", ncol=1) 
net
# Save plot as PNG
html_name <- tempfile(fileext = ".html")
visSave(net, html_name)

library(webshot); #webshot::install_phantomjs() #in case phantomjs was not installed 

webshot(html_name, zoom = 1, file = "/home/esaha/Aging_thyroid/results/plots/network_plot/HT_ESR_AR_subnetwork_pathColored.png")


