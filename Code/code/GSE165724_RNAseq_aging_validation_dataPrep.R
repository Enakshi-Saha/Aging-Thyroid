library(GEOquery)
library(ggplot2)
library(data.table)
library(dplyr)
library(gghalves)

# Download RNAseq expression data
datasetID = "GSE165724"
url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE165724&format=file&file=GSE165724%5Fcounts%5F74Samples%2Etsv%2Egz"
file_name = paste("~/Aging_thyroid/data/validation_data/",datasetID,"/",datasetID,"_expression.tsv.gz", sep = "")
# download.file(url = url, destfile = file_name)

# Import expression data
expr = data.frame(fread(file_name))
rownames(expr) = expr$Geneid
expr = expr[,-1]

data = get(load(paste("~/Aging_thyroid/data/validation_data/",datasetID,"/",datasetID,".RData", sep = "")))
dataName = paste(datasetID,"_series_matrix.txt.gz", sep = "")

# Count number of unique genes and Y genes
ngene = nrow(expr)
nYgene = length(which(rownames(expr) %in% Y_genes))

# Get normal thyroid and tumor adjacent normal thyroid samples
phenotypes = pData(data[[dataName]])
expr_normal = expr[,which(phenotypes$`tissue type:ch1` != "PTC tumor tissue")]
nsample = ncol(expr_normal)

gender = phenotypes$`gender:ch1`[which(phenotypes$`tissue type:ch1` != "PTC tumor tissue")]
  
# PCA of Y genes to identify sex
yexpr = expr_normal[which(rownames(expr_normal) %in% Y_genes),]
# Remove genes that are zero everywhere
zerogenes = length(which(colSums(yexpr) == 0))
if (zerogenes > 0){
  yexpr = yexpr[,-which(colSums(yexpr) == 0)]
}
pca <- prcomp(t(yexpr))
U <- pca$x
theme_set(bigstatsr::theme_bigstatsr(0.8))
p = qplot(U[, 1], U[, 2], col = gender) + coord_equal() + 
  ggtitle(paste("PCA of Y genes:",datasetID,sep = " ")) + xlab("PC 1") + ylab ("PC 2")
print(p)

######### Save expression and phenotypic data #######################
# Normalize by TPM
bplength = fdata.GTEx$bp_length[match(rownames(expr),fdata.GTEx$gene_id)]
# divide raw counts by length
expr_matrix = as.matrix(expr_normal)
expression <- expr_matrix / bplength
expression = na.omit(expression)
expression <- t(t(expression) * 1e6 / colSums(expression))

# For genes with multiple transcripts, keep transcripts with the highest variability
gene = fdata.GTEx$gene_name[match(rownames(expression),fdata.GTEx$gene_id)]
transcript_sd = unlist(apply(expression, MARGIN = 1, sd))
df = data.frame(cbind(gene, transcript_sd, names(transcript_sd)))

filtered = df %>% group_by(gene) %>% dplyr::slice(which.max(transcript_sd))
genes_filtered = filtered$gene
transcripts_filtered = filtered$V3
expression = expression[which(rownames(expression) %in% transcripts_filtered),]
rownames(expression) = genes_filtered[match(rownames(expression), transcripts_filtered)]
expression = data.frame(expression)

##### Remove genes with counts < 1 TPM in at least 10% samples
expression_cutoff = 1 # 1 TPM
percent_cutoff = 0.1
minSamples = percent_cutoff*ncol(expression) # at least 20% of samples
keep = rowSums(expression > expression_cutoff) >= minSamples
table(keep)
expression_filtered = expression[keep,]
paste0(nrow(expression), ":Total number of genes overlapped between TCGA and GTEx")
paste0(nrow(expression_filtered), ":Number of genes after filtering ", expression_cutoff," cpm in ", minSamples, " samples")

# Convert to log2(1+TPM)
expression = log2(1+expression_filtered)

# Construct normal-only expression file
colnames(expression) = rownames(phenotypes)[which(phenotypes$`tissue type:ch1` != "PTC tumor tissue")]
expression = cbind(rownames(expression), expression)
colnames(expression)[1] = "gene"
head(expression)[,1:5]

# Save expression data
# write.table(expression, paste("~/Aging_thyroid/data/validation_data/",datasetID,"/",datasetID,"_expression_normalThyroid.txt", sep = ""),
#            row.names = F, col.names = T, sep = '\t', quote = F)

# Prepare phenotypic file
phenotypes_subset = phenotypes[which(phenotypes$`tissue type:ch1` != "PTC tumor tissue"),
                        c("age:ch1","gender:ch1","tissue type:ch1")]
colnames(phenotypes_subset) = c("age", "gender", "tissue_type")
phenotypes_subset$tissue_type[which(phenotypes_subset$tissue_type != "Normal")] = "Normal_adjacent"

head(phenotypes_subset)

# Save phenotypic data
# write.table(phenotypes_subset, paste("~/Aging_thyroid/data/validation_data/",datasetID,"/",datasetID,"_phenotypes_normalThyroid.txt", sep = ""),
#            col.names = T, row.names = T)

# Prepare motif file to run BONOBO-PANDA
sample = rownames(phenotypes_subset)
female_prior = "/home/ubuntu/prior_large/prior_female_HGNC.txt"
male_prior = "/home/ubuntu/prior_large/prior_male_HGNC.txt"
prior = rep(0, nrow(phenotypes_subset))
prior[which(phenotypes_subset$gender == "M")] = male_prior
prior[which(phenotypes_subset$gender == "F")] = female_prior
prior_tab = data.frame(cbind(sample, prior))
head(prior_tab)

# write.table(prior_tab, paste("~/Aging_thyroid/data/validation_data/",datasetID,"/",datasetID,"_normalThyroid_motif.csv", sep = ""), row.names = FALSE, col.names = TRUE, sep = ",", quote = F)

# Plot histogram of age by gender and tissue type
fulldata = phenotypes_subset
fulldata$age = as.numeric(fulldata$age)
fulldata$tissue_type = factor(fulldata$tissue_type, levels = c("Normal_adjacent", "Normal"))
ggplot(fulldata, aes(x = gender, y = age, fill = tissue_type)) + geom_half_violin(side = "l") + geom_half_boxplot(side = "r") +theme_bw() + 
  ggtitle("GSE165724: age distribution by gender and tissue type") + xlab("Sex") + ylab("Age") 


ggplot(fulldata, aes(x = tissue_type, y = age, fill = tissue_type)) + geom_half_violin(side = "l") + geom_half_boxplot(side = "r") +theme_bw() + 
  ggtitle("GSE165724: age distribution by tissue type") + xlab("") + ylab("Age") +
  scale_fill_manual(values=c("purple", "darkgreen"))
