library(GEOquery)
library(ggplot2)
library(data.table)
library(dplyr)

# Load microarray data
data = get(load(paste("~/Aging_thyroid/data/validation_data/GSE33630/GSE33630.RData", sep = "")))
dataName = "GSE33630_series_matrix.txt.gz"
expr = exprs(data[[dataName]])
gene_info = fData(data[[dataName]])

# Get gene names
genes = unlist(lapply(strsplit(gene_info$`Gene Symbol`, split="///"), function(x){x[1]}))
gene = genes[match(rownames(expr),rownames(gene_info))]

# For genes with multiple transcripts, keep transcripts with the highest variability
transcript_sd = unlist(apply(expr, MARGIN = 1, sd))
df = data.frame(cbind(gene, transcript_sd, names(transcript_sd)))

filtered = df %>% group_by(gene) %>% dplyr::slice(which.max(transcript_sd))
genes_filtered = filtered$gene
transcripts_filtered = filtered$V3
expression = expr[match(transcripts_filtered, rownames(expr)),]
rownames(expression) = genes_filtered
expression = data.frame(expression)

# log transform expression
expression = log2(1+expression)
head(expression[,1:5])

#######################################################################
# Phenotypic data processed on the "microarray_identify_sex_with_Ygene.R"
phenotypes <- read.csv("~/Aging_thyroid/data/validation_data/GSE33630/GSE33630_phenotypes.txt", sep="")

# Prepare motif file to run BONOBO-PANDA
sample = rownames(phenotypes)
female_prior = "/home/ubuntu/prior_large/prior_female_HGNC.txt"
male_prior = "/home/ubuntu/prior_large/prior_male_HGNC.txt"
prior = rep(0, nrow(phenotypes))
prior[which(phenotypes$sex == "M")] = male_prior
prior[which(phenotypes$sex == "F")] = female_prior
prior_tab = data.frame(cbind(sample, prior))

expression = cbind(rownames(expression), expression)
colnames(expression)[1] = "gene"

head(expression[,1:5])
head(prior_tab)

# Save expression data
write.table(expression, "~/Aging_thyroid/data/validation_data/GSE33630/GSE33630_expression.txt",
            row.names = F, col.names = T, sep = '\t', quote = F)
# save motif prior table
write.table(prior_tab, paste("~/Aging_thyroid/data/validation_data/GSE33630/GSE33630_motif.csv", sep = ""), row.names = FALSE, col.names = TRUE, sep = ",", quote = F)


