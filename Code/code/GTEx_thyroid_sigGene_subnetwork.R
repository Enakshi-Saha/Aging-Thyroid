##### Vectorize lioness matrices of individual samples and combine into a single matrix ######
# Get list of significant genes
# Load limma tables
tb_male <- read.csv("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_aging_male.txt", sep="")
tb_female <- read.csv("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_aging_female.txt", sep="")

m_i = rownames(tb_male)[which((tb_male$logFC >= 1) & (tb_male$P.Value < 0.05))]
m_d = rownames(tb_male)[which((tb_male$logFC <= -1) & (tb_male$P.Value < 0.05))]

f_i = rownames(tb_female)[which((tb_female$logFC >= 1) & (tb_female$P.Value < 0.05))]
f_d = rownames(tb_female)[which((tb_female$logFC <= -1) & (tb_female$P.Value < 0.05))]

sigGenes = unique(c(m_i, m_d, f_i, f_d))

library(data.table)
setwd("/home/ubuntu/BonoboPanda_GTEx_thyroid/single_panda/")
files_ls = list.files(path=".", pattern=".txt", all.files=TRUE,
                      full.names=TRUE)

indices = unlist(lapply(gsub("\\/*","",files_ls), function(x){unlist(strsplit(x, split = ".", fixed = T))[2]}))
data = list()
i=1
new_data = fread(files_ls[1])
new_data = new_data[which(new_data$gene %in% sigGenes),]
TFs = new_data$tf
genes = new_data$gene

for (file_name in files_ls){
  new_data = fread(file_name)
  data[[i]] = new_data[which(new_data$gene %in% sigGenes),4]
  cat(paste(i,"\n"))
  i = i+1
}

outdegree = data.frame(do.call(cbind, data))

colnames(outdegree) = indices
outdegree$TF = TFs
outdegree$gene = genes

head(outdegree)
dim(outdegree)

setwd("/home/esaha/Aging_thyroid/data/network/")
write.table(outdegree, file="BonoboPanda_GTEx_thyroid_significantGene_subnetwork.txt", row.names=T, col.names=T)
