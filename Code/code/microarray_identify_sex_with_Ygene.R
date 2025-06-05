library(GEOquery)
library(ggplot2)

# Get names of Y genes
GTEx_thyroid = get(load("/home/esaha/BONOBO/data/GTEx_thyroid/GTEx_thyroid_RSE.RData"))
fdata.GTEx = as.data.frame(GTEx_thyroid@rowRanges)
Y_genes = fdata.GTEx$gene_name[which(fdata.GTEx$seqnames == "chrY")]

# Count number of genes and Y genes
result = list()
# pdf for saving PCA plots
datasetIDs = c("GSE29315", "GSE65144", "GSE33630", "GSE27155")

pdf("/home/esaha/Aging_thyroid/results/plots/PCAplots_microarray.pdf")
for (datasetID in datasetIDs){
  paste("running ", datasetID, "\n", sep ="")
  data = get(load(paste("~/Aging_thyroid/data/validation_data/",datasetID,"/",datasetID,".RData", sep = "")))
  dataName = paste(datasetID,"_series_matrix.txt.gz", sep = "")
  expr = t(exprs(data[[dataName]]))
  gene_info = fData(data[[dataName]])
  
  # Count number of unique genes and Y genes
  ngene = length(unique(gene_info$`Gene Symbol`))
  nYgene = length(which(gene_info$`Gene Symbol` %in% Y_genes))
  Ygene = rownames(gene_info)[which(gene_info$`Gene Symbol` %in% Y_genes)]
  nsample = nrow(expr)
  result[[datasetID]] = c(ngene, nYgene, nsample)
  
  # PCA of Y genes to identify sex
  yexpr = expr[,which(colnames(expr) %in% Ygene)]
  # Remove genes that are zero everywhere
  zerogenes = length(which(colSums(yexpr) == 0))
  if (zerogenes > 0){
    yexpr = yexpr[,-which(colSums(yexpr) == 0)]
  }
  pca <- prcomp(yexpr)
  U <- pca$x
  theme_set(bigstatsr::theme_bigstatsr(0.8))
  p1 = qplot(U[, 1], U[, 2]) + coord_equal() + 
    ggtitle(paste("PCA of Y genes:",datasetID,sep = " ")) + xlab("PC 1") + ylab ("PC 2")
  print(p1)
  
  # PCA plot colored by tissue type
  phenotypes = pData(data[[dataName]])
  tissue_type = phenotypes[,ncol(phenotypes)]
  p2 = qplot(U[, 1], U[, 2], col = tissue_type) + coord_equal() + 
    ggtitle(paste("PCA of Y genes:",datasetID,sep = " ")) + xlab("PC 1") + ylab ("PC 2")
  print(p2)
  
  # PCA plot of Y genes and XIST
  gene_names = gene_info$`Gene Symbol`[match(colnames(expr), rownames(gene_info))]
  xist = expr[,which(gene_names == "XIST")]
  if (length(xist) != 0){
    if (class(xist)[1] == "matrix"){
      xist = rowMeans(xist)
    }
    # inferred sex
    p3 = qplot(U[, 1], U[, 2], col = xist) + coord_equal() + 
      ggtitle(paste("PCA of Y genes (colored by XIST)",datasetID,sep = " ")) + xlab("PC1") + ylab ("PC2")
    print(p3)
  }
  # Identify sex
  sex = rep("NA",length(tissue_type))
  if (datasetID == "GSE29315"){
    sex[which(U[,1]<=-2.5)] = "M"
    sex[which(U[,1]>=-1)] = "F"
    p_GSE29315 = data.frame(cbind(sex, tissue_type))
    rownames(p_GSE29315) = rownames(phenotypes)
    write.table(p_GSE29315, paste("~/Aging_thyroid/data/validation_data/",datasetID,"/",datasetID,"_phenotypes.txt", sep = ""),
               col.names = T, row.names = T)
  }
  sex = rep("NA",length(tissue_type))
  if (datasetID == "GSE33630"){
    sex[which(U[,1]>=5)] = "M"
    sex[which(U[,1]<=5)] = "F"
    tissue_type =  phenotypes[,(ncol(phenotypes)-1)]
    p_GSE33630 = data.frame(cbind(sex, tissue_type))
    rownames(p_GSE33630) = rownames(phenotypes)
    write.table(p_GSE33630, paste("~/Aging_thyroid/data/validation_data/",datasetID,"/",datasetID,"_phenotypes.txt", sep = ""),
                col.names = T, row.names = T)
  }
  sex = rep("NA",length(tissue_type))
  if (datasetID == "GSE27155"){
    sex[which(U[,1]>=1)] = "M"
    sex[which(U[,1]<=0)] = "F"
    p_GSE27155 = data.frame(cbind(sex, tissue_type))
    rownames(p_GSE27155) = rownames(phenotypes)
    write.table(p_GSE27155, paste("~/Aging_thyroid/data/validation_data/",datasetID,"/",datasetID,"_phenotypes.txt", sep = ""),
                col.names = T, row.names = T)
  }
}
dev.off()
result



