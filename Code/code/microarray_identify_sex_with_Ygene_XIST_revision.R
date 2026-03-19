# Modified PCA plots colored by XIST expression,
# with "included in study" samples marked by a different symbol
library(GEOquery)
library(ggplot2)
library(recount3)

# Get names of Y genes
GTEx_thyroid = get(load("/home/esaha/thyroid_revision_data/restored_data/GTEx_thyroid_RSE.RData"))
fdata.GTEx = as.data.frame(GTEx_thyroid@rowRanges)
Y_genes = fdata.GTEx$gene_name[which(fdata.GTEx$seqnames == "chrY")]

# Count number of genes and Y genes
result = list()
# pdf for saving PCA plots
datasetIDs = c("GSE29315", "GSE65144", "GSE33630", "GSE27155")

library(ggplot2)

# Define which tissue types are "included in study" per dataset
included_tissues <- list(
  GSE29315 = c("hashimoto thyroiditis","hyperplasia"),
  GSE65144 = "anaplastic thyroid carcinoma (ATC)",
  GSE33630 = c("anaplastic thyroid carcinoma (ATC)","patient-matched non-tumor control"),
  GSE27155 = c("Anaplastic Thyroid Carcinoma", "Normal Thyroid")
)

pdf("/home/esaha/thyroid_revision_results/revision_plots/PCAplots_microarray_marked.pdf")

for (datasetID in datasetIDs){
  
  data = get(load(paste("/home/esaha/thyroid_revision_data/validation_data/",datasetID,"/",datasetID,".RData", sep = "")))
  dataName = paste(datasetID,"_series_matrix.txt.gz", sep = "")
  expr = t(exprs(data[[dataName]]))
  gene_info = fData(data[[dataName]])
  
  # Y gene PCA (same as before)
  Ygene = rownames(gene_info)[which(gene_info$`Gene Symbol` %in% Y_genes)]
  yexpr = expr[, which(colnames(expr) %in% Ygene)]
  zerogenes = which(colSums(yexpr) == 0)
  if (length(zerogenes) > 0){
    yexpr = yexpr[, -zerogenes]
  }
  pca <- prcomp(yexpr)
  U <- pca$x
  
  # Get XIST expression (same as before)
  gene_names = gene_info$`Gene Symbol`[match(colnames(expr), rownames(gene_info))]
  xist = expr[, which(gene_names == "XIST")]
  
  if (length(xist) != 0){
    if (class(xist)[1] == "matrix"){
      xist = rowMeans(xist)
    }
    
    # Get tissue types for this dataset
    phenotypes = pData(data[[dataName]])
    tissue_type = phenotypes[, ncol(phenotypes)]
    
    # Special case: GSE33630 uses second-to-last column for tissue type
    if (datasetID == "GSE33630"){
      tissue_type = phenotypes[, (ncol(phenotypes) - 1)]
    }
    
    # Create "included in study" grouping variable
    included_label <- included_tissues[[datasetID]]
    in_study <- ifelse(
      tolower(trimws(tissue_type)) %in% tolower(trimws(included_label)),
      "included in study",
      "other"
    )
    
    # Build plot dataframe
    plot_df <- data.frame(
      PC1   = U[, 1],
      PC2   = U[, 2],
      xist  = xist,
      group = in_study
    )
    
    # Plot: color = XIST, shape = included in study vs other
    # Plot: color = XIST, shape = included in study vs other
    p <- ggplot(plot_df, aes(x = PC1, y = PC2)) +
      # Red border layer for "included in study" points only
      geom_point(
        data  = subset(plot_df, group == "included in study"),
        aes(x = PC1, y = PC2),
        shape = 21, size = 3, color = "red", fill = NA, stroke = 1.2
      ) +
      geom_point(aes(color = xist, shape = group, size = group)) +
      coord_equal() +
      scale_color_continuous(name = "XIST") +
      scale_shape_manual(
        name   = "Sample group",
        values = c("included in study" = 17, "other" = 16),
        labels = c("included in study" = "included in study", "other" = "other")
      ) +
      scale_size_manual(
        values = c("included in study" = 3, "other" = 1),
        guide  = "none"   # hides the redundant size legend
      ) +
      ggtitle(paste("PCA of Y genes (colored by XIST)", datasetID, sep = " ")) +
      xlab("PC1") +
      ylab("PC2") +
      theme_bw()
    
    print(p)
  }
}

dev.off()
