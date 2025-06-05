library(ggplot2)
library(ggrepel)

# Load limma tables
tb <- read.csv("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/limma_GTEx_aging.txt", sep="")

# cutoff for logFC
FClim = 1
# add a column of NAs
tb$direction <- "Not significant"
# if log2Foldchange > 25 and -log10(pvalue) > 5, set as "UP" 
tb$direction[tb$logFC > FClim & tb$P.Value < 0.05] <- "increasing"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
tb$direction[tb$logFC < -FClim & tb$P.Value < 0.05] <- "decreasing"

mycolors <- c("red", "blue", "grey")
names(mycolors) <- c("increasing", "decreasing", "Not significant")

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
tb$delabel <- NA
tb$delabel[tb$direction != "Not significant"] <- tb$gene_name[tb$direction != "Not significant"]

# Remove gene names that copntain "-" or "." in name
tb$delabel[grep("-",tb$delabel)] <- NA
tb$delabel[grep(".",tb$delabel, fixed = T)] <- NA

# Re-plot but this time color the points with "diffexpressed"
ggplot(data=tb, aes(x=logFC, y=-log10(P.Value), col=direction, label = delabel)) +
  geom_point() + theme_minimal() +
  geom_vline(xintercept=c(-FClim, FClim), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  scale_colour_manual(values = mycolors) +
  ggtitle("Volcanoplot of Differentially Targeting TFs with Age") +
  geom_text_repel(size = 3)
