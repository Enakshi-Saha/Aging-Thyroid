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

# KEGG pathway superclass
pathways = get(load("/home/esaha/Aging_thyroid/results/KEGG_pathway_superClass_geneList.RData"))

# Get ESR-AR edges for all samples, males, females
ESR1_table = data.frame(fread(paste("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/ESR_AR/limma_GTEx_aging_ESR1.txt", sep = "")))
ESR1_male_table = data.frame(fread(paste("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/ESR_AR/limma_GTEx_M_aging_ESR1.txt", sep = "")))
ESR1_female_table = data.frame(fread(paste("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/ESR_AR/limma_GTEx_F_aging_ESR1.txt", sep = "")))

ESR2_table = data.frame(fread(paste("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/ESR_AR/limma_GTEx_aging_ESR2.txt", sep = "")))
ESR2_male_table = data.frame(fread(paste("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/ESR_AR/limma_GTEx_M_aging_ESR2.txt", sep = "")))
ESR2_female_table = data.frame(fread(paste("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/ESR_AR/limma_GTEx_F_aging_ESR2.txt", sep = "")))

AR_table = data.frame(fread(paste("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/ESR_AR/limma_GTEx_aging_AR.txt", sep = "")))
AR_male_table = data.frame(fread(paste("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/ESR_AR/limma_GTEx_M_aging_AR.txt", sep = "")))
AR_female_table = data.frame(fread(paste("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/ESR_AR/limma_GTEx_F_aging_AR.txt", sep = "")))

TFs = c(rep("ESR1", nrow(ESR1_table)), rep("ESR2", nrow(ESR2_table)), rep("AR", nrow(AR_table)))
all_sample_table = rbind(ESR1_table, ESR2_table, AR_table)
all_sample_table$TF = TFs
male_table = rbind(ESR1_male_table, ESR2_male_table, AR_male_table)
male_table$TF = TFs
female_table = rbind(ESR1_female_table, ESR2_female_table, AR_female_table)
female_table$TF = TFs

all_sample_table = all_sample_table[order(all_sample_table$adj.P.Val),]
male_table = male_table[order(all_sample_table$adj.P.Val),]
female_table = female_table[order(all_sample_table$adj.P.Val),]

# Top aging genes
aging_genes = all_sample_table$gene_name[which(all_sample_table$P.Value<0.05)]

# Create network object with only significant edges and connected to aging genes
# tb0 = male_table[which(male_table$gene_name %in% aging_genes),]
tb0 = female_table[which(female_table$gene_name %in% aging_genes),]

diffNet0 = data.frame(tb0$t)
colnames(diffNet0) = "force"
diffNet0$from = tb0$TF
diffNet0$to = tb0$gene_name
diffNet0$arrows = "to" 
diffNet0$pval = tb0$P.Value
head(diffNet0)

# Select only top 20% edges
diffNet = diffNet0

# Plot graph
library(visNetwork)

# webshot(html_name, zoom = 2, file = "/home/esaha/Aging_thyroid/results/plots/network_plot/HT_ESR_AR_subnetwork.png")

##############################
# Color genes by KEGG pathway
# Load aging fgsea results
aging_fgsea = get(load("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/GTEx_Thyroid_sexDifference_GSEA.RData"))
sig_paths_age = aging_fgsea$pathway[which(aging_fgsea$padj<0.05)]
sig_genes = aging_genes
# Pathways associated to differential genes
pathlist = lapply(pathways, function(x){intersect(x,sig_genes)})
pathlength = lapply(pathlist, length)

# pathways that are also significant in aging
# pathlist = pathlist[which(names(pathlist) %in% sig_paths_age)]

# Plot graph
library(visNetwork)

# Network Plot
nodes <- data.frame(id = unique(as.vector(unlist(diffNet[,c(2,3)]))) , 
                    label=unique(as.vector(unlist(diffNet[,c(2,3)]))))
nodes$group <- ifelse(nodes$id %in% diffNet$from, "TF", "gene")

node_genes = nodes$id[which(nodes$group == "gene")]
pathtoplot = rep("", length(node_genes))
for (i in 1:length(node_genes)){
  g = node_genes[i]
  y = unlist(lapply(pathlist, function(x){g %in% x}))
  if (sum(as.numeric(y)) > 0){
    pathtoplot[i] = names(which.max(pathlength[y]))
  }else {pathtoplot[i] == ""}
}
# pathtoplot = stringr::str_to_title(lapply(strsplit(pathtoplot, split = "_"), function(x){paste(x[-1], collapse = " ")}))
nodes$group[which(nodes$group == "gene")] = pathtoplot

# Color non-significant edges gray, significantly positive edges red, significantly negative edges blue
edge_col = rep("gray", nrow(diffNet))
edge_col[which(sign(diffNet$force) == -1 & diffNet$pval < 0.05)] = "blue"
edge_col[which(sign(diffNet$force) == 1 & diffNet$pval < 0.05)] = "red"

diffNet$color = edge_col

# Remove genes that are not in any pathway
nodes0 = nodes
if ("" %in% nodes$group){
  nodes = nodes0[-which(nodes$group == ""),]
}

# Adjust node size for better visibility
nodes$size <- rep(25, nrow(nodes))  # Increase node size
nodes$font.size <- rep(25, nrow(nodes))  # Make labels more readable

# Make gene groups square and different color by pathway
gene_groups = unique(nodes$group)
# Generate distinct colors for each group
set.seed(123)  # For reproducibility
library(RColorBrewer)
# Generate color blind friendly colors using the "Dark2" palette
# group_colors <- colorRampPalette(brewer.pal(n=9, name = "Set1"))(length(gene_groups)) # 6=length(gene_groups)
group_colors = c("white", "skyblue", "green", "yellow", "orange", "pink")

net <- visNetwork(nodes, diffNet, width = "100%", 
                  main = "GTEx females: top aging-associated edges connected to ESR1, ESR2, AR") 

gene_groups = c("TF", "cell signaling, cell growth and death",
                "endocrine, circulatory, nervous system and related diseases",
                "immune system and immune diseases",                          
                "metabolism",                                                 
                "cell communication")
for (i in 1:length(gene_groups)) {
  net <- visGroups(net, groupname = gene_groups[i], 
                   shape = "square", 
                   color = list(background = group_colors[i], border = "black"))
}
net <- visGroups(net, groupname = "TF", shape = "triangle",
                 color = list(background = "white", border="black"))

# Adjust physics settings to spread nodes out
net <- visPhysics(net, solver = "forceAtlas2Based", 
                  forceAtlas2Based = list(gravity = -20)) 
# %>% visLegend(stepY=10) # Negative gravity spreads nodes

# net <- visLegend(net, main="Legend", position="right", ncol=1) 
net
# Save plot as PNG
html_name <- tempfile(fileext = ".html")
visSave(net, html_name)

library(webshot)

# webshot::install_phantomjs() #in case phantomjs was not installed 

# webshot(html_name, zoom = 1, file = "/home/esaha/Aging_thyroid/results/plots/network_plot/GTEx_male_ESR_AR_subnetwork_pathColored.png")
# webshot(html_name, zoom = 1, file = "/home/esaha/Aging_thyroid/results/plots/network_plot/GTEx_female_ESR_AR_subnetwork_pathColored.png")

# Create a clean legend
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend = gene_groups[-1], pch=16, pt.cex=3, cex=1.5, bty='n',
       col = group_colors[-1])
mtext("Pathway Group", at=0.2, cex=2)
