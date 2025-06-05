################################################
# Color genes by KEGG pathway
# Pathways associated to differential genes
pathlist = lapply(pathways[which(names(pathways) %in% sig_pathsAB)], function(x){intersect(x,sig_genes)})
pathlist[which(names(pathlist) %in% sig_path_limma)]
pathlength = lapply(pathlist, length)

# Network Plot
nodes <- data.frame(id = unique(as.vector(as.matrix(diffNet[,c(2,3)]))) , 
                    label=unique(as.vector(as.matrix(diffNet[,c(2,3)]))))
nodes$group <- ifelse(nodes$id %in% diffNet$from, "miRNA", "gene")
node_genes = nodes$id[which(nodes$group == "gene")]

pathtoplot = node_genes
for (i in 1:length(node_genes)){
  g = node_genes[i]
  pathtoplot[i] = names(which.max(pathlength[unlist(lapply(pathlist, function(x){g %in% x}))]))
}
pathtoplot = stringr::str_to_title(lapply(strsplit(pathtoplot, split = "_"), function(x){paste(x[-1], collapse = " ")}))
nodes$group[which(nodes$group == "gene")] = pathtoplot

# Color male edges blue and female edges red
edge_col = rep("red", nrow(diffNet))
edge_col[which(sign(diffNet[,1]) == -1)] = "blue"

diffNet$color = edge_col

net <- visNetwork(nodes, diffNet, width = "100%") 
net <- visGroups(net, groupname = "miRNA", shape = "square",
                 color = list(background = "teal", border="black"))
net
#visLegend(net, main="Legend", position="right", ncol=1) 

# Save plot as PNG
html_name <- tempfile(fileext = ".html")
visSave(net, html_name)

library(webshot); #webshot::install_phantomjs() #in case phantomjs was not installed 

webshot(html_name, zoom = 5, file = "/home/esaha/BONOBO/plots/mRNA_miRNA_breastCancer/LumABsigDiff_network_pathcolor.png")

cbind(node_genes, pathtoplot)

