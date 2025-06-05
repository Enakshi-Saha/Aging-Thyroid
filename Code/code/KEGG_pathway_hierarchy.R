library(keggorthology)
library(fgsea)

# KEGG pathways
pathways <- gmtPathways("/home/ubuntu/c2.cp.kegg.v7.1.symbols.gmt")

# KEGG hierarchy graph
data(KOgraph)

# All nodes
nodes(KOgraph)

# KEGG pathway hierarchy
path_hierarchy = adj(KOgraph, nodes(KOgraph))
level1 = adj(KOgraph, nodes(KOgraph)[1])$KO.Feb10root

# Remove level 1 nodes
path_hierarchy = path_hierarchy[-which(names(path_hierarchy) %in% level1)]

# Remove single pathways, i.e. list elements with length 0
node_degree = lapply(path_hierarchy, length)

filtered_hierarchy = path_hierarchy[-which(node_degree == 0)]

# Create list of KEGG pathways with superclassses they belong to

find_superclass = function(pathname){
  pathname = gsub("_", "", gsub("KEGG_", "", pathname))
  position = unlist(lapply(filtered_hierarchy,function(x){pathname %in% gsub("[()]|/ |-| |'", "", toupper(x))}))
  return(names(filtered_hierarchy)[position])
}

kegg_path_graph = lapply(names(pathways),find_superclass)
names(kegg_path_graph) = names(pathways)

path_graph = list()
for (i in 1:length(filtered_hierarchy)){
  superclass = names(filtered_hierarchy)[i]
  x = filtered_hierarchy[[i]]
  paths = gsub("[()]|/ |-| |'", "", toupper(x))
  pathway_list = lapply(names(pathways),function(pathname){gsub("_", "", gsub("KEGG_", "", pathname))})
  flag = unlist(lapply(pathway_list,function(z){z %in% paths}))
  path_graph[[superclass]] = names(pathways)[flag]
}

save(kegg_path_graph, file = "/home/esaha/Aging_thyroid/data/KEGG_pathway_hierarchy_byPathway.RData")
save(path_graph, file = "/home/esaha/Aging_thyroid/data/KEGG_pathway_hierarchy_byFunction.RData")


