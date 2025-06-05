library(fgsea)

# Create Pathway super group gene list

# Load KEGG pathways
pathways <- gmtPathways("/home/ubuntu/c2.cp.kegg.v7.1.symbols.gmt")

# Load KEGG pathway hierarchy
pathway_hierarchy = get(load("/home/esaha/Aging_thyroid/data/KEGG_pathway_hierarchy_byFunction.RData"))

pathway_superclass = c("metabolism",
                       "genetic information processing",
                       "cell signaling, cell growth and death",
                       "cell communication",
                       "immune system and immune diseases",
                       "endocrine, circulatory, nervous system and related diseases",
                       "metabolic diseases",
                       "infectious diseases")

pathway_superclass_index = c("c(2,4)", "c(14,15,17)", "c(20,23)", "c(21, 24)", "c(25, 34)", "c(26, 27, 29, 33, 35)", "37", "38")
path_superclass = cbind(pathway_superclass, pathway_superclass_index)

pathway_superclass_genes = list()
for (i in 1:length(pathway_superclass)){
  pathway_class = pathway_superclass[i]
  
  j = eval(parse(text = path_superclass[i,2]))
  pathList = unlist(pathway_hierarchy[j])
  
  pathgenes = unname(unlist(c(pathways[which(names(pathways) %in% pathList)])))
  
  pathway_superclass_genes[[pathway_class]] = pathgenes
}


pathway_superclass_genes

save(pathway_superclass_genes, file = "/home/esaha/Aging_thyroid/results/KEGG_pathway_superClass_geneList.RData")

