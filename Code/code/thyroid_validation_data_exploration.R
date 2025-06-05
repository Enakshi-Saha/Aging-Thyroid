# Load the GEOquery R/Bioconductor package:
library(GEOquery)

########## Dataset 1: GSE165724
#To access GEO Dataset (GDS), use the function getGEO() which returns a list of ExpressionSets:
gse <- getGEO("GSE165724", GSEMatrix = TRUE)
show(gse)

save(gse, file = "/home/esaha/Aging_thyroid/data/validation_data/GSE165724/GSE165724.RData")

# phenotypic data
phenotypes = pData(gse$GSE165724_series_matrix.txt.gz)
tissue_type = phenotypes$`tissue type:ch1`
gender = phenotypes$`gender:ch1`

table(tissue_type, gender)

########## Dataset 2: GSE138198
#To access GEO Dataset (GDS), use the function getGEO() which returns a list of ExpressionSets:
gse <- getGEO("GSE138198", GSEMatrix = TRUE)
show(gse)

save(gse, file = "/home/esaha/Aging_thyroid/data/validation_data/GSE138198/GSE138198.RData")

# Feature data
GPL6244 <- getGEO("GPL6244", GSEMatrix = TRUE)
show(GPL6244)

GPL6244 = Table(GPL6244)
fdata.GSE138198 = GPL6244[,c("ID", "seqname", "gene_assignment")]
fdata.GSE138198$gene_assignment = unlist(lapply(strsplit(fdata.GSE138198$gene_assignment, split=" // "), function(x){x[2]}))

write.table(fdata.GSE138198, "/home/esaha/Aging_thyroid/data/validation_data/GSE138198/GSE138198_featuredata.txt", row.names = T, col.names = T)

# phenotypic data
phenotypes = pData(gse$GSE138198_series_matrix.txt.gz)
tissue_type = phenotypes$`sample type:ch1`
gender = phenotypes$`gender:ch1`

tissue_labels = levels(factor(tissue_type))
tissue_type[which(tissue_type == tissue_labels[1])] = "HT"
tissue_type[which(tissue_type == tissue_labels[2])] = "mPTC"
tissue_type[which(tissue_type == tissue_labels[3])] = "TN"
tissue_type[which(tissue_type == tissue_labels[4])] = "PTC-HT"
tissue_type[which(tissue_type == tissue_labels[5])] = "PTC"

table(tissue_type, gender)

# Get gene annotations: GPL6244
gpl6244 <- getGEO("GSE165724", GSEMatrix = TRUE)
show(gpl6244)

########## Dataset 3: GSE29315
#To access GEO Dataset (GDS), use the function getGEO() which returns a list of ExpressionSets:
gse <- getGEO("GSE29315", GSEMatrix = TRUE)
show(gse)

save(gse, file = "/home/esaha/Aging_thyroid/data/validation_data/GSE29315/GSE29315.RData")

# phenotypic data
phenotypes = pData(gse$GSE29315_series_matrix.txt.gz)
tissue_type = phenotypes$`clinical description:ch1`

table(tissue_type)

########## Dataset 4: GSE65144
#To access GEO Dataset (GDS), use the function getGEO() which returns a list of ExpressionSets:
gse <- getGEO("GSE65144", GSEMatrix = TRUE)
show(gse)

save(gse, file = "/home/esaha/Aging_thyroid/data/validation_data/GSE65144/GSE65144.RData")

# phenotypic data
phenotypes = pData(gse$GSE65144_series_matrix.txt.gz)
tissue_type = phenotypes$`tissue type:ch1`

table(tissue_type)

########## Dataset 5: GSE33630
#To access GEO Dataset (GDS), use the function getGEO() which returns a list of ExpressionSets:
gse <- getGEO("GSE33630", GSEMatrix = TRUE)
show(gse)

save(gse, file = "/home/esaha/Aging_thyroid/data/validation_data/GSE33630/GSE33630.RData")

# phenotypic data
phenotypes = pData(gse$GSE33630_series_matrix.txt.gz)
tissue_type = phenotypes$`pathological diagnostic:ch1`

table(tissue_type)

########## Dataset 6: GSE9340
#To access GEO Dataset (GDS), use the function getGEO() which returns a list of ExpressionSets:
gse <- getGEO("GSE9340", GSEMatrix = TRUE)
show(gse)

save(gse, file = "/home/esaha/Aging_thyroid/data/validation_data/GSE9340/GSE9340.RData")

# phenotypic data
phenotypes = pData(gse$GSE9340_series_matrix.txt.gz)
tissue_type = phenotypes$`DISEASE STATE:ch1`
gender = phenotypes$`SEX:ch1`

table(tissue_type, gender)

########## Dataset 7: GSE27155
#To access GEO Dataset (GDS), use the function getGEO() which returns a list of ExpressionSets:
gse <- getGEO("GSE27155", GSEMatrix = TRUE)
show(gse)

save(gse, file = "/home/esaha/Aging_thyroid/data/validation_data/GSE27155/GSE27155.RData")

# phenotypic data
phenotypes = pData(gse$GSE27155_series_matrix.txt.gz)
tissue_type = phenotypes$`tissue:ch1`

table(tissue_type)
