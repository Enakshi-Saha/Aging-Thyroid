library(data.table)

datasetIDs = c("GSE29315", "GSE33630", "GSE27155", "GSE65144")
datasetID = datasetIDs[4]

setwd(paste("/home/ubuntu/BonoboPanda_", datasetID, "/single_panda/", sep = ""))
files_ls = list.files(path=".", pattern=".txt", all.files=TRUE,
                      full.names=TRUE)

indices = unlist(lapply(gsub("\\/*","",files_ls), function(x){unlist(strsplit(x, split = ".", fixed = T))[2]}))
data = list()
i=1
new_data = fread(files_ls[1])
new_data = aggregate(new_data$force, by=list(Category=new_data$gene), FUN=sum)
genes = new_data$Category

for (file_name in files_ls){
  new_data = fread(file_name)
  data[[i]] = aggregate(new_data$force, by=list(Category=new_data$gene), FUN=sum)$x
  cat(paste(i,"\n"))
  i = i+1
}

lioness = do.call(cbind, data)

rownames(lioness) = genes
colnames(lioness) = indices

head(lioness)
dim(lioness)

setwd("/home/esaha/Aging_thyroid/networks/")
write.table(lioness, file=paste("BonoboPanda_", datasetID, ".txt", sep = ""), row.names=T, col.names=T)
