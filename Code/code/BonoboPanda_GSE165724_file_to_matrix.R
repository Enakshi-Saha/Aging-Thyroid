library(data.table)
setwd("/home/ubuntu/BonoboPanda_GSE165724/single_panda/")
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
write.table(lioness, file="BonoboPanda_GSE165724.txt", row.names=T, col.names=T)
