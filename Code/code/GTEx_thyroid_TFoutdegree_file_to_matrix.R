##### Vectorize lioness matrices of individual samples and combine into a single matrix ######
library(data.table)
setwd("/home/ubuntu/BonoboPanda_GTEx_thyroid/single_panda/")
files_ls = list.files(path=".", pattern=".txt", all.files=TRUE,
                      full.names=TRUE)

indices = unlist(lapply(gsub("\\/*","",files_ls), function(x){unlist(strsplit(x, split = ".", fixed = T))[2]}))
data = list()
i=1
new_data = fread(files_ls[1])
new_data = aggregate(new_data$force, by=list(Category=new_data$tf), FUN=sum)
TFs = new_data$Category

for (file_name in files_ls){
  new_data = fread(file_name)
  data[[i]] = aggregate(new_data$force, by=list(Category=new_data$tf), FUN=sum)$x
  cat(paste(i,"\n"))
  i = i+1
}

outdegree = do.call(cbind, data)

rownames(outdegree) = TFs
colnames(outdegree) = indices

head(outdegree)
dim(outdegree)

setwd("/home/esaha/Aging_thyroid/data/network/")
write.table(outdegree, file="BonoboPanda_GTEx_thyroid_TFoutdegree.txt", row.names=T, col.names=T)
