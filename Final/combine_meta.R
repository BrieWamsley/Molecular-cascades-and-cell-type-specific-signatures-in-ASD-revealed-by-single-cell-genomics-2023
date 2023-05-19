# combine all meta file into one
rm(list=ls())
files<-list.files("./",pattern="^meta",recursive = T)
files

data<-read.csv(files[1],head=T,stringsAsFactors = F)
head(data)
for(i in 2:length(files)){
  data2<-read.csv(files[i],head=T,stringsAsFactors = F)
  data<-rbind(data,data2)
}

write.csv(data,"all_meta.csv")

colnames(data)=c("barcodekey","n_genes1","n_counts1","percent_mito1","NotUse","scale1","louvain_labels1","scrublet_score1","doublet_call")

write.csv(data,"all_meta.csv")

