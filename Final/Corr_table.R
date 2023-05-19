# generate correlation table
rm(list=ls())

a=data.frame(cluster=rep("none",36),cor0.01=rep(0,36),cor0.05=rep(0,36),cor0.1=rep(0,36),corAll=rep(0,36))

for(i in 1:36){
  print(paste0("Cluster_",i))
  # read DEA for DUP15q dataset
  dat1.dup15=read.csv(paste0("./Complete_geneList_",i,".csv"),head=T,stringsAsFactors = F)
  # read DEA for whole dataset
  dat1.all=read.csv(paste0("DEA_results/Complete_geneList_",i,".csv"),head=T,stringsAsFactors = F)

  # merging two datasets
  dat.comb=merge(dat1.dup15,dat1.all,by="gene")
  dim(dat.comb)
  # rank each dataset
  dat.comb$rankDup15=rank(dat.comb$logFC_ASD_vs_CTL.x)
  dat.comb$rankAll=rank(dat.comb$logFC_ASD_vs_CTL.y)

  # DEGs only (FDR<0.1, in all data)
  dat.DEG=dat.comb[dat.comb$pVal_ASD_vs_CTL.y<0.01,];dim(dat.DEG)
  k0.01=cor(dat.DEG$rankDup15, dat.DEG$rankAll,method = c("pearson", "kendall", "spearman"))
  plot(dat.DEG$rankDup15,dat.DEG$rankAll)
  
  dat.DEG=dat.comb[dat.comb$pVal_ASD_vs_CTL.y<0.05,];dim(dat.DEG)
  k0.05=cor(dat.DEG$rankDup15, dat.DEG$rankAll,method = c("pearson", "kendall", "spearman"))
  dat.DEG=dat.comb[dat.comb$pVal_ASD_vs_CTL.y<0.1,];dim(dat.DEG)
  k0.1=cor(dat.DEG$rankDup15, dat.DEG$rankAll,method = c("pearson", "kendall", "spearman"))
  dat.DEG=dat.comb[dat.comb$pVal_ASD_vs_CTL.y<=1,];dim(dat.DEG)
  kall=cor(dat.DEG$rankDup15, dat.DEG$rankAll,method = c("pearson", "kendall", "spearman"))
  print(paste(k0.01,k0.05,k0.1,kall,sep=" "))
  a[i,1]=i;a[i,2]=k0.01;a[i,3]=k0.05;a[i,4]=k0.1;a[i,5]=kall
}

write.csv(a,"Cor_rank_ASDvsCTL_DUP15.csv")
