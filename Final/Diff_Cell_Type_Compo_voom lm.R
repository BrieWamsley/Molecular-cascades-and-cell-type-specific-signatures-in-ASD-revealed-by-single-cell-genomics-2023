rm(list=ls())

# Load necessary R-packages
library(GenomicFeatures)
library(heatmap.plus)
library(MASS)
library(DESeq)
library(edgeR)
library(biomaRt)
library(Biobase)
library(limma)
library(RColorBrewer)
library(gcrma)  
library(gplots)
library(annaffy)
library(plyr)
library(ggplot2)

# Load cell number matrix (cluster x sample)
rawCount<-read.csv("cell_composition_dim65res2.0.csv",head=T,row.names=1,stringsAsFactors = F);head(rawCount)
rownames(rawCount)<-paste("Cluster",rownames(rawCount))                   

# Load meta data
conds<-read.csv("target.csv",head=T,stringsAsFactors=F,row.name=1)
conds<-conds[colnames(rawCount),]

# assign a color for each trait in the meta data
cv <- as.matrix(conds)
dim(cv) <- c(1,prod(dim(cv)))
rownames(table(cv))
names <- rownames(table(cv))
col_default <- c("red ","blue","green","orange","bisque4 ","black ","brown ","cyan ","darkgreen ","darkgrey ","darkmagenta ","darkolivegreen ","darkorange ","darkred ","darkslateblue ","darkturquoise ","floralwhite ","greenyellow ","grey ","lightcyan ","lightcyan1 ","lightgreen ","lightsteelblue1 ","lightyellow ","magenta ","mediumpurple3 ","midnightblue ","paleturquoise ","pink ","plum1 ","plum2 ","royalblue ","saddlebrown ","salmon ","sienna3 ","skyblue ","skyblue3 ","steelblue ","tan ","thistle1 ","thistle2 ","turquoise ","violet ","white ","yellowgreen","grey60 ","orangered4 ","brown4 ","darkorange2 ","ivory ")
col <- rainbow(length(names))
names(col)<- names
clab <- matrix(col[as.matrix(conds)],nrow=nrow(conds),ncol=ncol(conds))
colnames(clab) <- colnames(conds) 
rownames(clab)<-rownames(conds)
for(j in 1:ncol(conds)){
  a=1
  for (i in unique(clab[,j])){
    clab[which(clab[,j]==i),j]=col_default[a]
    a=a+1
  }
}
clab

# Linear model function to test various models 
LM<-function(norm,design,model){
  # we apply limma voom to run linear model using TMM normalization method
  nf = calcNormFactors(rawCount, method = "TMM")
  group <- conds$Diagnosis
  group<-factor(group,levels=unique(group)) #group as a factor
  names(group) <- rownames(conds)
  voom.data = voom(rawCount, design =design, lib.size = colSums(rawCount) * nf)
  voom.data$Cluster = rownames(rawCount)
  voom.fitlimma = lmFit(voom.data, design = design)
  voom.fitbayes = eBayes(voom.fitlimma)
  voom.pvalues = voom.fitbayes$p.value[, 2]
  voom.adjpvalues = p.adjust(voom.pvalues, method = "BH")
  
  ##################################################################
  my.contrasts<-makeContrasts("ASD_vs_CTL"=ASD - CTL,
                              "CTL_vs_ASD"=CTL - ASD,
                              levels=design)
  ####################################################################
  
  fit2.anova<- contrasts.fit(voom.fitlimma, my.contrasts) 
  fitb<- eBayes(fit2.anova) 

  Voom.logFC<-fit2.anova$coefficients
  colnames(Voom.logFC)<-paste("logFC_",colnames(Voom.logFC),sep="")
  Voom.pvalues<-fitb$p.value
  colnames(Voom.pvalues)<-paste("pVal_",colnames(Voom.pvalues),sep="")
  Voom.adjpvalues<-Voom.pvalues
  for (i in 1:dim(my.contrasts)[2]){
    Voom.adjpvalues[,i]<-p.adjust(Voom.pvalues[,i],method="BH")
  }
  colnames(Voom.adjpvalues)<-paste("FDR_",colnames(Voom.adjpvalues),sep="")
  Rank <- rank(rowSums(-voom.data$E), ties.method="first")
  voom.data$E.count<-2^voom.data$E
  colnames(voom.data$E.count)<-paste(colnames(voom.data$E.count),"_voom_count",sep="")
  Voom.count<-rawCount
  colnames(Voom.count)<-paste(colnames(Voom.count),"_raw_count",sep="")
  Complete<-as.data.frame( cbind(Voom.logFC,Voom.pvalues,Voom.adjpvalues))
  head(Complete)
  tail(Complete)
  cpComp1<-Complete
  write.csv(Complete, file= paste0("Cell_composition_analysis_",norm,"_",model,".csv"))
}

group <- conds$Diagnosis
group<-factor(group,levels=unique(group)) #group as a factor
names(group) <- rownames(conds)
design<-model.matrix(~0+group)
colnames(design)<-gsub("group","",colnames(design))

model="Diagnosis"
LM("TMM",design,model)

BrainBank=factor(conds$BrainBank,levels=unique(conds$BrainBank));names(BrainBank)=rownames(conds)

design<-model.matrix(~0+group+BrainBank)
colnames(design)<-gsub("group|UM-BTB","",colnames(design))

model="Diagnosis+BrainBank"
LM("TMM",design,model)

rnaBatch=factor(conds$rnaBatch,levels=unique(conds$rnaBatch));names(rnaBatch)=rownames(conds)

design<-model.matrix(~0+group+rnaBatch)
colnames(design)<-gsub("group|rna","",colnames(design))

model="Diagnosis+rnaBatch"
LM("TMM",design,model)

###################

design<-model.matrix(~0+group+BrainBank+rnaBatch)
colnames(design)<-gsub("group|UM-BTB|rna","",colnames(design))

model="Diagnosis+BrainBank+rnaBatch"
LM("TMM",design,model)

###################
Age=conds$Age

design<-model.matrix(~0+group+BrainBank+rnaBatch+Age);design
colnames(design)<-gsub("group|UM-BTB|rna","",colnames(design))

model="Diagnosis+BrainBank+rnaBatch+Age"
LM("TMM",design,model)
###################
PMI=factor(paste0("X",conds$PMI),levels=unique(paste0("X",conds$PMI)));names(PMI)=rownames(conds)

PMI
design<-model.matrix(~0+group+BrainBank+rnaBatch+Age+PMI);design
colnames(design)<-gsub("group|UM-BTB|rna|PMI","",colnames(design))

model="Diagnosis+BrainBank+rnaBatch+Age+PMI"
LM("TMM",design,model)

################## adding covariate one by one #######################

design<-model.matrix(~0+group+rnaBatch);design
colnames(design)<-gsub("group|rna","",colnames(design))

model="Diagnosis+rnaBatch"
LM("TMM",design,model)
################## 
design<-model.matrix(~0+group+Age);design
colnames(design)<-gsub("group","",colnames(design))

model="Diagnosis+Age"
LM("TMM",design,model)
################## 
design<-model.matrix(~0+group+PMI);design
colnames(design)<-gsub("group|PMI","",colnames(design))

model="Diagnosis+PMI"
LM("TMM",design,model)
################## 
Chemistry10X=factor(conds$Chemistry10X,levels=unique(conds$Chemistry10X));names(PMI)=rownames(conds)

design<-model.matrix(~0+group+Chemistry10X);design
colnames(design)<-gsub("group|Chemistry10X","",colnames(design))

model="Diagnosis+Chemistry10X"
LM("TMM",design,model)

################## 

design<-model.matrix(~0+group+Chemistry10X+Age);design
colnames(design)<-gsub("group|Chemistry10X","",colnames(design))

model="Diagnosis+Chemistry10X+Age"
LM("TMM",design,model)
################## 
design<-model.matrix(~0+group+Chemistry10X+Age+rnaBatch);design
colnames(design)<-gsub("group|Chemistry10X|rna","",colnames(design))

model="Diagnosis+Chemistry10X+Age+rnaBatch"
LM("TMM",design,model)

colnames(conds)

################## 
nGenes=conds$nGenes
design<-model.matrix(~0+group+Chemistry10X+Age+nGenes);design
colnames(design)<-gsub("group","",colnames(design))

model="Diagnosis+Chemistry10X+Age+nGenes"
LM("TMM",design,model)

################## 
conds$Race[is.na(conds$Race)]="notKnown"
Race=factor(conds$Race)
design<-model.matrix(~0+group+Chemistry10X+Age+Race);design
colnames(design)<-gsub("group","",colnames(design))

model="Diagnosis+Chemistry10X+Age+Race"
LM("TMM",design,model)

################## 
conds$Race[is.na(conds$Race)]="notKnown"
Race=factor(conds$Race)
design<-model.matrix(~0+group+Chemistry10X+Age+Race);design
colnames(design)<-gsub("group","",colnames(design))

model="Diagnosis+Chemistry10X+Age+Race"
LM("TMM",design,model)
