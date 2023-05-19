# pseudo-bulk analysis of ASD snRNAseq
rm(list=ls())
library(Seurat)
library(GenomicFeatures)
library(heatmap.plus)
library(MASS)
library(DESeq2)
library(edgeR)
library(biomaRt)
library(Biobase)
library(limma)
library(RColorBrewer)
library(gcrma)  
library(gplots)
library(ggplot2)
library(plyr)
library(WGCNA)

f<-list.files("../Seurat_Obj/")
f

for(k in c(1)){
  s1<-readRDS(paste0("../Seurat_Obj/",f[k]))
  samples<-unique(s1@meta.data$Channel)
  
  # generate count table with gene x sample
  dim(s1@assays$RNA)
  rawC<-matrix(nrow=dim(s1@assays$RNA)[1],ncol=length(unique(s1@meta.data$Channel)));dim(rawC)
  rawC<-as.data.frame(rawC)
  
  for(i in 1:length(samples)){
    print(samples[i])
    dat.sub<-s1@assays$RNA[,grep(samples[i],colnames(s1@assays$RNA))]
    print(dim(dat.sub))
    rawC[,i]=apply(dat.sub,1,sum);colnames(rawC)[i]=samples[i]
  }
  rownames(rawC)<-rownames(s1@assays$RNA)
  
  # Loading meta data
  conds<-read.csv("Sample_Meta_ASDscRNAseq.csv",head=T,row.names=1,stringsAsFactors = F)
  conds$sample<-rownames(conds)
  table(colnames(rawC)%in%conds$sample)
  conds<-conds[colnames(rawC),]
  conds$nGenes=scale(conds$nGene)
  conds$PMI=scale(conds$PMI)
  conds$nCells=scale(conds$nCells)
  conds$totalReads=scale(conds$totalReads)
  conds$RIN=scale(conds$RIN)
  conds$medPercMito=0
  for(i in conds$sample){
    print(i)
    conds[i,"medPercMito"]=median(s1@meta.data$percent_mito[grep(i,rownames(s1@meta.data))])
  }
  conds$medPercMito=scale(conds$medPercMito)
  
  # Assign colors for each covariate
  col_default <- c("green","blue","red","brown4 ","black ","brown ","cyan ","darkgreen ","darkgrey ","darkmagenta ","darkolivegreen ","darkorange ","darkred ","darkslateblue ","darkturquoise ","floralwhite ","green ","greenyellow ","grey ","lightcyan ","lightcyan1 ","lightgreen ","lightsteelblue1 ","lightyellow ","magenta ","mediumpurple3 ","midnightblue ","orange ","paleturquoise ","pink ","plum1 ","plum2 ","royalblue ","saddlebrown ","salmon ","sienna3 ","skyblue ","skyblue3 ","steelblue ","tan ","thistle1 ","thistle2 ","turquoise ","violet ","white ","yellowgreen","grey60 ","orangered4 ","bisque4 ","darkorange2 ","ivory ")
  
  cv <- as.matrix(conds)
  dim(cv) <- c(1,prod(dim(cv)))
  rownames(table(cv))
  names <- rownames(table(cv))
  col <- rainbow(length(names))
  names(col)<- names
  clab <- matrix(col[as.matrix(conds)],nrow=nrow(conds),ncol=ncol(conds))
  colnames(clab) <- colnames(conds) 
  rownames(clab)<-rownames(conds)
  for(j in 2:ncol(conds)){
    a=1
    for (i in unique(clab[,j])){
      clab[which(clab[,j]==i),j]=col_default[a]
      a=a+1
    }
  }

  for(i in 1:ncol(rawC)){rawC[,i]<-as.integer(rawC[,i])}
  
  # Filter low expressed gene
  temp<-rawC;temp$flag<-0
  for (i in 1:ncol(rawC)){
    temp$flag<-temp$flag+(temp[,i]>=5)
  }
  table(temp$flag>=5) #22208 clusters
  rawC<-rawC[temp$flag>=5,]
  
  conds$group<-conds$Diagnosis
  cds <- DESeqDataSetFromMatrix(countData=rawC,colData=conds,design=~Diagnosis)
  cds <- estimateSizeFactors( cds )
  cds <- estimateDispersions( cds  )
  vsd <- getVarianceStabilizedData( cds )
  
  # Generate correlation table pre-adjustment
  pdf(paste0("PCA_corr_",gsub(".rds","",f[k]),".pdf"),height=20,width=24)
  thisdat.HTSC <- t(scale(t(vsd),scale=F))
  PC.HTSC <- prcomp(thisdat.HTSC,center=F);
  TopPC1 <- PC.HTSC$rotation[,1:5];
  varexp <- (PC.HTSC$sdev)^2 / sum(PC.HTSC$sdev^2)
  topvar <- varexp[1:5]
  colnames(TopPC1) <- paste("Exp.PC\n",colnames(TopPC1)," (",signif(100*topvar[1:5],2),"%)",sep="")
  
  pairsdat <- data.frame(group=as.numeric(factor(conds$group)),
                         BrainBrank=as.numeric(factor(conds$BrainBrank)),
                         Race=as.numeric(factor(conds$Race)),
                         SexChromosome=as.numeric(factor(conds$SexChromosome)),
                         Clinical=as.numeric(factor(conds$Clinical)),
                         rnaBatch=as.numeric(factor(conds$rnaBatch)),
                         Chemistry10X=as.numeric(factor(conds$Chemistry10X)),
                         BrainRegion=as.numeric(factor(conds$BrainRegion)),
                         PMI=as.numeric(conds$PMI),
                         RIN=as.numeric(conds$RIN),
                         nCells=as.numeric(conds$nCells),
                         totalReads=as.numeric(conds$ totalReads),
                         nGenes=as.numeric(conds$nGenes),
                         medPercMito=as.numeric(conds$medPercMito)
              
    )
  rownames(pairsdat)<-conds$sample
  
  cond=labels2colors(as.numeric(factor(conds$sample)))  ## colors
  
  panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) { ## Useful function for comparing multivariate data
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y,use="pairwise.complete.obs",method="pearson"))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
  }
  
  
  panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) { ## Useful function for comparing multivariate data
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    if (class(x) == "numeric" & class(y) == "numeric") {
      r <- abs(cor(x, y,use="pairwise.complete.obs",method="pearson"))
    } else {
      lmout <- lm(y~x)
      r <- sqrt(summary(lmout)$adj.r.squared)
    }
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
  }
  
  pairs(cbind(TopPC1,pairsdat),col= cond,pch=19,upper.panel = panel.cor,main="Covariates and Normalized CellType Composition Comparison -- |Spearman's rho| correlation values")
  
  dev.off()
  
  
  #### Calculate median beta by bootstrapping
  
  group=as.numeric(factor(conds$group))
  Clinical=as.numeric(factor(conds$Clinical))
  rnaBatch=as.numeric(factor(conds$rnaBatch))
  Chemistry10X=as.numeric(factor(conds$Chemistry10X))
  BrainRegion=as.numeric(factor(conds$BrainRegion))
  PMI=as.numeric(conds$PMI)
  RIN=as.numeric(conds$RIN)
  nCells=as.numeric(conds$nCells)
  medPercMito=as.numeric(conds$medPercMito)
  
  regvars <- as.data.frame(cbind(group,Chemistry10X,rnaBatch,BrainRegion,nCells,medPercMito))
  varnames <- c("group","Chemistry10X","Clinical","rnaBatch","BrainRegion","nCells","medPercMito")
  
  
  X <-  model.matrix(~group+Chemistry10X+rnaBatch+BrainRegion+nCells+medPercMito,data = regvars)
  rownames(X) <- colnames(vsd)
  Y <- vsd
  
  # perform 500 samplings with replacement of the samples
  beta_samplings <- list()
  set.seed(8675309)
  original_sample=colnames(vsd)
  for (j in 1:500) {
    print(j)
    keepSamp <- sample(colnames(Y),length(colnames(Y)),replace=T)
    print(table(original_sample%in%keepSamp))
    X_sample <- X[keepSamp,]
    Y_sample <- Y[,keepSamp]
    beta <- (solve(t(X_sample)%*%X_sample)%*%t(X_sample))%*%t(Y_sample)
    beta_samplings[[j]] <- beta
  }

  # identify median beta value for each covariate for each gene across the bootstraps
  a <- array(unlist(beta_samplings), c(nrow(beta), ncol(beta), length(beta_samplings)))
  beta_median <- apply(a, 1:2, median)
  
  # regress out median beta values from bootstraps
  to_regress <- (as.matrix(X[,3:ncol(X)]) %*% (as.matrix(beta_median[3:ncol(X),])))
  normExpr.reg <- vsd - t(to_regress)
  
  ## This datExpr.reg is now a technical variable corrected matrix.
  rownames(normExpr.reg) <- rownames(vsd)
  colnames(normExpr.reg) <- colnames(vsd)
  
  dim(normExpr.reg)
  write.csv(normExpr.reg,paste0("adj_normExpr_exp",gsub(".rds","",f[k]),".csv"))
  
  pdf(paste0("PCA_corr_postAdj_",gsub(".rds","",f[k]),".pdf"),height=20,width=24)
  thisdat.HTSC <- t(scale(t(normExpr.reg),scale=F))
  PC.HTSC <- prcomp(thisdat.HTSC,center=F);
  TopPC2 <- PC.HTSC$rotation[,1:5];
  varexp <- (PC.HTSC$sdev)^2 / sum(PC.HTSC$sdev^2)
  topvar <- varexp[1:5]
  colnames(TopPC2) <- paste("Exp.PC\n",colnames(TopPC2)," (",signif(100*topvar[1:5],2),"%)",sep="")
  
  pairs(cbind(TopPC2,pairsdat[rownames(TopPC2),]),col= cond,pch=19,upper.panel = panel.cor,main="Covariates and HTSeqCount Comparison -- |Spearman's rho| correlation values")
  
  dev.off()
  
  ### DEA with limma VSD 

  adj.log.norm<-normExpr.reg[,conds$sample]
  DESeq.vst = adj.log.norm
  group <- conds$group
  group<-factor(group,levels=unique(group)) #group as a factor
  names(group) <- rownames(conds)
  design<-model.matrix(~0+group)
  colnames(design)<-gsub("group","",colnames(design))
  DESeq.vst.fitlimma = lmFit(DESeq.vst, design = design)
  ####################################################
  sizeText=0.7
  pdf(paste0("MDSplot_",gsub(".rds","",f[k]),".pdf"), width=10,height=10)
  data <- adj.log.norm
  ldat <- dist(t(data))
  fit <- isoMDS(ldat, k=2)
  x <- fit$points[,1]
  y <- fit$points[,2]
  plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="PCA corrected normalized counts", type="n",cex=1,col=clab[,grep("Diagnosis",colnames(clab))],xlim=c(min(x)*1.4,max(x)*1.4))
  text(x, y, labels = row.names(t(data)), cex=sizeText, col=clab[,grep("Diagnosis",colnames(clab))])
  legend("topright", legend=unique(conds$group), fill=unique(clab[,grep("Diagnosis",colnames(clab))]),cex=sizeText)
  dev.off()
  ####################################################
  my.contrasts<-makeContrasts("ASD_vs_CTL"=ASD - CTL,
                              "CTL_vs_ASD"=CTL - ASD,
                              levels=design)
  ####################################################
  
  fit2.anova<- contrasts.fit(DESeq.vst.fitlimma, my.contrasts) 
  fitb<- eBayes(fit2.anova) 
  
  #selecting the statistical cutoff
  padj=0.1
  ########################################
  # Get FC & p-values & save the output  #
  ########################################
  DESeq.vst.FC<-fit2.anova$coefficients
  colnames(DESeq.vst.FC)<-paste("logFC_",colnames(DESeq.vst.FC),sep="")
  DESeq.vst.pvalues<-fitb$p.value
  colnames(DESeq.vst.pvalues)<-paste("pVal_",colnames(DESeq.vst.pvalues),sep="")
  DESeq.vst.adjpvalues<-DESeq.vst.pvalues
  for (i in 1:dim(my.contrasts)[2]){
    DESeq.vst.adjpvalues[,i]<-p.adjust(DESeq.vst.pvalues[,i],method="BH")
  }
  colnames(DESeq.vst.adjpvalues)<-paste("FDR_",colnames(DESeq.vst.adjpvalues),sep="")
  Rank <- rank(rowSums(-DESeq.vst), ties.method="first")
  DESeq.vst.count<-DESeq.vst
  colnames(DESeq.vst.count)<-paste(colnames(DESeq.vst.count),"_vst_count",sep="")
  DESeq.count<-2^adj.log.norm
  colnames(DESeq.count)<-paste(colnames(DESeq.count),"adj_count",sep="")
  Complete<-as.data.frame(cbind(as.data.frame(rownames(adj.log.norm)),Rank,DESeq.vst.FC,DESeq.vst.pvalues,DESeq.count))
  colnames(Complete)[1]="gene"
  dim(Complete)
  cpComp2<-Complete
  write.csv(Complete, file= paste0("Complete_geneList_",gsub(".rds","",f[k]),".csv"))
  
  pdf(paste0("ContrastAnalysis_",gsub(".rds","",f[k]),".pdf"))
  lfc <- DESeq.vst.FC
  lfc <- cbind(lfc, DESeq.vst.adjpvalues)
  ups <- NA
  downs <- NA
  numComp<-dim(DESeq.vst.FC)[2]
  for(i in 1:numComp)
  {
    downs <- c(downs, -table(lfc[,i]<0&lfc[,i+numComp]<padj)[2])
    ups <- c(ups, table(lfc[,i]>0&lfc[,i+numComp]<padj)[2])   
  }
  ups[is.na(ups)]=0
  downs[is.na(downs)]=0
  ups[1]=NA
  downs[1]=NA
  print(ups)
  print(downs)
  names(ups)=c("","")
  names(downs)=c("","")
  mx <- max(ups[is.na(ups)==FALSE])
  mn <- min(downs[is.na(downs)==FALSE])
  bp1 <- barplot(downs,horiz=TRUE,xlim=c(mn,mx), col="green", ) 
  bp2 <- barplot(ups, horiz=TRUE,xlim=c(mn,mx), col="red",add=TRUE,axes=TRUE)
  axis(2, at=bp1[2:length(bp1)],tick=FALSE,labels=downs[is.na(downs)==FALSE],las=1)
  axis(4,at=bp2[2:length(bp2)],tick=FALSE,labels=ups[is.na(ups)==FALSE], las=1,line=-2)
  labs <- colnames(my.contrasts)
  labs <- labs[is.na(labs)==FALSE]
  text(x=0,y=bp2[2:length(bp2)], labels=gsub("logFC_","",labs))
  title(main=paste("Differentially Expressed Genes @ FDR<",padj,sep=""))
  dev.off()
  
  ######## generate rnk files ###################################################
  tempRNK<-cpComp2[,c(1,grep("^pVal_",colnames(cpComp2)),grep("logFC_",colnames(cpComp2)))]
  library(plyr)
  avg<-function(x){
    return(mean(x[,2]))
  }
  
  tempRNKsub<-tempRNK[,c(1,2,4)]
  # convert p-value to directional inversed p-value
  tempRNKsub[,2]<-log10(tempRNK[,2])*((tempRNK[,3]>0)*(-1)+(tempRNK[,3]<0))
  
  
  tempRNKsub<-tempRNKsub[tempRNKsub$gene!="",]
    # remove duplicate by applying average
  tempRNKfin<-ddply(tempRNKsub[,c(1,2)],.(gene),avg)
  tempRNKfin<-tempRNKfin[order(-tempRNKfin$V1),]
  dim(tempRNKfin)
  tempRNKfin$gene<-toupper(tempRNKfin$gene)
  fname<-gsub("pVal_","Signed_log10pVal_",colnames(tempRNK)[i],)
  tempRNKfin<-tempRNKfin[complete.cases(tempRNKfin),]
  write.table(tempRNKfin,file=paste0(fname,gsub(".rds","",f[k]),".rnk"),col.names=FALSE,sep="\t",row.name=FALSE,quote=FALSE)
}

# Make a list of DEG numbers
f2<-list.files("./",pattern="Complete")
f2

pVal=0.001
m1<-matrix(nrow=2,ncol=length(f2))
for(i in 1:length(f2)){
  dat.temp<-read.csv(paste0("./Counts5sample5/",f2[i]),head=T,row.names=1,stringsAsFactors = F)
  upGenes<-dat.temp$logFC_ASD_vs_CTL>0&dat.temp$pVal_ASD_vs_CTL<pVal
  downGenes<-dat.temp$logFC_ASD_vs_CTL<0&dat.temp$pVal_ASD_vs_CTL<pVal
  m1[1,i]=length(upGenes[upGenes=="TRUE"])
  m1[2,i]=length(downGenes[downGenes=="TRUE"])
}
m1<-as.data.frame(m1)
colnames(m1)=gsub(".csv","",gsub("Complete_geneList_","",f2))
rownames(m1)=c(paste0("up_pval<",pVal),paste0("down_pval<",pVal))
pdf(paste0("DEG-barplot_",pVal,".pdf"),height=15,width=6)
barplot(as.matrix(m1),cex.names=0.6,horiz=TRUE,legend = rownames(m1), beside=TRUE)
dev.off()
