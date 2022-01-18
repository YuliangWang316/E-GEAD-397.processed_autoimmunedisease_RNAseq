setwd("C:/Users/xjmik/Downloads/E-GEAD-397.processed_autoimmunedisease_RNAseq/")
a = list.files("counts/")
dir = paste("./counts/",a,sep="")
n = length(dir)  
library(dplyr)
library(tidyverse)
library(DESeq2)

for ( k in c(1:2,4,7,11:16,18:22,24:length(a))) {
  b<-read.table(file = dir[k],sep="\t",row.names=1,header = TRUE)
  b<-b[match(unique(b[,1]),b$Gene_name), ]
  rownames(b)<-b[,1]
  b<-b[,-1]
  c<-as.data.frame(t(read.table("C:/Users/xjmik/Downloads/E-GEAD-397.processed_autoimmunedisease_RNAseq/clinical_diagnosis_age_sex_v2.txt",sep = "\t",header = TRUE,row.names = 1)))
  c<-c[,match(colnames(b),colnames(c))]
  d<-rbind(c,b)
  for (j in 1:length(colnames(d))) {
    colnames(d)[j] <- paste(d[1,j],colnames(d)[j],sep = "-")
  }
  n<-as.factor(unique(d[1,]))
  d<-d[-1,]
  d<-d[-1,]
  d<-d[-1,]
  d<-d[-1,]
#  for (p in 1:100) {
#    d<-d[d[,p] >=5,]
#  }
  e<-select(d,starts_with("HC"))
  for (o in c(1:3,5:length(levels(n)))) {
    f<-select(d,starts_with(levels(n)[o]))
    h<-cbind(e,f)
    l<-rownames(h)
    h<-apply(h,MARGIN = 2,as.numeric)
    rownames(h)<-l
    condition<-factor(c(rep("HC",length(colnames(e))),rep(levels(n)[o],length(colnames(f)))),levels = c("HC",levels(n)[o]))
    colData<-data.frame(row.names = colnames(h),condition)
    dds <- DESeqDataSetFromMatrix(h, colData, design= ~ condition)
    dds <- DESeq(dds)
    res= results(dds)
    res = res[order(res$pvalue),]
    head(res)
    summary(res)
    write.csv(res,file= paste("new1",a[k],"All_results_",levels(n)[o],".csv",sep = "_"))
  }
  
  
}

for ( k in c(3,5,6,8,9,10,17,23)) {
  b<-read.table(file = dir[k],sep="\t",row.names=1,header = TRUE)
  b<-b[match(unique(b[,1]),b$Gene_name), ]
  rownames(b)<-b[,1]
  b<-b[,-1]
  c<-as.data.frame(t(read.table("C:/Users/xjmik/Downloads/E-GEAD-397.processed_autoimmunedisease_RNAseq/clinical_diagnosis_age_sex_v2.txt",sep = "\t",header = TRUE,row.names = 1)))
  c<-c[,match(colnames(b),colnames(c))]
  d<-rbind(c,b)
  for (j in 1:length(colnames(d))) {
    colnames(d)[j] <- paste(d[1,j],colnames(d)[j],sep = "-")
  }
  n<-as.factor(unique(d[1,]))
  d<-d[-1,]
  d<-d[-1,]
  d<-d[-1,]
  d<-d[-1,]
#  for (p in 1:length(colnames(d))) {
#    d<-d[d[,p] >=5,]
#  }
  e<-select(d,starts_with("HC"))
  for (o in c(2:length(levels(n)))) {
    f<-select(d,starts_with(levels(n)[o]))
    h<-cbind(e,f)
    l<-rownames(h)
    h<-apply(h,MARGIN = 2,as.numeric)
    rownames(h)<-l
    condition<-factor(c(rep("HC",length(colnames(e))),rep(levels(n)[o],length(colnames(f)))),levels = c("HC",levels(n)[o]))
    colData<-data.frame(row.names = colnames(h),condition)
    dds <- DESeqDataSetFromMatrix(h, colData, design= ~ condition)
    dds <- DESeq(dds)
    res= results(dds)
    res = res[order(res$pvalue),]
    head(res)
    summary(res)
    write.csv(res,file= paste("new1",a[k],"All_results_",levels(n)[o],".csv",sep = "_"))
  }
  
  
}
