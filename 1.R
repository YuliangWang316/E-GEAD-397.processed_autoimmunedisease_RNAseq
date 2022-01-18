setwd("C:/Users/xjmik/Downloads/E-GEAD-397.processed/")
a = list.files("counts/")
dir = paste("./counts/",a,sep="")
n = length(dir)  
library(dplyr)
library(tidyverse)
library(DESeq2)

for ( k in c(1:12,14:length(a))) {
  b<-read.table(file = dir[k],sep="\t",row.names=1,header = TRUE)
  b<-b[match(unique(b[,1]),b$Gene_name), ]
  rownames(b)<-b[,1]
  b<-b[,-1]
  c<-as.data.frame(t(read.table("C:/Users/xjmik/Downloads/E-GEAD-397.processed/clinical_diagnosis_age_sex_v2.txt",sep = "\t",header = TRUE,row.names = 1)))
  c<-c[,match(colnames(b),colnames(c))]
  d<-rbind(c,b)
  for (j in 1:length(colnames(d))) {
    colnames(d)[j] <- paste(d[1,j],colnames(d)[j],sep = "-")
  }
  d<-d[-1,]
  d<-d[-1,]
  d<-d[-1,]
  d<-d[-1,]
  e<-select(d,starts_with("HC"))
  f<-select(d,starts_with("SLE"))
  g<-select(d,starts_with("RA"))
  h<-cbind(e,f)
  i<-cbind(e,g)
  l<-rownames(h)
  m<-rownames(i)
  h<-apply(h,MARGIN = 2,as.numeric)
  i<-apply(i,MARGIN = 2,as.numeric)
  rownames(h)<-l
  rownames(i)<-m
  condition<-factor(c(rep("HC",length(colnames(e))),rep("SLE",length(colnames(f)))),levels = c("HC","SLE"))
  colData<-data.frame(row.names = colnames(h),condition)
  
  dds <- DESeqDataSetFromMatrix(h, colData, design= ~ condition)
  dds <- DESeq(dds)
  
  res= results(dds)
  res = res[order(res$pvalue),]
  head(res)
  summary(res)
  write.csv(res,file= paste(a[k],"All_results_SLE.csv",sep = "_"))
  condition<-factor(c(rep("HC",length(colnames(e))),rep("RA",length(colnames(g)))),levels = c("HC","RA"))
  colData<-data.frame(row.names = colnames(i),condition)
  
  dds <- DESeqDataSetFromMatrix(i, colData, design= ~ condition)
  dds <- DESeq(dds)
  
  res= results(dds)
  res = res[order(res$pvalue),]
  head(res)
  summary(res)
  write.csv(res,file= paste(a[k],"All_results_RA.csv"))
  
}

b<-read.table(file = dir[13],sep="\t",row.names=1,header = TRUE)
b<-b[match(unique(b[,1]),b$Gene_name), ]
rownames(b)<-b[,1]
b<-b[,-1]
c<-as.data.frame(t(read.table("C:/Users/xjmik/Downloads/E-GEAD-397.processed/clinical_diagnosis_age_sex_v2.txt",sep = "\t",header = TRUE,row.names = 1)))
c<-c[,match(colnames(b),colnames(c))]
d<-rbind(c,b)
for (j in 1:length(colnames(d))) {
  colnames(d)[j] <- paste(d[1,j],colnames(d)[j],sep = "-")
}
d<-d[-1,]
d<-d[-1,]
d<-d[-1,]
d<-d[-1,]
e<-select(d,starts_with("HC"))
f<-select(d,starts_with("SLE"))
g<-select(d,starts_with("RA"))
h<-cbind(e,f)
i<-cbind(e,g)
l<-rownames(h)
m<-rownames(i)
h<-apply(h,MARGIN = 2,as.numeric)
i<-apply(i,MARGIN = 2,as.numeric)
rownames(h)<-l
rownames(i)<-m
#condition<-factor(c(rep("HC",length(colnames(e))),rep("SLE",length(colnames(f)))),levels = c("SLE","HC"))
#colData<-data.frame(row.names = colnames(h),condition)

#dds <- DESeqDataSetFromMatrix(h, colData, design= ~ condition)
#dds <- DESeq(dds)

#res= results(dds)
#res = res[order(res$pvalue),]
#head(res)
#summary(res)
#write.csv(res,file= paste(a[k],"All_results_SLE.csv",sep = "_"))
condition<-factor(c(rep("HC",length(colnames(e))),rep("RA",length(colnames(g)))),levels = c("HC","RA"))
colData<-data.frame(row.names = colnames(i),condition)

dds <- DESeqDataSetFromMatrix(i, colData, design= ~ condition)
dds <- DESeq(dds)

res= results(dds)
res = res[order(res$pvalue),]
head(res)
summary(res)
write.csv(res,file= paste(a[13],"All_results_RA.csv"))
