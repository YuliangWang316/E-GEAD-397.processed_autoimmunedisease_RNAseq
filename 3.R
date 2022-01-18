library(dplyr)
library(tidyverse)
b<-read.table(file = "c:/Users/xjmik/Downloads/E-GEAD-397.processed_autoimmunedisease_RNAseq/counts/USM_B_count.txt",sep="\t",row.names=1,header = TRUE)
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
e<-select(d,starts_with("HC"))
f<-select(d,starts_with("RA"))
h<-cbind(e,f)
write.table(h,file = "c:/Users/xjmik/Desktop/USMB/noGOOD/HC_RA_counts.txt",sep = "\t")
