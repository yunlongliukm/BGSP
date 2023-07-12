 ###(a) Normalization of relative expression levels of the A, B, and C subgenomes
 #bash
 awk '{if(($2+$3+$4) >=0.5)print $1, $2/($2+$3+$4), $3/($2+$3+$4), $4/($2+$3+$4)}' >TPM_logTPM_normolized.txt

##(b) Definition of homoeologous expression bias categories for Normalization of relative expression
#R version 4.1.2
data=read.csv("TPM_logTPM_normolized.txt",header=F, sep="\t",row.names = 1)
dist=as.matrix(dist(data[,1:3],method = "euclidean"))[,1:7]
write.table(t(dist),file="dist", quote=F)

#bash
for i in {9..3518};do awk '{print $"'${i}'", $1}' dist  |sort -n |head -1 ;done > dist_category

##(c) plot the ternary diagrams using the R package ggtern
#R version 4.1.2
library(ggplot2)
library(ggtern)
data=read.table("dist_category.csv",header=T,check.names=F,sep=",",quote="",dec=".")
ggtern(data=data,aes(x=A,y=B,z=C,color=Group ) ) + theme_rgbw() + geom_point(aes(fill=Group),size=1,shape=21)
ggsave("Mhu_combined_Ternary_Plot1.pdf")
