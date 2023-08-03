# Bamboo Genome Sequencing Project

Bamboo Genome Sequencing Project is focused on the sequencing of 11 bamboo genomes. 
This repository contains workflows and scripts used for the identification and comparative analysis of bamboo subgenomes.

## Reconstruction of "perfect-copy" syntenic gene data set

   - Input files: `sp.cds`, `sp.bed`, `sp`
   - Script: `Perfect-copy_syntenic_gene_identification.sh`
   - Usage: `bash Perfect-copy_syntenic_gene_identification.sh Alu 2 1`

## Reconstruction of single "perfect-copy" gene matrix and concatenated supermatrix for each syntenic block

   - Input files: `cds.fasta`, `pep.fasta`, `gene#.id`
   - Script: `2.syntenic_gene_matrix.nf`
   - Usage: `nextflow 2.syntenic_gene_matrix.nf`

## Time tree estimation

   - Input files: `11bg.trees`, `11bg_alignment.txt`, `11bg_mcmctree.ctl`
   - Software: [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html)
   - Usage: `mcmctree 11bg_mcmctree.ctl`

## PhyloNet analyses

   - Input files: `InferNetwork_MPL_430gene_r2.nex`
   - Software: `PhyloNet_3.8.0.jar`
   - Usage: `java -Xmx70g -jar $PHYLONET PhyloNet_3.8.0.jar InferNetwork_MPL_430gene_r2.nex`

## Gene tree topologies calculation

   - Input files: `11bg_430gene.tre`
   - Script: `newick_utils`, `Tree_topologies_calculation.sh`
   - Usage: `bash Tree_topologies_calculation.sh`

## Subgenome dominance of hexaploid

- Normalization of relative expression levels of the A, B, and C subgenomes

```bash
awk '{if(($2+$3+$4) >=0.5)print $1, $2/($2+$3+$4), $3/($2+$3+$4), $4/($2+$3+$4)}' >TPM_logTPM_normolized.txt
```

- Definition of homoeologous expression bias categories for Normalization of relative expression (with M. baccifera as example).
```R
#R version 4.1.2
data=read.csv("TPM_logTPM_normolized.txt",header=F, sep="\t",row.names = 1)
dist=as.matrix(dist(data[,1:3],method = "euclidean"))[,1:7]
write.table(t(dist),file="dist", quote=F)
```
```bash
for i in {9..3518};do awk '{print $"'${i}'", $1}' dist  |sort -n |head -1 ;done > dist_category
```

- Plot the ternary diagrams using the R package ggtern
```R
#R version 4.1.2
library(ggplot2)
library(ggtern)
data=read.table("dist_category.csv",header=T,check.names=F,sep=",",quote="",dec=".")
ggtern(data=data,aes(x=A,y=B,z=C,color=Group ) ) + theme_rgbw() + geom_point(aes(fill=Group),size=1,shape=21)
ggsave("Mhu_combined_Ternary_Plot1.pdf")
```
