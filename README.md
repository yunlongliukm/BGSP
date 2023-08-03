# Bamboo Genome Sequencing Project

Bamboo Genome Sequencing Project is focused on the sequencing of 11 bamboo genomes. 
This repository contains workflows and scripts used for the identification and comparative analysis of bamboo subgenomes.

## Reconstruction of "perfect-copy" syntenic gene data set
   - Dependencies: [JCVI v1.1.17](https://github.com/tanghaibao/jcvi)
   - Input files: `sp.cds`, `sp.bed`, `sp`
   - Script: `Perfect-copy_syntenic_gene_identification.sh`
   - Usage: `bash Perfect-copy_syntenic_gene_identification.sh Alu 2 1`
```bash
sp=$1
a=$2
b=$3
python -m jcvi.compara.catalog ortholog --quota=$a:$b $sp Osa --no_strip_names
python -m jcvi.compara.synteny mcscan --iter=1 $sp.bed $sp.Osa.anchors -o $sp.Osa.blocks
cat $sp.Osa.blocks|sort -k2|groupBy -g 2 -c 1,1 -o count,distinct|awk -v num=$a '{ if( $2==num) print $0}'|cut -f 1,3 |sed -e 's/,/\t/g' >$sp.block
```
After running for all species, use the following command to reconstruct "perfect-copy" syntenic gene data set for all species
```bash
python -m jcvi.formats.base join $sp.block >perfect-copy_syntenic_gene.txt
```

## Reconstruction of single "perfect-copy" gene matrix and concatenated supermatrix for each syntenic block
   - Dependencies:
      * [samtools v1.10](https://samtools.sourceforge.net/)
      * [mafft v7.471](https://mafft.cbrc.jp/alignment/software/)
      * [pal2nal v14](http://www.bork.embl.de/pal2nal/)
      * [bioawk](https://github.com/lh3/bioawk)
      * [raxml v8.2.12](https://github.com/stamatak/standard-RAxML)
      * [AMAS](https://github.com/marekborowiec/AMAS)
   - Input files: `cds.fasta`, `pep.fasta`, `gene#.id`
   - Script: `syntenic_gene_matrix.nf`
   - Usage: `nextflow syntenic_gene_matrix.nf`

## Time tree estimation

   - Input files: `11bg.trees`, `11bg_alignment.txt`, `11bg_mcmctree.ctl`
   - Dependencies: [PAML v4.9](http://abacus.gene.ucl.ac.uk/software/paml.html)
   - Usage: `mcmctree 11bg_mcmctree.ctl`

## PhyloNet analyses

   - Input files: `InferNetwork_MPL_430gene_r2.nex`
   - Dependencies: [PhyloNet_3.8.0](https://wiki.rice.edu/confluence/display/PHYLONET/Home)
   - Usage: `java -Xmx70g -jar $PHYLONET PhyloNet_3.8.0.jar InferNetwork_MPL_430gene_r2.nex`

## Gene tree topologies calculation
   - Dependencies: [newick_utils](https://github.com/tjunier/newick_utils)
   - Input files: `11bg_430gene.tre`
   - Script:  `Tree_topologies_calculation.sh`
   - Usage: `bash Tree_topologies_calculation.sh`
```bash Tree_topologies_calculation.sh
nw_reroot 11bg_430gene.tre osa > 11bg_430gene_reroot.tre
nw_condense 11bg_430gene_reroot.tre > 11bg_430gene_reroot_condense.tre
nw_topology 11bg_430gene_reroot_condense.tre > 11bg_430gene_reroot_condense_topology.tre
nw_order 11bg_430gene_reroot_condense_topology.tre |sort|groupBy -g 1 -c 1 -o count
```

## Subgenome dominance of hexaploid

- Normalization of relative expression levels of the A, B, and C subgenomes

```bash
awk '{if(($2+$3+$4) >=0.5)print $1, $2/($2+$3+$4), $3/($2+$3+$4), $4/($2+$3+$4)}' ABC_111_TPM.csv >TPM_logTPM_normolized.txt
```

- Definition of homoeologous expression bias categories for Normalization of relative expression (with M. baccifera as example).
```R
#R version 4.1.2
data=read.csv("TPM_logTPM_normolized.txt",header=F, sep="\t",row.names = 1)
dist=as.matrix(dist(data[,1:3],method = "euclidean"))[,1:7]
write.table(t(dist),file="dist", quote=F)
```
```bash
for i in {9..3518};do awk '{print $"'${i}'", $1}' dist  |sort -n |head -1 ;done > dist_category.csv
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
