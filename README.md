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


