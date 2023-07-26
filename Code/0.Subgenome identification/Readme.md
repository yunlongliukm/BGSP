# Subgenome identification


## 1. Reconstruction of "perfect-copy" syntenic gene data set

   - Input files: `sp.cds`, `sp.bed`, `sp`
   - Script: `Perfect-copy_syntenic_gene_identification.sh`
   - Usage: `bash Perfect-copy_syntenic_gene_identification.sh Alu 2 1`

## 2. Reconstruction of single "perfect-copy" gene matrix and concatenated supermatrix for each syntenic block

   - Input files: `sp.perfect-copy_gene.cds.fa`
   - Script: `gene_matrix_buliding.pl`
   - Usage: `perl gene_matrix_buliding.pl`

## 3. Time tree estimation

   - Input files: `11bg.trees`, `11bg_alignment.txt`, `11bg_mcmctree.ctl`
   - Software: [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html)
   - Usage: `mcmctree 11bg_mcmctree.ctl`

## 4. PhyloNet analyses

   - Input files: `InferNetwork_MPL_430gene_r2.nex`
   - Software: `PhyloNet_3.8.0.jar`
   - Usage: `java -Xmx70g -jar $PHYLONET PhyloNet_3.8.0.jar InferNetwork_MPL_430gene_r2.nex`

## 5. Gene tree topologies calculation

   - Input files: `11bg_430gene.tre`
   - Script: `newick_utils`, `Tree_topologies_calculation.sh`
   - Usage: `bash Tree_topologies_calculation.sh`
