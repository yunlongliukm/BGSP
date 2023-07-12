## Comparative Genome Analysis Codes

The comparative genome analysis involves several steps, each of which is performed by a specific script in the `code` directory. Here is a brief description of each script:

1. `orthologous_gene_clustering.py`: This script performs orthologous gene clustering using the OrthoMCL algorithm.

2. `phylogenetic_tree_construction.py`: This script constructs a phylogenetic tree using the Maximum Likelihood method implemented in the PhyML software.

3. `divergence_time_estimation.py`: This script estimates divergence times using the MCMCTree program in the PAML package.

4. `positive_selection_analysis.py`: This script performs positive selection analysis using the branch-site model in the PAML package.

To run the comparative genome analysis, use the following command:
# To run the comparative genome analysis
./run_comparative_genome_analysis.sh
Please refer to the individual README files in the respective directories for more detailed instructions on how to run the codes. 
