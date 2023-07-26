nw_reroot 11bg_430gene.tre osa > 11bg_430gene_reroot.tre
nw_condense 11bg_430gene_reroot.tre > 11bg_430gene_reroot_condense.tre
nw_topology 11bg_430gene_reroot_condense.tre > 11bg_430gene_reroot_condense_topology.tre
nw_order 11bg_430gene_reroot_condense_topology.tre |sort|groupBy -g 1 -c 1 -o count