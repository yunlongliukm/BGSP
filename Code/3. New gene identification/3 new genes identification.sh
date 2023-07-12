###new genes identification

orthofinder -f newgene -a 90
awk '{if($2>0)print $1}' Orthogroups.GeneCount.tsv > Alu_C
awk '{if($3>0)print $1}' Orthogroups.GeneCount.tsv > Alu_D
awk '{if($4>0)print $1}' Orthogroups.GeneCount.tsv > Bam_A
awk '{if($5>0)print $1}' Orthogroups.GeneCount.tsv > Bam_B
awk '{if($6>0)print $1}' Orthogroups.GeneCount.tsv > Bam_C
awk '{if($7>0)print $1}' Orthogroups.GeneCount.tsv > Dsi_A
awk '{if($8>0)print $1}' Orthogroups.GeneCount.tsv > Dsi_B
awk '{if($9>0)print $1}' Orthogroups.GeneCount.tsv > Dsi_C
awk '{if($10>0)print $1}' Orthogroups.GeneCount.tsv > Gan_B
awk '{if($11>0)print $1}' Orthogroups.GeneCount.tsv > Gan_C
awk '{if($12>0)print $1}' Orthogroups.GeneCount.tsv > Hca_C
awk '{if($13>0)print $1}' Orthogroups.GeneCount.tsv > Hca_D
awk '{if($14>0)print $1}' Orthogroups.GeneCount.tsv > Mhu_A
awk '{if($15>0)print $1}' Orthogroups.GeneCount.tsv > Mhu_B
awk '{if($16>0)print $1}' Orthogroups.GeneCount.tsv > Mhu_C
awk '{if($17>0)print $1}' Orthogroups.GeneCount.tsv > Ogl_B
awk '{if($18>0)print $1}' Orthogroups.GeneCount.tsv > Ogl_C
awk '{if($19>0)print $1}' Orthogroups.GeneCount.tsv > P10_Ola
awk '{if($20>0)print $1}' Orthogroups.GeneCount.tsv > P10_Rgu
awk '{if($21+$22+$23+$24+$25>0)print $1}' Orthogroups.GeneCount.tsv > P1
awk '{if($26+$27+$28+$29+$30>0)print $1}' Orthogroups.GeneCount.tsv > P2
awk '{if($31+$32+$33>0)print $1}' Orthogroups.GeneCount.tsv > P3
awk '{if($34>0)print $1}' Orthogroups.GeneCount.tsv > P4
awk '{if($35>0)print $1}' Orthogroups.GeneCount.tsv > P5
awk '{if($36+$37+$38+$39+$40+$41+$42>0)print $1}' Orthogroups.GeneCount.tsv > P6
awk '{if($43+$44+$45+$46>0)print $1}' Orthogroups.GeneCount.tsv > P7
awk '{if($47+$48+$49+$50+$51>0)print $1}' Orthogroups.GeneCount.tsv > P8
awk '{if($52+$53+$54+$55+$57+$58>0)print $1}' Orthogroups.GeneCount.tsv > P9

###
cat Alu_C Alu_D |sort |uniq > Alu
cat Bam_* |sort |uniq > Bam
cat Dsi_* |sort |uniq > Dsi
cat Gan_* |sort |uniq > Gan
cat Hca_* |sort |uniq > Hca
cat Mhu_* |sort |uniq > Mhu
cat Ogl_* |sort |uniq > Ogl
cat Ped_C Ped_D |sort |uniq > Ped
cat Rhi_* |sort |uniq > Rhi

###gene age
cat Alu Bam Dsi Gan Hca Mhu Ogl Ped Rhi |sort |uniq  > WB_all
cat  WB_all P1 |sort |uniq -d |wc -l

cat  WB_all P1 |sort |uniq -d > WB_PS/PS1 
cat WB_all WB_PS/PS1 |sort |uniq -u > WB_PS/PS1_left

cat WB_PS/PS1_left P2 |sort |uniq -d > WB_PS/PS2
cat WB_PS/PS1_left  WB_PS/PS2 |sort |uniq -u > WB_PS/PS2_left

cat WB_PS/PS2_left P3 |sort |uniq -d > WB_PS/PS3
cat WB_PS/PS2_left  WB_PS/PS3 |sort |uniq -u > WB_PS/PS3_left

cat WB_PS/PS3_left P4 |sort |uniq -d > WB_PS/PS4
cat WB_PS/PS3_left  WB_PS/PS4 |sort |uniq -u > WB_PS/PS4_left

cat WB_PS/PS4_left P5 |sort |uniq -d > WB_PS/PS5
cat WB_PS/PS4_left  WB_PS/PS5 |sort |uniq -u > WB_PS/PS5_left

cat WB_PS/PS5_left P6 |sort |uniq -d > WB_PS/PS6
cat WB_PS/PS5_left  WB_PS/PS6 |sort |uniq -u > WB_PS/PS6_left

cat WB_PS/PS6_left P7 |sort |uniq -d > WB_PS/PS7
cat WB_PS/PS6_left  WB_PS/PS7 |sort |uniq -u > WB_PS/PS7_left

cat WB_PS/PS7_left P8 |sort |uniq -d > WB_PS/PS8
cat WB_PS/PS7_left  WB_PS/PS8 |sort |uniq -u > WB_PS/PS8_left

cat WB_PS/PS8_left P9 |sort |uniq -d > WB_PS/PS9
cat WB_PS/PS8_left  WB_PS/PS9 |sort |uniq -u > WB_PS/PS9_left

cat WB_PS/PS9_left P10 |sort |uniq -d > WB_PS/PS10
cat WB_PS/PS9_left  WB_PS/PS10 |sort |uniq -u > WB_PS/PS10_left

cat Alu Bam Dsi Gan Hca Mhu Ogl Ped Rhi |sort |uniq -d > P11_WB2sp
cat WB_PS/PS10_left P11_WB2sp |sort |uniq -d > WB_PS/PS11
cat WB_PS/PS10_left  WB_PS/PS11 |sort |uniq -u > WB_PS/PS12


###WB new genes
cat Alu Hca Ped |sort |uniq > WB_PS/TWB_cat
cat Bam Dsi Mhu |sort |uniq > WB_PS/PWB_cat
cat Gan Ogl Rhi |sort |uniq > WB_PS/NWB_cat

cat TWB_cat PWB_cat NWB_cat |sort |uniq -c |awk '{if($1==3)print $2}' > 3lineage_comm
cat 3lineage_comm PS11 |sort |uniq -d > WB_new_gene_families

