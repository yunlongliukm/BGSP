
#Identification of positive selected genes using (I) PWB-subA as exapmle

##Identificaton of single genes
orthofinder -f /home/jinguihua/11BG_Analysis/kaks_A  -a 90

##CDS sequence extration of single copy genes
mkdir single_all 
cd single_all
cat Orthogroups_SingleCopyOrthologues.txt |while read i; do mv ../../Orthogroup_Sequences/${i}.fa ./ ;done

##Translate CDS to protein
/home/jinguihua/biosoft/cds2prot.pl 
#Usage: perl cds2prot.pl cds.fa > prot.fa
cat Orthogroups_SingleCopyOrthologues.txt |while read i; do perl cds2prot.pl ${i}.fa > ${i}.pep;done

##Alignment of protein sequences
cat Orthogroups_SingleCopyOrthologues.txt |while read i; do mafft --auto --thread 30 ${i}.pep > ${i}.aln.pep;done

##The coding sequences were aligned by PAL2NAL70 from the corresponding aligned protein sequences.
cat Orthogroups_SingleCopyOrthologues.txt |while read i; do perl pal2nal.pl ${i}.aln.pep ${i}.fa -output fasta > ./single_anligment/${i}.fa ;done

#rename genes 
cd ./single_anligment
cat Orthogroups_SingleCopyOrthologues.txt |while read i; do sed 's/Bam/bamA|Bam/g' ${i}.fa |sed 's/Dsi/dsiA|Dsi/g' |sed  's/Mhu/mhuA|Mhu/g' |sed  's/Ol/ola|Ol/g' |sed  's/LOC/osa|LOC/g' |sed  's/Chr/osa|Chr/g' |sed  's/Rgu/rgu|Rgu/g' > ${i}.fa1;done

rm *.fa
rename .fa1 .fa *

##caculate ka/ks with callCodeml2.py (https://github.com/byemaxx/callCodeml)
##singall.tre: ((mhuA,(dsiA,bamA))#1,(osa,(ola,rgu)));
 python3 single_anligment singall.tre

##positive selected genes with P value <= 0.05
awk '{if($7<=0.05)print $1}' result.tsv > PWB_subA_PSG_OG

