sp=$1
a=$2
b=$3
python -m jcvi.compara.catalog ortholog --quota=$a:$b $sp Osa --no_strip_names
python -m jcvi.compara.synteny mcscan --iter=1 $sp.bed $sp.Osa.anchors -o $sp.Osa.blocks
cat $sp.Osa.blocks|sort -k2|groupBy -g 2 -c 1,1 -o count,distinct|awk -v num=$a '{ if( $2==num) print $0}'|cut -f 1,3 |sed -e 's/,/\t/g' >$sp.block
#after running for all species, use the following command to reconstruct "perfect-copy" syntenic gene data set for all species
#python -m jcvi.formats.base join $sp.block >perfect-copy_syntenic_gene.txt