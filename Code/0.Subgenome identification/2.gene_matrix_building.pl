#script from Ke-Cheng Qian
use strict;
no strict 'refs';
#system ("sed -i 's/^M//g' *");
mkdir "dataset1" if(! -d "dataset1");
mkdir "dataset2" if(! -d "dataset2");
my @fasta_files=glob("*.fas");
my %gene_name_hash_seq;
foreach my $fasta_file(@fasta_files)
{
	open(FASTA,"$fasta_file")||die "$!";
	my @name_seq_hash=<FASTA>;
	for(my $i=0;$i<$#name_seq_hash;$i++)
	{
		my $j=$i+1;
		chomp $name_seq_hash[$i];
		chomp $name_seq_hash[$j];
		if ($name_seq_hash[$i]=~/>(.*-gene\d+)/)
		{
			my $gene_id="$1";
			#print "$gene_id\n";
			$gene_name_hash_seq{$gene_id}=$name_seq_hash[$j];
		}
	}
	close FASTA;
}
my $gene_num=0;
foreach my $gene_id(sort keys %gene_name_hash_seq)
{
	my @spls=split(/-/,$gene_id);
	$spls[1]=~/gene(\d+)/;
	my $gene_mark=$1;
	if ($gene_mark>$gene_num)
	{
		$gene_num=$gene_mark;
	}
}
for(my $i=1;$i<=$gene_num;$i++)
{
	my $marker_for_dataset1="gene$i";
	open(GENEOUT,">$marker_for_dataset1\_cds.fas")||die "$!";
	foreach my $gene_id(sort keys %gene_name_hash_seq)
	{
		$gene_id=~/(gene\d+)/;
		my $gene_eee=$1;
		if($gene_eee eq $marker_for_dataset1)
		{
			print GENEOUT ">$gene_id\n$gene_name_hash_seq{$gene_id}\n";
		}
	}
	close GENEOUT;
}
system("mv gene*_cds.fas dataset1/");


my %dataset2_seq;
for(my $i=1;$i<=$gene_num;$i++)
{
	foreach my $gene_id(sort keys %gene_name_hash_seq)
	{
		my @spls=split(/-/,$gene_id);
		my $script_name;
		$script_name=$spls[0];#
		$spls[1]=~/gene(\d+)/;
		my $gene_mark=$1;
		if($gene_mark==$i)
		{
			$dataset2_seq{$script_name} .="$gene_name_hash_seq{$gene_id}";
		}
	}
}
open(SPEOUT,">super_matrix_dataset_cds.fas")||die "$!";
foreach my $matrix(sort keys %dataset2_seq)
{
	print SPEOUT ">$matrix\n$dataset2_seq{$matrix}\n";
}
close SPEOUT;
system("mv super_matrix_dataset_cds.fas ./dataset2/");
