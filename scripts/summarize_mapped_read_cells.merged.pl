#!/usr/bin/perl
use strict;
use warnings;

my $s = $ARGV[0];
my $genomes = $ARGV[1];
my $sampledir = $ARGV[2];
my $conditions = $ARGV[3];
my $fastq = "/bi/home/wangy/projects/ctr_seq/02.trimmed/$sampledir/$s\_BC_cov_trimmed.merged\.$conditions\.fq.gz";
my $bam = "$s\_$genomes\_sorted_merged\.$conditions\.bam";

my $out = "./$s\_$genomes\_merged\.$conditions\.summary.xls" or die $!;

my %hash_all;
my %hash_50;
my %unique_all;
my %unique_100_all;
my %unique_50;
my %unique_100_50;
my %n_reads_all;
my %n_reads_50;
my %n_mapped_all;
my %n_mapped_50;
my %n_uniquely_mapped_all;
my %n_uniquely_mapped_50;
my %n_fragments_all;
my %n_fragments_all_100;
my %n_fragments_50;
my %n_fragments_50_100;


open IN, "zcat $fastq|" or die $!;
while(<IN>){
	my $line1 = $_;
	my $line2 = <IN>;
	my $line3 = <IN>;
	my $line4 = <IN>;
	chomp($line1);
	chomp($line3);
	my $cell_id = substr($line1, -20, 11);
	my $umi = substr($line1, -8, 8);
	$hash_all{$cell_id}{"raw"} = 0 if not exists $hash_all{$cell_id}{"raw"};
	$hash_all{$cell_id}{"raw"}++;
	next if length($line2)<20;
	$hash_50{$cell_id}{"raw"} = 0 if not exists $hash_50{$cell_id}{"raw"};
	$hash_50{$cell_id}{"raw"}++;
}
close IN;


open IN, "samtools view $bam|" or die $!;
while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	my $umi = substr($tmp[0], -8, 9);
	my $cell_id = substr($tmp[0], -20, 11);
	my $chr = $tmp[2];
my $pos = $tmp[3];
	my $pos100 = int($pos/50);
	## all reads
	$hash_all{$cell_id}{"mapped"} = 0 if not exists $hash_all{$cell_id}{"mapped"};
	$hash_all{$cell_id}{"mapped"}++;
	if(substr($tmp[12], 0, 2) ne "XS"){
		$hash_all{$cell_id}{"unique_mapped"} = 0 if not exists $hash_all{$cell_id}{"unique_mapped"};
		$hash_all{$cell_id}{"unique_mapped"}++;
}
	my $umi_all = $umi."_".$chr."_".$pos;
	my $umi_100 = $umi."_".$chr."_".$pos100;
	$unique_all{$cell_id}{$umi_all} = 1;
	$unique_100_all{$cell_id}{$umi_100} = 1;
	my $length = length($tmp[9]);

	## reads < 20
	next if $length < 20;
	$hash_50{$cell_id}{"mapped"} = 0 if not exists $hash_50{$cell_id}{"mapped"};
	$hash_50{$cell_id}{"mapped"}++;
	if(substr($tmp[12], 0, 2) ne "XS"){
		$hash_50{$cell_id}{"unique_mapped"} = 0 if not exists $hash_50{$cell_id}{"unique_mapped"};
		$hash_50{$cell_id}{"unique_mapped"}++;
	}
	$unique_50{$cell_id}{$umi_all} = 1;
	$unique_100_50{$cell_id}{$umi_100} = 1;
}
close IN;

open OUT, ">$s\_$genomes\_merged\.$conditions\.summary.xls"  or die $!;
foreach my $cell_id (keys %hash_50){
	my $n_reads_all = $hash_all{$cell_id}{"raw"};
	my $n_reads_50 = $hash_50{$cell_id}{"raw"};
	my $n_mapped_all = $hash_all{$cell_id}{"mapped"};
	my $n_mapped_50 = $hash_50{$cell_id}{"mapped"};
	my $n_uniquely_mapped_all = $hash_all{$cell_id}{"unique_mapped"};
	my $n_uniquely_mapped_50 = $hash_50{$cell_id}{"unique_mapped"};
	my $n_fragments_all = keys %{$unique_all{$cell_id}};
	my $n_fragments_all_100 = keys %{$unique_100_all{$cell_id}};
	my $n_fragments_50 = keys %{$unique_50{$cell_id}};
	my $n_fragments_50_100 = keys %{$unique_100_50{$cell_id}};
	my $output = "$cell_id\t$n_reads_all\t$n_mapped_all\t$n_uniquely_mapped_all\t$n_fragments_all\t$n_fragments_all_100\t";
	   $output .= "$n_reads_50\t$n_mapped_50\t$n_uniquely_mapped_50\t$n_fragments_50\t$n_fragments_50_100\n";
## cellid, all reads, all mapped reads, unique mapped reads (mapping),unique reads (cellid, chr,umi,start position), unique reads (cellid, chr,umi, start position not in the same 50bp bin)
	print OUT $output;
}
close OUT;
