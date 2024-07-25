#!/usr/bin/perl
use strict;
use warnings;

open IN, "zcat $ARGV[0]|" or die $!;
#my $sample_dir = $ARGV[1];
#my $name = substr($ARGV[0],0,length($ARGV[0])-4);
my @name = split /.bam./, $ARGV[0];
open OUT, "|gzip - > ./$name[0]\.reads.split.bc.tsv.gz";

while(<IN>){
#  next if m/^@/;
	my @tmp = split/\s+/, $_;
#   next if $tmp[2] eq '*';
	my $chr = $tmp[0];
	my $start = $tmp[1];
	my $end = $tmp[2];
	my @sp = split/:/, $tmp[3]; 
#	my $l = length($sp);
#	my $read_id = substr($sp, 0, $l-18);
	my $cell_id = $sp[1].":".$sp[2].":".$sp[3]; 
	my $bc4 = $sp[4];
	my $umi	= $sp[5];
#	my $sample = "sample_".$bc4;
#	my $lib_type = substr($name[1],0,1)."NA";
#                  foo<-foo[,.(chr,start,end,read.id,bc4,umi,cellid,sub.lib, sample,lib.type)]
	print OUT "$chr\t$start\t$end\t$umi\t$cell_id\t$bc4\n";
}	
close IN;
close OUT;
#


