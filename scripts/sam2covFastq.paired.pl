#!/usr/bin/perl
use strict;
use warnings;

open IN, "$ARGV[0]" or die $!;
my $name = substr($ARGV[0],0,length($ARGV[0])-4);
open OUT1, "|gzip - > $name\_cov.R1.fq.gz";
open OUT2, "|gzip - > $name\_cov.R2.fq.gz";

while(<IN>){
  next if m/^@/;
  my @tmp = split/\s+/, $_;
  next if $tmp[2] eq '*';
  my $cell_id = $tmp[2];
  my @sp = split/:/, $tmp[0];
 # my $readname = "@".$sp[0].":".$cell_id.":".$sp[1];
  my $readname = "@".$sp[0].".".$sp[3].".".$sp[4].".".$sp[5].".".$sp[6].":".$cell_id.":".$sp[7];
  my $read1 = $sp[8];
  my $R1len = length($sp[8]);
  my $R2len = $ARGV[1]-146;
  my $qual1 = substr($tmp[0],length($tmp[0])-$R2len-1-$R2len-1-$R1len,$R1len);
  my $read2 = substr($tmp[0],length($tmp[0])-$R2len-1-$R2len,$R2len);
  my $qual2 = substr($tmp[0],length($tmp[0])-$R2len,$R2len);
#  my $l = length($read);
  my $mark = "+";
#  my $qual = substr($tmp[0], -$l, $l);
  print OUT1 "$readname\n$read1\n$mark\n$qual1\n";
  print OUT2 "$readname\n$read2\n$mark\n$qual2\n";

}
close IN;
close OUT1;
close OUT2;
