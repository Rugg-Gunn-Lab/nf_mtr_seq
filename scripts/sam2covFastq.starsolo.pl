#!/usr/bin/perl
use strict;
use warnings;

open IN, "$ARGV[0]" or die $!;

my $name = $ARGV[0];
$name =~ s/\.[^.]+$//;

open OUT1, "|gzip - > ${name}_cov.R1.fq.gz" or die $!;
open OUT2, "|gzip - > ${name}_cov.R2.fq.gz" or die $!;

while(<IN>){

  chomp;
  next if m/^@/;

  my @tmp = split /\t/, $_;

  # SAM columns:
  # tmp[0] = read name
  # tmp[2] = reference / cell ID in your custom file
  # tmp[9] = query sequence
  # tmp[10] = query quality

  my $cell_id;
  if ($tmp[2] eq '*') {
    $cell_id = "unmapped";
  } else {
    $cell_id = $tmp[2];
  }

  my @sp = split /:/, $tmp[0];

  # Assume sp[7] is UMI
  my $umi = $sp[7];

  # Short read name to avoid SAM/BAM query-name-too-long problem
  my $readname = "@".$sp[0].".".$sp[3].".".$sp[4].".".$sp[5].".".$sp[6].":".$cell_id.":".$sp[7];

  # Original R1 sequence from encoded read name
  my $read1 = $sp[8];
  my $R1len = length($read1);

  my $R2len = $ARGV[1] - 146;

  # Original R1/R2 qualities and R2 sequence encoded in read name
  my $qual1 = substr($tmp[0], length($tmp[0]) - $R2len - 1 - $R2len - 1 - $R1len, $R1len);
  my $read2 = substr($tmp[0], length($tmp[0]) - $R2len - 1 - $R2len, $R2len);
  my $qual2 = substr($tmp[0], length($tmp[0]) - $R2len, $R2len);

  # Query sequence and quality from SAM alignment
  my $query_seq  = $tmp[9];
  my $query_qual = $tmp[10];

  if ($query_seq eq '*') {
    $query_seq = "";
  }

  if ($query_qual eq '*') {
    $query_qual = "";
  }

  # UMI quality is not stored separately, so use high-quality placeholder
  my $umi_qual = "I" x length($umi);

  # Put query sequence + UMI at the front of read2
  # This prepares R2 as: query_seq + UMI + original_read2
  $read2 = $query_seq . $umi . $read2;
  $qual2 = $query_qual . $umi_qual . $qual2;

  my $mark = "+";

  print OUT1 "$readname\n$read1\n$mark\n$qual1\n";
  print OUT2 "$readname\n$read2\n$mark\n$qual2\n";
}

close IN;
close OUT1;
close OUT2;
