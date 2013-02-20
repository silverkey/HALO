#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Bio::Tools::GFF;
use Bio::SearchIO;

my $USAGE = "\nperl $0 [blast output file in format 0] [source] [organism1] [organism2]\n\n".
            "N.B.: organism1 is the query organism2 is the target\n\n";

die $USAGE unless $ARGV[0];
die $USAGE unless -e $ARGV[0];
die $USAGE unless $ARGV[1];
die $USAGE unless $ARGV[2];

my $ID = 1;
my $blast_file = $ARGV[0];
my $source = $ARGV[1];
my $organism1 = $ARGV[2];
my $organism2 = $ARGV[3];

my $gff_file1 = "$blast_file";
$gff_file1 =~ s/\.[\w\d\_]+$//;
$gff_file1 .= "\_$organism1\.gff";

my $gff_file2 = "$blast_file";
$gff_file2 =~ s/\.[\w\d\_]+$//;
$gff_file2 .= "\_$organism2\.gff";

open(OUT1,">$gff_file1");
open(OUT2,">$gff_file2");

my $in = Bio::SearchIO->new(-file => $blast_file,
                            -format => 'blast',
                            -report_format => 0);

while(my $res = $in->next_result) {
  while(my $hit = $res->next_hit) {
    while(my $hsp = $hit->next_hsp) {

      my $f1 = $hsp->feature1;
      my $f2 = $hsp->feature2;

      my $type = 'match';

      my $score = $hsp->percent_identity;
      $score =~ s/^(\d+\.\d)\d+$/$1/;

      my $strand1 = $f1->strand;
      $strand1 = '+' if $strand1 eq '1';
      $strand1 = '-' if $strand1 eq '-1';

      my $strand2 = $f2->strand;
      $strand2 = '+' if $strand2 eq '1';
      $strand2 = '-' if $strand2 eq '-1';

      my $phase = '0';

      # The string of the attributes can be defined
      # ID => UNIQUE ID FOR THE FEATURE (predefined gff3 attribute)
      # Target => TARGET  (predefined gff3 attribute)
      # E => E-VALUE
      # C => CIGAR STRING
      # O => ORGANISM

      my $attribute1 = 'ID='.$ID.';'.'Target='.$f2->seq_id.' '.$f2->start.' '.$f2->end.' '.$strand2.';'.'E='.$hsp->evalue.';'.'O='.$organism2.';'.'C='.$hsp->cigar_string;
      print OUT1 $f1->seq_id."\t".$source."\t".$type."\t".$f1->start."\t".$f1->end."\t".$score."\t".$strand1."\t".$phase."\t".$attribute1."\n";

      my $attribute2 = 'ID='.$ID.';'.'Target='.$f1->seq_id.' '.$f1->start.' '.$f1->end.' '.$strand1.';'.'E='.$hsp->evalue.';'.'O='.$organism1.';'.'C='.$hsp->cigar_string;
      print OUT2 $f2->seq_id."\t".$source."\t".$type."\t".$f2->start."\t".$f2->end."\t".$score."\t".$strand2."\t".$phase."\t".$attribute2."\n";

      $ID ++;
    }
  }
}

close(OUT1);
close(OUT2);
