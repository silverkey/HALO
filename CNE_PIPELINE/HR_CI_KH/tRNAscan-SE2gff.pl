#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

use lib '/home/remo/src/bioperl-live';
use lib '.';
use tRNAscanSE;

my $usage = "\n\tperl $0 [fasta]\n\n";
die $usage unless scalar(@ARGV) == 1;
die $usage unless -e $ARGV[0];

my $file = $ARGV[0];

my $parser = tRNAscanSE->new(-file => $file);

# parse the results
while(my $gene = $parser->next_prediction) {
  $gene->gff_format(Bio::Tools::GFF->new(-gff_version => 3));
  print $gene->gff_string;
  print "\n";
  foreach my $exon($gene->get_SeqFeatures) {
    $exon->gff_format(Bio::Tools::GFF->new(-gff_version => 3));
    print $exon->gff_string;
    print "\n";
  }
}
