#!/usr/bin/perl;
use strict;
use warnings;
use Bio::SeqIO;

my $usage = "\n\tUSAGE: perl $0 [fasta] [length]\n\n";
my $fasta = $ARGV[0];
my $length = $ARGV[1];
die $usage unless scalar(@ARGV) == 2;
die $usage unless -e $fasta;
die $usage unless $length > 0;

my $outname = "$fasta";
$outname =~ s/\..+$//;
$outname .= '_min_'.$length.'bp';
$outname =~ s/\..+$/.fa/;
$outname .= '.fa';

my $seqin = Bio::SeqIO->new(-file => $fasta,
                            -format => 'fasta');

my $seqout = Bio::SeqIO->new(-file => ">$outname",
                            -format => 'fasta');

while(my $seq = $seqin->next_seq) {
  if($seq->length >= $length) {
    $seqout->write_seq($seq);
  }
}
