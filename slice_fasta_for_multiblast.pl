#!/usr/bin/perl;
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Seq;

my $usage = "\n\tUSAGE: perl $0 [fasta] [processes]\n\n";
my $fasta = $ARGV[0];
my $processes = $ARGV[1];
die $usage unless scalar(@ARGV) == 2;
die $usage unless $processes > 0;

my $outname = "$fasta";
$outname =~ s/\..+$//;
$outname .= '_slice';

my $total_bases = 0;
calculate_total_bases();
my $bases_per_slice = int($total_bases/$processes)-1;

print "\nTotal bases in $fasta: $total_bases\n";
print "The average number of bases in the $processes slices will be: $bases_per_slice\n\n";

my $seqin = Bio::SeqIO->new(-file => $fasta,
                            -format => 'fasta');

my $written = 0;
my $c = 1;

my $name = "$outname\_$c.fasta";
my $seqout;
$seqout = Bio::SeqIO->new(-file => ">$name",
                          -format => 'fasta');

while(my $seq = $seqin->next_seq) {
  reset_seqout() if $written > $bases_per_slice;
  $seqout->write_seq($seq);
  $written += $seq->length;
}

sub calculate_total_bases {
  my $seqin = Bio::SeqIO->new(-file => $fasta,
                              -format => 'fasta');
  while(my $seq = $seqin->next_seq) {
    $total_bases += $seq->length;
  }
  $seqin->close;
}

sub reset_seqout {
  $c++;
  $written = 0;
  $seqout->close;
  my $name = "$outname\_$c.fasta";
  $seqout = Bio::SeqIO->new(-file => ">$name",
                            -format => 'fasta');
}
