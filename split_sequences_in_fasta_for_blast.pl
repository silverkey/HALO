#!/usr/bin/perl;
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Seq;

my $usage = "\n\tUSAGE: perl $0 [fasta] [fragments length] [fragments overlap]\n\n";
my $fasta = $ARGV[0];
my $frag_length = $ARGV[1];
my $frag_overlap = $ARGV[2];
die $usage unless scalar(@ARGV) == 3;
die $usage unless -e $fasta;
die $usage unless $frag_length > 0;
die $usage unless $frag_overlap >= 0;

my $outname = "$fasta";
$outname =~ s/\..+$//;
$outname .= '_split_'.$frag_length;
$outname .= '_overlap_'.$frag_overlap;
$outname =~ s/\..+$/.fa/;
$outname .= '.fa';

my $seqin = Bio::SeqIO->new(-file => $fasta,
                            -format => 'fasta');

my $seqout = Bio::SeqIO->new(-file => ">$outname",
                            -format => 'fasta');

while(my $seq = $seqin->next_seq) {
  my $seq_length = $seq->length;
  my $start = 1;
  my $end = $frag_length;
  $end = $seq_length if $end > $seq_length;

  while($start < $seq_length) {
    my $fragid = $seq->id.'_'.$start.'_'.$end;
    my $fragseq = $seq->subseq($start,$end);
    my $fragdesc = "original seq length $seq_length";

    my $frag = Bio::Seq->new(-id => $fragid,
                             -seq => $fragseq,
                             -desc => $fragdesc);

    $seqout->write_seq($frag);

    $start = $start + $frag_length - $frag_overlap;
    $end = $start + $frag_length - 1;
    $end = $seq_length if $end > $seq_length;
  }
}
