#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;

# VERSION: 2

my $processes = 10;
my $cmd = '/root/src/ncbi-blast-2.2.27+/bin/blastn ';
$cmd .= '-db /root/BLASTDB/CI_KH_rm_1000bp ';
$cmd .= '-task blastn -evalue 0.01 -word_size 4 -culling_limit 1 -num_threads 1 ';

open(LOG,'>MULTIBLAST.LOG');

my $pm = new Parallel::ForkManager($processes);

$pm->run_on_finish(
  sub {
    my($pid,$exit_code) = @_;
    print "** Just got out of the pool ".
          "with PID $pid and exit code: $exit_code\n";
  }
);

my @fasta = glob('*.fasta');
foreach my $fasta (@fasta) {
  my $out = "$fasta";
  $out =~ s/.fasta$/.blast_out/;
  my $run = $cmd.'-query '.$fasta.' -out '.$out;
  sleep 1;
  # Forks and returns the pid for the child:
  my $pid = $pm->start and next;

  # Here is the parallelized block
  # -----------
  print "$pid ---> running on $fasta\n";
  print LOG "Launching $run\n";
  system($run);
  # -----------

  # Terminates the child process
  $pm->finish;
}
