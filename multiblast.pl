#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;

my $cmd = '/home/remo/src/ncbi-blast-2.2.27+/bin/blastn ';
$cmd .= '-db ../GENOMES/MASKED/CI_KH_rm_1000bp ';
$cmd .= '-task blastn -evalue 0.01 -word_size 4 -culling_limit 1 -num_threads 1 ';

open(LOG,'>MULTIBLAST.LOG');

my $pm = new Parallel::ForkManager(20);

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


__END__
nohup /home/remo/src/ncbi-blast-2.2.27+/bin/blastn -query ../GENOMES/MASKED/HR_rm_1000bp.fa -db ../GENOMES/MASKED/CI_KH_rm_1000bp -task blastn -out HR_vs_CI_ghost_rm.out -evalue 
0.01 -word_size 4 -culling_limit 1 -num_threads 20 &

