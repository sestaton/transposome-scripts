#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use File::Basename;
use Transposome::Run::Blast;
use Transposome::PairFinder;
use Getopt::Long;

my $infile;
my $format;
my $threads;
my $seqnum;
my $cpus;
my $out_dir;
my $help;

GetOptions(
	   'i|infile=s'   => \$infile,
	   'f|format=s'   => \$format,
           't|threads=i'  => \$threads,
           'n|seqnum=i'   => \$seqnum,
           'c|cpus=i'     => \$cpus,
	   'o|outdir=s'   => \$out_dir,
	   'h|help'       => \$help,
	   );

usage() and exit(0) if $help;
usage() and exit(1) if !$infile or !$out_dir;
my ($iname, $ipath, $isuffix) = fileparse($infile, qr/\.[^.]*/);

$format  //= 'fasta';
$threads //= 4;
$seqnum  //= 25_000;
$cpus    //= 2;

my $rep = $iname."_blast_report.txt";
my $blast = Transposome::Run::Blast->new( file      => $infile,
					  format    => $format,
                                          dir       => $out_dir,
                                          threads   => $threads,
                                          cpus      => $cpus,
                                          seq_num   => $seqnum,
                                          report    => $rep );

my $blastdb = $blast->run_allvall_blast;

for my $pid (qw(90.0 85.0 80.0)) {
    for my $fcov (qw(0.60 0.55 0.50)) {
        my $samp_out = $out_dir."/".$iname."_PID$pid"."_COV$fcov"."_out";

	say "Calculating graph edges from pairwise matches at $pid percent identity and $fcov coverage.";

        my $blast_res = Transposome::PairFinder->new( file              => $blastdb,  
                                                      dir               => $samp_out,
                                                      in_memory         => 1,
                                                      percent_identity  => $pid,
                                                      fraction_coverage => $fcov );

        my ($idx_file, $int_file, $hs_file) = $blast_res->parse_blast;
    }
}

exit;
## methods
sub usage {
    my $script = basename($0);
    print STDERR <<END

USAGE: $script [-i] [-n] [-h]

Required:
 -i|infile           :       The sequence file to sample.
 -o|outdir           :       The name of a directory to store the results
                             of each parameter set.

Options:
 -n|seqnum           :       The number of sequences to default (Default: 25000)
 -t|threads          :       The number of threads for each simulation (Default: 4).
 -c|cpus             :       The number of CPUs to use for each thread (Default: 2).
 -f|format           :       The input sequence format (Default: FASTA).
 -h|help             :       Print a usage statement.

END
}
