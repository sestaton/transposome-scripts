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
my $out_dir;
my $help;

GetOptions(
	   'i|infile=s'   => \$infile,
	   'f|format=s'   => \$format,
	   'o|outdir=s'   => \$out_dir,
	   'h|help'       => \$help,
	   );

usage() and exit(0) if $help;
usage() and exit(1) if !$infile or !$out_dir;
my ($iname, $ipath, $isuffix) = fileparse($infile, qr/\.[^.]*/);

$format //= 'fasta';

my $rep = $iname."_blast_report.txt";
my $blast = Transposome::Run::Blast->new( file      => $infile,
					  format    => $format,
                                          dir       => $out_dir,
                                          threads   => 4,
                                          cpus      => 2,
                                          seq_num   => 25_000,
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
# methods
sub usage {
    my $script = basename($0);
    print STDERR <<END
USAGE: $script [-i] [-n] [-h]

Required:
 -i|infile           :       The sequence file to sample.
 -o|outdir           :       The name of a directory to store the results
                             of each parameter set.

Options:
 -f|format           :       The input sequence format (Default: FASTA).
 -h|help             :       Print a usage statement.

END
}
