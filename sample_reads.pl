#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use File::Basename;
use Transposome::SeqUtil;
use Getopt::Long;

my $infile;
my $format;
my $sample_size;
my $help;

GetOptions(
	   'i|infile=s'      => \$infile,
	   'f|format=s'      => \$format,
	   'n|sample_size=i' => \$sample_size,
	   'h|help'          => \$help,
	   );

usage() and exit(0) if $help;
usage() and exit(1) if !$infile or !$sample_size;
my ($iname, $ipath, $isuffix) = fileparse($infile, qr/\.[^.]*/);

$format //= 'fasta';

say STDERR "Sampling $infile for sample size: ", $sample_size;
my $sequtil = Transposome::SeqUtil->new( file        => $infile, 
					 format      => $format,
					 sample_size => $sample_size, 
					 no_store    => 1 );

my $out = $iname."_$sample_size".".fasta";

{ 
    local *STDOUT; 
    open STDOUT, '>', $out; 
    $sequtil->sample_seq; 
}

exit;

# methods
sub usage {
    my $script = basename($0);
    print STDERR <<END
USAGE: $script -i seqs.fasta -n 100000 

Required:
 -i|infile           :       The sequence file to sample.
 -n|sample_size      :       The number of sequences to sample (integer).
                             (NB: the value may be as '1000000' or with underscores
                              as in '1_000_000').
Options:
 -f|format            :       The input sequence format (Default: FASTA).
 -h|help              :       Print a usage statement.

END
}
