#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use File::Basename;
use Transposome::SeqUtil;
use Getopt::Long;

my $infile;
my $outfile;
my $format;
my $list;
my $help;

GetOptions(
	   'i|infile=s'  => \$infile,
	   'f|format=s'  => \$format,
	   'o|outfile=s' => \$outfile,
	   'l|list=s'    => \$list,
	   'h|help'      => \$help,
	   );

usage() and exit(0) if $help;
usage() and exit(1) if !$infile || !$outfile || !$list;

$format //= 'fasta';

my $seq_obj = Transposome::SeqUtil->new( file      => $infile,
					 format    => $format,
					 in_memory => 0 ); 

my ($seqs, $seqct) = $seq_obj->store_seq;

say "There are $seqct sequences in $infile.";

open my $out, '>', $outfile;
open my $l, '<', $list;
my %ids = map { $_ => 1 } <$l>;
close $l;

for my $id (keys %ids) {
    if (exists $seqs->{$id}) {
	say $out join "\n", ">$id", $seqs->{$id};
    }
}
close $out;

exit;
# methods
sub usage {
    my $script = basename($0);
    print STDERR <<END
USAGE: $script [-i] [-o] [-l] [-h]

Required:
-i|infile           :       The sequence file to pull sequences from.
-o|outfile          :       The file to write the extracted sequences to.
-l|list             :       A file containing IDs (one per line) that should be written
                            to the output file.

Options:
-f|format           :       The input sequence format (Default: FASTA).
-h|help             :       Print a usage statement.

END
}
