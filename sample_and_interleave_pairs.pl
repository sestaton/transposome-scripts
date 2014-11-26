#!/usr/bin/env perl

use 5.012;
use strict;
use warnings;
use autodie qw(open);
use File::Basename;
use Getopt::Long;
use Transposome::SeqUtil;

my $pair1;
my $pair2;
my $sample_size;
my $outfile;
my $help;

GetOptions(
           'f|pair1=s'       => \$pair1,
           'r|pair2=s'       => \$pair2,
           'n|sample_size=i' => \$sample_size,
	   'o|outfile=s'     => \$outfile,
	   'h|help'          => \$help,
           );

usage() and exit(0) if $help;
usage() and exit(1) if !$pair1 || !$pair2 || !$sample_size || !$outfile;

open my $out, '>', $outfile;
my $sequtil_pair1 = Transposome::SeqUtil->new( file => $pair1, sample_size => $sample_size, seed => 100, in_memory => 1 );
my $sequtil_pair2 = Transposome::SeqUtil->new( file => $pair2, sample_size => $sample_size, seed => 100, in_memory => 1 );
my ($seqs1, $seqs1_ct) = $sequtil_pair1->sample_seq;
my ($seqs2, $seqs2_ct) = $sequtil_pair2->sample_seq;

for (my ($seqs1_id, $seqs1_seq) = each%$seqs1) { # each reduces the memory usage, compared to keys
    my $seqs2_id = $seqs1_id; 
    $seqs2_id =~ s/1$/2/;
    if (exists $seqs2->{$seqs2_id}) {
        say $out join "\n", ">".$seqs1_id, $seqs1_seq;
        say $out join "\n", ">".$seqs2_id, $seqs2->{$seqs2_id};
	delete $seqs1->{$seqs1_id};
	delete $seqs2->{$seqs2_id};
    }
}
close $out;

exit;
#methods
sub usage {
    my $script = basename($0);
    print STDERR <<END
USAGE: $script [-f] [-r] [-o] [-n] [-h]

Required:
-f|pair2           :       The forwared sequence file to sample.
-r|pair2           :       The reverse sequence file to sample.
-o|outfile         :       The number of sequences to sample (integer).
-n|sample_size     :       (NB: the value may be given as '1000000' or with underscores
                            as in '1_000_000').

Options:
-h|help            :       Print a usage statement.

END
}
