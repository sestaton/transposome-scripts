#!/usr/bin/env perl

use 5.012;
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Transposome::Cluster;
use Transposome::SeqUtil;

#my $usage = "$0 -s seqfile.fas -int int_file -idx idx_file -o output_directory\n";
my $report = 'cluster_merge_permute_report.txt';
my $int_file;
my $idx_file;
my $seq_file;
my $out_dir;
my $help;

GetOptions( 
	   's|seqfile=s'     => \$seq_file,
           'int|intfile=s'   => \$int_file,
           'idx|idxfile=s'   => \$idx_file,
           'o|outdir=s'      => \$out_dir,
	   'h|help'          => \$help,
	  );

usage() and exit(1) if $help;
usage() and exit(1) if !$int_file or !$idx_file or !$seq_file or !$out_dir;

for my $merge_thresh (qw(100 500 1000)) { # the exact threshold will depend on the size of the data set
    my $cluster = Transposome::Cluster->new( file            => $int_file, 
                                             dir             => $out_dir,
                                             merge_threshold => $merge_thresh,
                                             cluster_size    => 10 );

    my $comm = $cluster->louvain_method;
    my $cluster_file = $cluster->make_clusters($comm, $idx_file);

    my ($read_pairs, $vertex, $uf) = $cluster->find_pairs($cluster_file, $report);
    my $memstore = Transposome::SeqUtil->new( file => $seq_file, in_memory => 1 );
    my ($seqs, $seqct) = $memstore->store_seq;

    my ($cls_dir_path, $cls_with_merges_path, $cls_tot) = $cluster->merge_clusters($vertex, $seqs, $read_pairs, $report, $uf);
}

exit;
# methods
sub usage {
    my $script = basename($0);
    print STDERR <<END
USAGE: $script [-i] [-o] [-int] [-idx] [-h]

Required:
-s|seqfile            :       The sequence file to pull sequences from.
-o|outdir             :       A directory name for the results.
-int|intfile          :       The integer mapping file of sequence IDs generated from the 'louvain_method' method
                              of the Transposome::Cluster class.
-idx|idxfile          :       The index file of matches generated from the 'louvain_method' method
                              of the Transposome::Cluster class.

Options:
-h|help             :       Print a usage statement.

END
}
