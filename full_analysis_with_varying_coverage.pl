#!/usr/bin/env perl

## pragmas and library imports
use 5.012;
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Transposome::PairFinder;
use Transposome::Cluster;
use Transposome::SeqUtil;
use Transposome::Annotation;
use Transposome::Run::Blast;

## lexical vars
my $seq_file;
my $seq_format;
my $out_dir;
my $rep_db;
my $cpu;
my $thread;
my $seq_num;
my $cls_size;
my $bls_pid;
my $bls_fcov;
my $help;

GetOptions( 's|sequence_file=s'        => \$seq_file,
	    'f|sequence_format=s'      => \$seq_format,
            'o|output_dir=s'           => \$out_dir,
            'repdb|repeat_database=s'  => \$rep_db,
            'c|cpus=i'                 => \$cpu,
            't|threads=i'              => \$thread,
            'n|sequence_num=i'         => \$seq_num,
            'pid|percent_identity=i'   => \$bls_pid,
            'fcov|fraction_coverage=f' => \$bls_fcov,
            'cls|cluster_size=i'       => \$cls_size,
            'h|help'                   => \$help,
            );

## check input
usage() and exit(0) if $help;

if (!$seq_file || !$out_dir || !$cpu || !$thread) {
    say "\n[ERROR]: Command line not parsed correctly. Check input and refer to the usage. Exiting.\n";
    usage();
    exit(1);
}


$seq_num  //= 25000;
$cls_size //= 100;
$bls_pid  //= 90;
$bls_fcov //= 0.55;
$seq_format //= 'fasta';
my $merge_thresh = 100;
my ($iname, $ipath, $isuffix) = fileparse($seq_file, qr/\.[^.]*/);

for my $sample_size (qw(100_000 200_000 300_000 400_000 500_000)) {
    my $sequtil = Transposome::SeqUtil->new( file        => $seq_file,
					     format      => $seq_format,
					     sample_size => $sample_size, 
                                             no_store    => 1 );

    my $samp_seq = $iname."_$sample_size".".fasta";

    { 
	local *STDOUT; 
	open STDOUT, '>', $samp_seq; 
	$sequtil->sample_seq; 
    }

    my $samp_rep = $iname."_$sample_size"."_report.txt";
    my $samp_dir = $out_dir."_$sample_size";

    # perform the all vs. all blast
    my $blast = Transposome::Run::Blast->new( file      => $samp_seq,
                                              dir       => $samp_dir,
                                              threads   => $thread,
                                              cpus      => $cpu,
                                              seq_num   => $seq_num,
                                              report    => $samp_rep );

    my $blastdb = $blast->run_allvall_blast;

    # parse mglblast results to find best scoring pairwise matches
    my $blast_res = Transposome::PairFinder->new( file              => $blastdb,  
                                                  dir               => $samp_dir,
                                                  in_memory         => 1,
                                                  percent_identity  => $bls_pid,
                                                  fraction_coverage => $bls_fcov );

    my ($idx_file, $int_file, $hs_file) = $blast_res->parse_blast;

    # cluster sequences and analyze groupings
    my $cluster = Transposome::Cluster->new( file            => $int_file,
                                             dir             => $samp_dir,
                                             merge_threshold => $merge_thresh,
                                             cluster_size    => $cls_size );

    my $comm = $cluster->louvain_method;
    my $cluster_file = $cluster->make_clusters($comm, $idx_file);
    my ($read_pairs, $vertex, $uf) = $cluster->find_pairs($cluster_file, $samp_rep);
    my $memstore = Transposome::SeqUtil->new( file => $samp_seq, in_memory => 1 );
    my ($seqs, $seqct) = $memstore->store_seq;

    my ($cls_dir_path, $cls_with_merges_path, $cls_tot) 
       = $cluster->merge_clusters($vertex, $seqs, $read_pairs, $samp_rep, $uf);

    # annotate clusters and generate whole-genome summary of results
    my $annotation = Transposome::Annotation->new( database  => $rep_db,
                                                   dir       => $samp_dir,
                                                   file      => $samp_rep );

    my ($anno_rp_path, $anno_sum_rep_path, $total_readct, $rep_frac, $blasts, $superfams) 
       = $annotation->annotate_clusters($cls_dir_path, $seqct, $cls_tot);

    $annotation->clusters_annotation_to_summary($anno_rp_path, $anno_sum_rep_path, $total_readct,
                                                $seqct, $rep_frac, $blasts, $superfams, $samp_rep);

    $merge_thresh += 100;
}
exit;

## methods
sub usage {
    my $script = basename($0);
    print STDERR <<END
USAGE: $script -s seqs.fastq -o my_cluster_report.txt -n 25000 -c 2 -t 12 [-h]

Required:
 -s|sequence_file       :       A Fasta/q file to analyze for repeats
 -o|output_dir          :       A name for the output to be written.
 -c|cpus                :       The number of  CPUs to use for each blast process.
 -t|threads             :       The number of parallel blast processes to run.
                                (NB: threads X cpus should be less than or equal to the number
                                of CPUs on your machine.)
 -repdb|repeat_database :       A sequence file of repeats to be used for annotation.

Options:
-f|sequence_format      :       The input sequence format (Default: FASTA).
-pid|percent_identity   :       Percent identity between pairwise matches in all vs. all blast (Default: 90).
-fcov|fraction_coverage :       The fraction coverage between pairwise matches in all vs. all blast (Default: 0.55).
-cls|cluster_size       :       The minimum size of a cluster to be used for annotation (Default: 100).
-n|sequence_num         :       The number of sequences to submit to each blast process (Default: 25000).
-h|help                 :       Print a usage statement.

END
}
