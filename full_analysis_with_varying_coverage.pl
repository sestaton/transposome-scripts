#!/usr/bin/env perl

use 5.014; # do not lower; required for non-destructive substitution
use strict;
use warnings;
use autodie;
use File::Basename;
use File::Find;
use File::Copy          qw(copy);
use File::Path          qw(make_path remove_tree);
use IPC::System::Simple qw(system);
use Cwd                 qw(getcwd);
use Regexp::Common;
use Sort::Naturally;
use Getopt::Long;

my %opts;

GetOptions(\%opts,
    'accession|id=s',
    'sequence_format|sf=s',  
    'output_dir|p=s',        
    'repeat_database|db=s',  
    'threads|t=i',           
    'percent_identity|pid=i',
    'coverage|cov=f',        
    'cluster_size|cls=i',    
    'genome_size|gs=i',
    'species|s=s',
    'h|help',                
);

## check input
usage() and exit(0) if $opts{help};

if (!$opts{accession} || !$opts{output_dir} || !$opts{repeat_database} || !$opts{genome_size}) {
    say "\n[ERROR]: Command line not parsed correctly. Check input and refer to the usage. Exiting.\n";
    usage();
    exit(1);
}

unless (-e $opts{output_dir}) { 
    make_path($opts{output_dir}, {verbose => 0, mode => 0711,});
}

my $sumdir = File::Spec->catdir($opts{output_dir}, 'annotation_summaries');
unless (-e $sumdir) {
    make_path($sumdir, {verbose => 0, mode => 0711,});
}

my $outfile = File::Spec->catfile($opts{output_dir}, $opts{accession}.'_all_transposome_results.tsv');
open my $out, '>>', $outfile;
say $out join "\t", 'x_coverage', 'number_clustered', 'number_unclustered', 'repeat_fraction', 
    'singleton_fraction', 'total_repeat_percentage', 'total_annotated_repeat_percentage', 
    'family_count', 'elapsed_time';

## set defaults
$opts{threads} //= 1;
$opts{cluster_size} //= 100;
$opts{percent_identity}  //= 90;
$opts{coverage} //= 0.55;
$opts{sequence_format} //= 'fastq';
my $merge_thresh = 100;

## try to get the sequence files
my ($forward, $reverse) = ($opts{accession}.'_1.fastq', $opts{accession}.'_2.fastq');
unless (-e $forward && -e $reverse) {
    my $fgzip = $opts{accession}.'_1.fastq.gz';
    my $rgzip = $opts{accession}.'_2.fastq.gz';
    my $fbzip = $opts{accession}.'_1.fastq.bz2';
    my $rbzip = $opts{accession}.'_2.fastq.bz2';

    $forward = -e $fgzip ? $fgzip : $fbzip;
    $reverse = -e $rgzip ? $rgzip : $rbzip;

    unless (-e $forward && -e $reverse) {
	say STDERR "\n[ERROR]: Could not find sequence files. Make sure they exist and end with '.fastq', '.fastq.gz' or '.fastq.bz2'. Exiting.\n";
	exit(1);
    }
}

my ($iname, $ipath, $isuffix) = fileparse($forward, qr/\.[^.]*/);
my $base = $iname =~ s/_1.*//r;

my %cvalues;
my $cvalue = $opts{genome_size} * 1e6;

for my $x (qw(0.01 0.03 0.05 0.07 0.09 0.10)) { 
    my $reads = $x*$cvalue / 100;
    my $bases = $x*$cvalue;
    my $pairs = $reads / 2;
    $cvalues{$x} =  { bases => $bases, reads => $reads, pairs => $pairs };
}

for my $xcov (sort { $a <=> $b } keys %cvalues) {
    say STDERR "=====> Working on $xcov...";
    $merge_thresh += 100;

    my $intseq = join_pairs($cvalues{$xcov}{reads}, $forward, $reverse);

    my $outdir   = File::Spec->catdir($opts{output_dir}, $base."_${xcov}X");
    my $logfile  = $base.'_transposome_'."${xcov}X".'.log';       
    my $clogfile = $base.'_transposome_clustering_'."${xcov}X".'.log';

    # for getting family counts
    my $sumfile  = File::Spec->catfile($outdir, $base.'_transposome_clustering_'."${xcov}X".'_annotations_summary.tsv'); 
    my $savfile  = File::Spec->catfile($sumdir, $base.'_transposome_clustering_'."${xcov}X".'_annotations_summary.tsv'); 

    my %run_opts = (
	coverage        => $xcov,
	sequence_file   => $intseq,
	sequence_format => $opts{sequence_format}, 
	threads         => $opts{threads},
	outdir          => $outdir,
	repeatdb        => $opts{repeat_database},
	cluster_size    => $opts{cluster_size},
	logfile         => $logfile,
	cluster_logfile => $clogfile);

    my $config = write_config(\%run_opts);
    system('transposome', '-c', $config) == 0 or die $!;

    my $reslog = File::Spec->catfile($outdir, $logfile);
    write_results($xcov, $cvalues{$xcov}{reads}, $reslog, $sumfile, $out);
    copy $sumfile, $savfile or die "[ERROR]: Copy failed: $!"; 
    remove_tree($outdir, { safe => 1 });
    unlink $intseq, $config;
}
say STDERR "=====> Done running Transposome. Results can be found in: $opts{output_dir}";

exit;
#
# Methods
#
sub join_pairs {
    my ($sample_size, $forward, $reverse) = @_;

    my $pair_size = sprintf("%.0f", $sample_size/2);
    my $intseq = File::Spec->catfile($opts{output_dir}, $base."_$sample_size".'_interl.fastq.gz');

    my $cmd = "pairfq joinpairs -f <(seqtk sample -s 11 $forward $pair_size) -r <(seqtk sample -s 11 $reverse $pair_size)";
    $cmd .= " -o $intseq -c gzip";

    system_bash($cmd);

    return $intseq;
}

sub write_results {
    my ($xcov, $sample_size, $logfile, $sumfile, $out) = @_;

    open my $in, '<', $logfile;

    my ($clstot, $uclstot, $repfrac, $singfrac, $reptot, $annotot, $time);
    while (my $line = <$in>) {
	chomp $line;
	if ($line =~ /Results - Total sequences clustered:\s+($RE{num}{real})/) {
	    $clstot = $1;
	}
	if ($line =~ /Results - Total sequences unclustered:\s+($RE{num}{real})/) {
	    $uclstot = $1;
	}
	if ($line =~ /Results - Repeat fraction from clusters:\s+($RE{num}{real})/) {
	    $repfrac = $1;
	}
	if ($line =~ /Results - Singleton repeat fraction:\s+($RE{num}{real})/) {
	    $singfrac = $1;
	}
	if ($line =~ /Results - Total repeat fraction \(theoretical\):\s+($RE{num}{real})/) {
	    $reptot = $1;
	} 
	if ($line =~ /Results - Total repeat fraction from annotations \(biological\):\s+($RE{num}{real})/) {
	    $annotot = $1;
	}
	if ($line =~ /======== Transposome completed at:.*Elapsed time: (.*). ========/) {
	    $time = $1;
	}
    }

    open my $sum, '<', $sumfile;
    
    my $famct = 0;
    while (my $line = <$sum>) {
	chomp $line;
	next if $line =~ /^ReadNum/;
	$famct++;
    }
    close $sum;

    say $out join "\t", "${xcov}X ($sample_size reads)", $clstot, $uclstot, $repfrac, $singfrac, $reptot, $annotot, $famct, $time;

    return;
}

sub write_config {
    my ($run_opts) = @_;

    my $config = 
"blast_input:
  - sequence_file:      $run_opts->{sequence_file}
  - sequence_format:    $run_opts->{sequence_format}
  - thread:             $run_opts->{threads}
  - output_directory:   $run_opts->{outdir}
clustering_options:
  - in_memory:          1
  - percent_identity:   90
  - fraction_coverage:  0.55
annotation_input:
  - repeat_database:    $run_opts->{repeatdb}
annotation_options:
  - cluster_size:       $run_opts->{cluster_size}
output:
  - run_log_file:       $run_opts->{logfile}
  - cluster_log_file:   $run_opts->{cluster_logfile}";
    
    my $config_file = "transposome_config_$run_opts->{coverage}X.yml";
    open my $out, '>', $config_file;
    say $out $config;
    close $out;

    return $config_file;
}

sub system_bash {
    my @args = ( 'bash', '-c', shift );
    system([0..5], @args);
}


sub usage {
    my $script = basename($0);
    print STDERR <<END

USAGE: $script -a FILEBASE -o transposome_results_out -t 24 -db sunflower_tephra_annotations_repbase-fmt.fasta -gs 3600

Required:
 -a|accession            :       An SRA accession (or simply the file basename). There should be two files found 
                                 in the working directory that end in '_1.fastq' and '_2.fastq'. In the example above,
                                 'FILEBASE_1.fastq.gz' and 'FILEBASE_2.fastq.gz' could be the actual file names. Note that
                                 the files may be compressed with gzip or bzip2, so the ending can be '.gz' or '.bzip2' but
				 uncompressed files with the actual '.fastq' ending works fine as well.
 -o|output_dir           :       A name for the output directory to hold results.
 -db|repeat_database     :       A sequence file of repeats to be used for annotation. This file must be formatted for
                                 use with Transposome, see the following link: 
				 https://github.com/sestaton/Transposome/wiki/Preparing-sequence-data-for-Tranposome
 -gs|genome_size         :       The 1C genome size expressed in Mbp (e.g., use 2500 for a 2.5 Gbp genome). 

Options:
 -sf|sequence_format     :       The input sequence format (Default: FASTQ).
 -pid|percent_identity   :       Percent identity between pairwise matches in all vs. all blast (Default: 90).
 -fcov|fraction_coverage :       The fraction coverage between pairwise matches in all vs. all blast (Default: 0.55).
 -cls|cluster_size       :       The minimum size of a cluster to be used for annotation (Default: 100).
 -t|threads              :       The number of parallel blast processes to run (Default: 1).
 -h|help                 :       Print a usage statement.

END
}

