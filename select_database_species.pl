#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use File::Basename;
use Getopt::Long;

my $infile;
my $outfile; 
my $genus;
my $species;
my $help;

GetOptions(
           'i|infile=s'  => \$infile,
           'o|outfile=s' => \$outfile,
           'g|genus=s'   => \$genus,
           's|species=s' => \$species,
           'h|help'      => \$help,
          );

usage() and exit(0) if $help;

if (!$infile || !$outfile || !$genus || !$species) {
    say "\n[ERROR]: Command line not parsed correctly. Check input and refer to the usage. Exiting.\n";
    usage();
    exit(1);
}

my ($seq, $seqid, @seqparts);

open my $in, '<', $infile;
open my $out, '>', $outfile;

{
    local $/ = '>';

    while (my $line = <$in>) {
	chomp $line;
	($seqid, @seqparts) = split /\n/, $line;
	$seq = join '', @seqparts;
	next unless defined $seqid && defined $seq;

	# format ID
	my ($id, $type, @epithet)  = split /\s+/, $seqid;
	next unless defined $type && @epithet;

	my ($gen, $sp) = @epithet;
	next unless defined $gen && defined $sp;

	if (lc($gen) eq lc($genus) && lc($sp) eq lc($species)) {
	    # format sequence
	    $seq =~ s/.{60}\K/\n/g;
	    
	    say $out join "\n", ">".$seqid, $seq;
	}
    }
}
close $in;
close $out;

exit;
## methods
sub usage {
    my $script = basename($0);
    print STDERR <<END

USAGE: $script [-i] [-o] [-g] [-s] [-h]

Required:
 -i|infile              :       A Fasta file of repeats to be formatted for Transposome.
 -o|outfile             :       A name for the selected sequences.
 -g|genus               :       The source genus of the sequences.
 -s|species             :       The source species of the sequences.

Options:
 -h|help                :       Print a usage statement.

END
}
