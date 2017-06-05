#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use Cwd     qw(getcwd);
use File::Basename;
use Transposome::Annotation;
use Transposome::SeqFactory;
use Getopt::Long;

my $infile; 
my $outfile; 
my $genus;
my $species;
my $quiet;
my $help;

GetOptions(
           'i|infile=s'  => \$infile,
           'o|outfile=s' => \$outfile,
           'g|genus=s'   => \$genus,
           's|species=s' => \$species,
           'q|quiet'     => \$quiet,
           'h|help'      => \$help,
          );

usage() and exit(0) if $help;

if (!$infile || !$outfile || !$genus || !$species) {
    say "\n[ERROR]: Command line not parsed correctly. Check input and refer to the usage. Exiting.\n";
    usage();
    exit(1);
}

my $outdir  = getcwd();
my $epithet = $genus.q{ }.$species;
my $report  = 'tmp_transposome_format_rep.txt';

# these class attributes are being faked because we just need one method
my $annotation = Transposome::Annotation->new(
      database => $infile, 
      dir      => $outdir,
      file     => $report,
      threads  => 1,
      cpus     => 1,
      verbose  => 0,
);

open my $out, '>', $outfile;
my $seqio = Transposome::SeqFactory->new( file => $infile )->make_seqio_object;

while (my $seq = $seqio->next_seq) {
    # format ID
    my $id = $seq->get_id;
    $id =~ s/\-|\*/_/g;

    # format sequence
    my $nt = $seq->get_seq;
    $nt =~ s/.{60}\K/\n/g;

    my $superfamily = $annotation->map_superfamily_name($id);

    if ($superfamily) { 
	say $out join "\n", ">".$id."\t".$superfamily."\t".$epithet, $nt;
    }
    else {
	if ($id =~ /#(\w+)\/?(\w+(\d+)?)?/) {
	    my $name = defined $2 ? $2 : $1;
	    $id =~ s/\s+.*//g; # this assumes IDs are unique up to the first space
	    say $out join "\n", ">".$id."\t".$name."\t".$epithet, $nt;
	    next;
	}

	say STDERR "\n[WARNING]: ",$seq->get_id," does not seem to match known TE superfamilies\n"
	    unless $quiet;
	$id =~ s/\s+.*//g; # this assumes IDs are unique up to the first space
	say $out join "\n", ">".$id."\tUnknown_repeat\t".$epithet, $nt;
    }
}

close $out;

exit;
## methods
sub usage {
    my $script = basename($0);
    print STDERR <<END
USAGE: $script [-i] [-o] [-g] [-s] [-h]

Required:
 -i|infile              :       A FASTA file of repeats to be formatted for Transposome.
 -o|outfile             :       A output filename for the formatted database.
 -g|genus               :       The source genus of the sequences.
 -s|species             :       The source species of the sequences.

Options:
 -q|quiet               :       Silence warnings when known superfamilies cannot be matched. Note
                                that some warnings will still come through with regards to the
                                ID format, but these are just warnings and can be ignored.
 -h|help                :       Print a usage statement.

END
}
