#!/usr/bin/env perl

use 5.012;
use strict;
use warnings;
use autodie qw(open);
use File::Basename;
use Transposome::SeqFactory;
use Getopt::Long;

my $infile = shift or die $usage;
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

my $epithet = $genus.q{ }.$species;

open my $out, '>', $outfile;
my $seqio = Transposome::SeqFactory->new( file => $infile )->make_seqio_object;

while (my $seq = $seqio->next_seq) {
    # format ID
    my $id = $seq->get_id;
    $id =~ s/-|*/_/g;

    # format sequence
    my $nt = $seq->get_seq;
    $nt =~ s/.{60}\K/\n/g;

    if ($id =~ /^RLG/) {
        say $out join "\n", ">".$id."\t"."Gypsy"."\t".$epithet, $nt;
    }
    elsif ($id =~ /^RLC/) {
        say $out join "\n", ">".$id."\t"."Copia"."\t".$epithet, $nt;
    }
    elsif ($id =~ /^DHH/) {
        say $out join "\n", ">".$id."\t"."Helitron"."\t".$epithet, $nt;
    }
    elsif ($id =~ /^DTA/) {
        say $out join "\n", ">".$id."\t"."hAT"."\t".$epithet, $nt;
    }
    elsif ($id =~ /^DTC/) {
        say $out join "\n", ">".$id."\t"."CACTA"."\t".$epithet, $nt;
    }
    elsif ($id =~ /^DTH/) {
        say $out join "\n", ">".$id."\t"."PIF/Harbinger"."\t".$epithet, $nt;
    }
    elsif ($id =~ /^DTM/) {
        say $out join "\n", ">".$id."\t"."Mutator"."\t".$epithet, $nt;
    }
    elsif ($id =~ /^DTT/) {
        say $out join "\n", ">".$id."\t"."Tc1/Mariner"."\t".$epithet, $nt;
    }
    elsif ($id =~ /^PPP/) {
        say $out join "\n", ">".$id."\t"."Penelope"."\t".$epithet, $nt;
    }
    elsif ($id =~ /^RIL/) {
        say $out join "\n", ">".$id."\t"."L1"."\t".$epithet, $nt;
    }
    elsif ($id =~ /^RIT/) {
        say $out join "\n", ">".$id."\t"."RTE"."\t".$epithet, $nt;
    }
    elsif ($id =~ /^RLX/) {
        say $out join "\n", ">".$id."\t"."Unknown_LTR"."\t".$epithet, $nt;
    }
    elsif ($id =~ /^RST/) {
        say $out join "\n", ">".$id."\t"."tRNA"."\t".$epithet, $nt;
    }
    else {
        # should never get here, but the data may be malformed
        say STDERR "\n[ERROR]: ",$seq->get_id," does not seem to match known TE superfamilies";
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
-i|infile              :       A Fasta file of repeats to be formatted for Transposome.
-o|outfile             :       A name for the formatted database.
-g|genus               :       The source genus of the sequences.
-s|species             :       The source species of the sequences.

Options:
-h|help                :       Print a usage statement.

END
}
