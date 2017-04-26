transposome-scripts
===================

Utility scripts using for working with Transposome

**DESCRIPTION**

This repository contains user-supplied scripts to help with many aspects of genome analysis using [Transposome](https://github.com/sestaton/Transposome). A lisitng of the included scripts is provided below.

<table>
<tr><th>Tool</th><th>Description</th></tr>
<tr><td><a href="https://github.com/sestaton/transposome-scripts/blob/master/sample_reads.pl">sample_reads.pl</a></td><td>Randomly sample a FASTA or FASTQ file for a desired number of sequences.</td></tr>
<tr><td><a href="https://github.com/sestaton/transposome-scripts/blob/master/get_seqs_from_list.pl">get_seqs_from_list.pl</a></td><td>Pull sequences from a file given a file of sequence IDs.</td></tr>
<tr><td><a href="https://github.com/sestaton/transposome-scripts/blob/master/sample_and_interleave_pairs.pl">sample_and_interleave_pairs.pl</a></td><td>Sample reads from both pairs of a paired-end read and interleave them for input to Transposome.</td></tr>
<tr><td><a href="https://github.com/sestaton/transposome-scripts/blob/master/permute_graph_edges.pl">permute_graph_edges.pl</a></td><td>Test a range of parameters for filtering pairwise matches to graph edges.</td></tr>
<tr><td><a href="https://github.com/sestaton/transposome-scripts/blob/master/permute_cluster_joining.pl">permute_cluster_joining.pl</a></td><td>Vary the threshold for joining clusters.</td></tr>
<tr><td><a href="https://github.com/sestaton/transposome-scripts/blob/master/format_database.pl">format_database.pl</a></td><td>Format a FASTA file of repeat sequences for use with Transposome.</td></tr>
<tr><td><a href="https://github.com/sestaton/transposome-scripts/blob/master/select_database_species.pl">select_database_species.pl</a></td><td>Select all sequences from a given species from a repeat database.</td></tr>
<tr><td><a href="https://github.com/sestaton/transposome-scripts/blob/master/full_analysis_with_varying_coverage.pl">full_analysis_with_varying_coverage.pl</a></td><td>Cluster and annotate sequences with varying genome coverage.</td></tr>
</table>

**USAGE**

All scripts will generate a usage statement if executed with no arguments, so to find out the proper usage of script simply type `perl` followed by the name of the script. For example,

    $ perl sample_reads.pl 
    USAGE: sample_reads.pl -i seqs.fasta -n 100000 
    
    Required:
     -i|infile            :       The sequence file to sample.
     -n|sample_size       :       The number of sequences to sample (integer).
                                  (NB: the value may be as '1000000' or with underscores
                                   as in '1_000_000').
    Options:
     -f|format            :       The input sequence format (Default: FASTA).
     -h|help              :       Print a usage statement.

**CONTRIBUTING**

If you write a script that uses Transposome please fork this repo, add your script and a description of it to the table above, and then send a pull request. I'll be happy to incorporate any script that looks useful. Please try to follow the coding style of the provided scripts as much as possible. At minimum, make sure the scripts handle options correctly and produce some kind of usage statement.

**CITATION**

If you use any of these tools in your work, please use the following citation: 

    Staton SE, and Burke JM. 2015. Transposome: A toolkit for annotation of transposable element families from unassembled sequence reads
        Bioinformatics, doi: 10.1093/bioinformatics/btv059

**LICENSE AND COPYRIGHT**

Copyright (C) 2014-2017 S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: http://www.opensource.org/licenses/mit-license.php