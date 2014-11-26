transposome-scripts
===================

Utility scripts using for working with Transposome

**DESCRIPTION**
This repository contains user-supplied scripts to help with many aspects of genome analysis using Transposome. A lisitng of the included scripts is provided below.

<table>
<tr><th>Tool</th><th>Description</th></tr>
<tr><th><a href="https://github.com/sestaton/transposome-scripts/blob/master/sample_reads.pl">sample_reads.pl</a></th><td>Randomly sample a FASTA or FASTQ file for a desired number of sequences.</td></tr>
<tr><th><a href="https://github.com/sestaton/transposome-scripts/blob/master/get_seqs_from_list.pl">get_seqs_from_list.pl</a></th><td>Pull sequences from a file given a file of sequence IDs.</td></tr>
<tr><th><a href="https://github.com/sestaton/transposome-scripts/blob/master/sample_and_interleave_pairs.pl">sample_and_interleave_pairs.pl</a></th><td>Sample reads from both pairs of a paired-end read and interleave them for input to Transposome.</td></tr>
</table>

**CONTRIBUTING**

If you write a script that uses Transposome please fork this repo, add your script and a description of it to the table above, and then send a pull request. I'll be happy to incorporate any script that looks useful. Please try to follow the coding style of the provided scripts as much as possible. At minimum, make sure the scripts handle options correctly and produce some kind of usage statement.

**CITATION**

If you use any of these tools in your work, please use the following citation: 
[![DOI](https://zenodo.org/badge/5565/sestaton/Transposome.png)](http://dx.doi.org/10.5281/zenodo.11303)

    Staton SE, and Burke JM. 2014. Transposome: Annotation of transposable element families
         from unassembled sequence reads. doi:10.5281/zenodo.11303

**LICENSE AND COPYRIGHT**

Copyright (C) 2014 S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: http://www.opensource.org/licenses/mit-license.php