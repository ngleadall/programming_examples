# Programming Examples

This repo contains examples of python and r code created by Nicholas Gleadall.
nicholas.gleadall@googlemail.com
Below is a brief explanation of each tool.  

NOTE: Most of these are modified examples from a larger analysis framework. There
may be dependencies which are not installed on your system.


NGS_data_analysis
=================
<body>
<b>bam_ref_find.py</b><BR>

Usage = $ python bam_ref_find.py "bamfile"

Script which streams a .bamfile, 'soft clips' reads below the
desired quality and then writes the remaining reads to a new HQ bamfile.

It then streams the new HQ bamfile, building a dictionary of reference ID' and
counts how many reads are mapped to each one. It then returns the translated name
of most "mapped to" reference strain.

This tool can be easily adapted for many uses e.g. selection of most appropriate alternate sequences
in the latest human genome build HG38/GRCh38.

<b> reservoir_sample_fastq.py </b><BR>

Usage = $ python reservoir_sample_fastq.py "sample name" "number of reads desired"

Streams reads from a pair of fastq files (produced by Illumina sequencing), uses a reservoir sampling technique
to write a set number of random read pairs to two new fastq files.

This tool is for subsetting random reads from samples as required for bespoke analysis.

</body>

R_data_plotting
===============
<body>
<b>genome_coverage.R</b><BR>
Usage = $ Rscript genome_coverage.R "depth_file"

Simple plotting of data contained in .depth files generated by samtools.

"HIV_full_genome.png" is an example of the ouput, this shows number of reads (x) at a given genomic position (y).
y values can be edited as desired.

This was used as a quick and easy way of visualising pre sequencing, PCR amplicon coverage performance.

<b>Boxplot_freqs.R</b>

Example of an R analysis on more complex data. The aim was to determine the
reliability of mutation prevalences obtained from a HIV-1 sample by NGS.

The "Boxplot_example.png" is an example output of this analysis, it shows the mutations observed in 10 repeats
of a single sample and the prevelance at which each mutation was observed at.
</body>

analysis_pipelines
==================
<body>
Here are two examples of two bespoke automated analysis pipelines for assembly of next generation sequencing data obtained from the Illumina platform.

<b>HIV_subtyping_pipeline.py</b><BR>

Usage = $ python HIV_subtyping_pipeline.py "sample name"

The aim of this pipeline is to identify the subtype of a given HIV-1 sample.

The pipeline links together previous tools such as "reservoir_sample_fastq" and "bam_ref_find"
in order (1) quickly identify which reference strain a subset of the reads are mapping to then (2)
report this information and then (3) remap all original reads to this reference strain alone.

<b>read_unclip_analysis.py</b><BR>

Usage = python read_unclip_analysis.py "sample name"

This pipeline produces two sequence alignments from raw read data (fastq's), one which is a "standard" alignment
and another which all "soft-clipped/hidden" reads are "unclipped/shown".

This can be helpful when designing alignment algorithms to visualise the level of quality filtering being applied by
an alignment tool (here SMALT is used).

<b> pipeliners.py </b><BR>

This script was included as it is a commonly used tool. Its function is to be called within other scripts and
provides a function for creating subprocesses.

This script has the added advantage of manually setting verbose level (for reporting of errors during de-bugging)
new tools.

<hr>

I would like to thank <b>Dr. Kim Brugger</b> for creative input, guidance and supervision during the creation of these tools.
</body>
