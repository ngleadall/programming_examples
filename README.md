# Programming Examples

This repo contains examples of python and r code created by Nicholas Gleadall.
Below is a brief explanation of each tool.  

NGS_data_analysis
=================
<body>
<b>bam_ref_find.py</b><BR>

Usage = $ python bam_ref_find.py "bamfile"

Script which streams a .bamfile, 'soft clips' reads below the
desired quality and then writes the reamining reads to a new HQ bamfile.

It then streams the new HQ bamfile, bulding a dictionary of reference ID' and
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
