# Programming Examples

This repo contains examples of python and r code created by Nicholas Gleadall.
Below is a brief explanation of each tool.  

NGS_data_analysis
=================
<body>
<b>bam_ref_find.py</b><BR>
Script which streams a .bamfile, 'soft clips' reads below the
desired quality and then writes the reamining reads to a new HQ bamfile.

It then streams the new HQ bamfile, bulding a dictionary of reference ID' and
counts how many reads are mapped to each one. It then returns the translated name
of most "mapped to" reference strain.

The original inteniton of this script was to subtype HIV-1 strains using NGS data,
however it can be adapted for other uses (e.g. selection of most appropriate alt sequences
  in the latest human genome build HG38/GRCh38)
</body>
