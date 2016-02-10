# Pipeline for un-clipping hidden reads in a bamfile and writing them to new output
# for quality control analysis.
# 23rd April 2015
# Nicholas S. Gleadall - ng384@cam.ac.uk
# Kim Brugger - kim@brugger.dk

# print "Hello World!"

#------------------------------------------------------------
# Importing Modules

import sys
import subprocess
import os
import re
sys.path.append("/Users/Nick/github/pipeliners/modules")

import pipeliners

#------------------------------------------------------------
# Setting globals / functions

def cutadapt( infile, outfile):
	# Removes sequencing adapters
    return "/software/packages/cutadapt-1.1/bin/cutadapt -b TGTAGAACCATGTCGTCAGTGT -b AGACCAAGTCTCTGCTACCGT "  + infile  + " | gzip -c  > " + outfile

def smalt( infile1, infile2, outfile):
    # Align the reads to the reference (run smalt-0.7.6 to see all options)
    cmd  = "cd " + analysis_dir + ";"
    cmd += "/software/bin/smalt_0.7.6 map -f samsoft /refs/HIV/K03455_s1k6 " + infile1 + " " + infile2 + " > " + outfile
    return cmd

def sam_to_bam(infile, outfile):
	# Converts samfile to bamfile
	return "/software/bin//samtools view -Sb " + infile + " -o " + outfile

def sort_bamfile(infile, outfile):
	# Sorts bamfile
	return "/software/bin/samtools sort " + infile + " " + outfile

def deduplicate_bamfile(infile, outfile, outcsv):
	# Marks and soft clips duplicate reads
	return "/software/bin/picard -T MarkDuplicates I= " + infile + " O= " + outfile + " AS=true M= " + outcsv

def index_bamfile(infile):
	# Index's the bam file so it can be viewed in IGV later on
	return "/software/bin/samtools index " + infile

def HIV_alignment_fix(infile, outfile):
	# Fix's HIV specific alignment issues
	return "/software/packages/HIV-pipeline/scripts/bam_fix_indels.pl " + infile + " " + outfile

def unclip_bamfile(infile):
	# Run's kims unclipping script
	return "/software/bin/scripts/bam_unclip_bases.pl " + infile


#------------------------------------------------------------
# Define fastq's and get sample names.
# Gets location of fastq
if not (len(sys.argv) == 2):
    print "Need's sample name as input"
    exit()

sample = sys.argv[2]
sample_name = re.sub(r"(.*).1.fq.gz", r"\1", sample)
sample_name_uc = sample_name + "_unclipped"
# Setup new analysis directory



# Preparation of bamfile + unclipped bamfile
pipeliners.system_call('Cut Adapters', cutadapt(sample_name + ".1.fq.gz", analysis_dir + "/" + sample_name + "_ar.1.fq.gz"))
pipeliners.system_call('Cut Adapters', cutadapt(sample_name + ".2.fq.gz", analysis_dir + "/" + sample_name + "_ar.2.fq.gz"))

pipeliners.system_call('SMALT alignment', smalt(sample_name + "_ar.1.fq.gz", sample_name + "_ar.2.fq.gz", sample_name + ".sam"))
pipeliners.system_call('Convert samfile to bamfile', sam_to_bam(sample_name + ".sam", sample_name + ".bam"))

pipeliners.system_call('Unclipping bamfile', unclip_bamfile(sample_name + ".bam"))


# Bamfile processing pipeline
pipeliners.system_call('Sort the bamfile', sort_bamfile(sample_name+ ".bam", sample_name + "_sorted"))
pipeliners.system_call('Deduplicate bamfile', deduplicate_bamfile(sample_name + "_sorted.bam", sample_name + "_rmdups.bam", sample_name + "_rmdup.csv"))
pipeliners.system_call('Index bamfile', index_bamfile(sample_name + "_rmdups.bam"))
pipeliners.system_call('Fix HIV alignment', HIV_alignment_fix(sample_name + "_rmdups.bam", sample_name + "_fixed.bam"))
pipeliners.system_call('Index bamfile', index_bamfile(sample_name + "_fixed.bam"))

# Bamfile_unclipped processing pipeline
pipeliners.system_call('Sort the unclipped bamfile', sort_bamfile(sample_name_uc+ ".bam", sample_name_uc + "_sorted"))
pipeliners.system_call('Deduplicate unclipped bamfile', deduplicate_bamfile(sample_name_uc + "_sorted.bam", sample_name_uc + "_rmdups.bam", sample_name_uc + "_rmdup.csv"))
pipeliners.system_call('Index unclipped bamfile', index_bamfile(sample_name_uc + "_rmdups.bam"))
pipeliners.system_call('Fix HIV alignment (unclipped)', HIV_alignment_fix(sample_name_uc + "_rmdups.bam", sample_name_ + "_fixed.bam"))
pipeliners.system_call('Index unclipped fixed bamfile', index_bamfile(sample_name_uc + "_fixed.bam"))


exit()
