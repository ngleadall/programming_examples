# Script for HIV-1 virus subtyping using high quality mapped reads for reference selection.
#
# Nick gleadall - nick.gleadall@googlemail.com r - 29/07/15


# -------------------------------------------------------------------------------------
# Imports
# -------------------------------------------------------------------------------------

import sys
import pipeliners as pipe
import os
import re
import bam_ref_find as RF
import sample_fq as SF

# -------------------------------------------------------------------------------------
# Define Functions
# -------------------------------------------------------------------------------------
# Pipeline commands

def gzip( infile ):
	# gzips stuff
	return "gzip " + infile

def cutadapt( infile, outfile):
	# Removes sequencing adapters
    return "/software/packages/cutadapt-1.1/bin/cutadapt -b TGTAGAACCATGTCGTCAGTGT -b AGACCAAGTCTCTGCTACCGT "  + infile  + " | gzip -c  > " + outfile

def smalt(refname, infile1, infile2, outfile):
    # Align the reads to the reference (run smalt-0.7.6 to see all options)
    return "/software/bin/smalt_0.7.6 map -f samsoft /refs/HIV/HIV_subtype_refs/"+ refname + " " + infile1 + " " + infile2 + " > " + outfile

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

# -------------------------------------------------------------------------------------
# Data I/O
# -------------------------------------------------------------------------------------

if (len(sys.argv) == 1):
    print "Needs sample name for analysis"
    exit()

sample_name = sys.argv[1]
# Get sample name and confirm to user which one is being worked on.
print "Working on sample \" %s \"" % sample_name
sample_fq1 = sample_name + ".1.fq.gz"
sample_fq2 = sample_name + ".2.fq.gz"


# Cut sequencing adapters from fq1 and fq2.
pipe.system_call("Removing sequencing adapters from fastq 1", cutadapt(sample_fq1, sample_name + "_cut.1.fq.gz"))
pipe.system_call("Removing sequencing adapters from fastq 2", cutadapt(sample_fq2, sample_name + "_cut.2.fq.gz"))

# Sample 10,000 random reads from fq1 and fq 2, write these reads to two files
# Make two variables for easy naming of the fastq's containing 10k subset of reads
fq_subset_1 = sample_name + "_cut_subset.1.fq"
fq_subset_2 = sample_name + "_cut_subset.2.fq"

print "Sampling random 10,000 reads from %s fastq's and printing them into \" %s_subset.1/2.fq.gz \"" % (sample_name, sample_name)
SF.res_sample(sample_name+"_cut", 10000)
# gzip the subsets - NEEDS ADDING TO THE SF script protperly!!!
pipe.system_call("zipping _subset fastq file", gzip(fq_subset_1))
pipe.system_call("zipping _subset fastq file", gzip(fq_subset_2))


# Rename these variable to account for their new .gz status!
fq_subset_1 = sample_name + "_cut_subset.1.fq.gz"
fq_subset_2 = sample_name + "_cut_subset.2.fq.gz"

# Map subset fastq to all HIV references using SMALT
# Create new variable for the sam and bam subset files
subset_samfile = sample_name + "_subset.sam"
subset_bamfile = sample_name + "_subset.bam"

pipe.system_call("Mapping subset of reads from sample \" %s \" to Los Alamos subtype reference panel", smalt("HIV_subtypes_all_s1k6", fq_subset_1, fq_subset_2, subset_samfile))
#Convert subset sam to subset bam
pipe.system_call('Converting samfile to bamfile', sam_to_bam(subset_samfile, subset_bamfile))

# Quality filter the bamfile and return the best reference name.
print "Making a higher quality bamfile and finding the best subtype reference for this sample"
best_ref_name = RF.ref_find(subset_bamfile)
print "Best reference for this sample is %s " % best_ref_name

# Map the origional fastq's to the best reference using smalt
# Define a variable for the sam and bam full sample files
full_sam = sample_name + ".sam"
full_bam = sample_name + ".bam"

pipe.system_call("Mapping all reads from " + sample_name + " to " + best_ref_name, smalt(best_ref_name + "_s1k6", sample_fq1, sample_fq2, full_sam ))

#Convert full sam to full bam
pipe.system_call('Converting samfile to bamfile', sam_to_bam(full_sam, full_bam))

print "%s has been made using %s as a reference. Further anlysis to be conducted at user discretion." %s full_bam
