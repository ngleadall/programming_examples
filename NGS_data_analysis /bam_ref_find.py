# Script which looks through a .bamfile, filters reads for high quality
# writes these high quality reads to a _hq.bamfile and gives back the name of reference
# which most reads are aligned to.
# Nick Gleadall - nick.gleadall@googlemail.com - 27/07/2015

#print "Hello World!"

import pysam
import sys
import re

# Filters a bamfile for high quality reads, edit MIN_ALIGN_LENGTH and
# MIN_ALIGNMENT_SCORE as appropriate.
def make_HQ_bam(inbam, outbam):

    MIN_ALIGN_LENGTH    = 70
    MIN_ALIGNMENT_SCORE = 20


    infile  = pysam.AlignmentFile(inbam, "rb" )
    outfile = pysam.AlignmentFile(outbam, "wb", template=infile)

    for line in infile.fetch(until_eof=True):

        alignment_length = line.query_alignment_end - line.query_alignment_start + 1
        alignment_score  = line.mapping_quality

        if ( alignment_length < MIN_ALIGN_LENGTH or alignment_score < MIN_ALIGNMENT_SCORE ):
            line.is_qcfail = 1

        mate = ''

        if line.mate_is_unmapped == True:
            line.is_qcfail = 1

        #print line
        outfile.write( line );

    infile.close()
    outfile.close()


# Creates a dictionay of reference id's present in a bamfile, counts the number
# of reads which each reference and finally prints the true name of the reference.
def find_best_ref(in_hq_bam):

    subtype_scores = {}

    infile = pysam.AlignmentFile(in_hq_bam, "rb" )

    for line in infile.fetch(until_eof=True):

        if line.is_qcfail == 1:
            continue

        else:
            if line.is_qcfail == 0:
                if line.reference_id not in subtype_scores.keys():
                    subtype_scores[line.reference_id] = 1
            else:
                if line.reference_id in subtype_scores.keys():
                    subtype_scores[line.reference_id] += 1

    #print subtype_scores

    best_ref = max(subtype_scores, key=subtype_scores.get)
    best_ref_name = infile.getrname(best_ref)
    print "Name of reference with most HQ reads mapped in %s to is %s" % (in_hq_bam, best_ref_name)
    return best_ref_name

# Main loop function.
def ref_find(inbam):
    outbam = re.sub(r'(.*).bam', r'\1_HQ.bam', inbam)
    print "ref_find module has made a high quality bamfile called %s" % outbam
    make_HQ_bam(inbam, outbam)
    return find_best_ref(outbam)
