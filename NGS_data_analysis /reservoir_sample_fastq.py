# Streams a pair of fastq files and reservoir samples 10k reads, writing them to a new fastq for QC and
# analysis purposes.
# Nick Gleadall & Kim Brugger - 14/04/15


#----------------------------------------------------------------------------------
# Importing Modules
#----------------------------------------------------------------------------------

import subprocess
import sys
import random

#----------------------------------------------------------------------------------
# Define Functions
#----------------------------------------------------------------------------------

# Reads in a fastq entry from a filehandle
# Every entry consists of 4 lines.
# if at the end of the file return an empty array
def read_fq_entry( process ):

    while ( True ):
    	name   = process.stdout.readline()
    	seq    = process.stdout.readline()
    	strand = process.stdout.readline()
    	qual   = process.stdout.readline()

    	if ( not name or name == '') :
        	return []
        	break

    	return [ name, seq, strand, qual ]


# System calls zcat on two .gz files and streams the output, selects n reads randomly and writes these to a new file.
def stream_fq(sample_fq1, sample_fq2, out_fq1, out_fq2, Nr_of_reads):

	# Variable for storing random reads from each fastq in the pair.
	random_reads_1 = []
	random_reads_2 = []

	# Read index counter - sampling counter
	read_index = 0

	proc1 = subprocess.Popen("zcat < " + sample_fq1, stdout=subprocess.PIPE, shell=True)
	proc2 = subprocess.Popen("zcat < " + sample_fq2, stdout=subprocess.PIPE, shell=True)

	while ( True ):

		fq_1_entry = read_fq_entry(proc1)
		fq_2_entry = read_fq_entry(proc2)

		#print fq_1_entry
		#print fq_2_entry
		read_index += 1

		if (read_index <= Nr_of_reads):
			random_reads_1.append( fq_1_entry )
			random_reads_2.append( fq_2_entry )

		else:
			random_position = random.randrange(0, read_index, 1)

			if ( random_position < Nr_of_reads -1):
				random_reads_1[ random_position ] = fq_1_entry
				random_reads_2[ random_position ] = fq_2_entry

		if fq_1_entry == [''] or fq_1_entry == []:
			break

	#print "printing reads \n"
	#print random_reads_1
	#print "printing readset 2 \n"
	#print random_reads_2

	for read in random_reads_1:
		out_fq1.write(''.join(read))
	for read2 in random_reads_2:
		out_fq2.write(''.join(read2))



#----------------------------------------------------------------------------------
# Main Loop
#----------------------------------------------------------------------------------

def res_sample(insample_id, read_number):

    # Generate fastq 1/2 names from insample_id
    sample_fq1 = insample_id + ".1.fq.gz"
    sample_fq2 = insample_id + ".2.fq.gz"

    # Open two subset files for lines to go into
    outfq_name_1 = insample_id + "_subset.1.fq"
    outfq_name_2 = insample_id + "_subset.2.fq"
    out_fq1 = open(insample_id + "_subset.1.fq", 'w')
    out_fq2 = open(insample_id + "_subset.2.fq", 'w')

    stream_fq(sample_fq1, sample_fq2, out_fq1, out_fq2, read_number)

    out_fq1.close()
    out_fq2.close()


if len(sys.argv) == 3:
    res_sample( str(sys.argv[1]) , str(sys.argv[2]) )
else:
    print "Please specify a sample name and desired number of reads as arguments"
