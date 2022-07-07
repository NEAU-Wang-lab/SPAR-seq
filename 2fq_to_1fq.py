import sys, regex, pysam, mappy
from itertools import islice

def revComp(s):
        revCDict = {'G':'C','C':'G','T':'A','A':'T','N':'N','-':'-'}
        return ''.join([revCDict[x] for x in s])[::-1]

fq_2_dict = {}
infile_2 = sys.argv[2]
infile_2 = mappy.fastx_read(infile_2, read_comment=True)
for read in infile_2:
	name = read[0][:-2]
	seq  = read[1]
	seq  = revComp(seq)
	qual = read[2]
	qual = qual[::-1]
	fq_2_list[name] = (seq, qual)

fq_merge = sys.argv[3] + '_merge.fq'
fq_merge = open(fq_merge, 'w')
infile_1 = sys.argv[1]
infile_1 = mappy.fastx_read(infile_1, read_comment=True)
for read in infile_1:
	name   = read[0][:-2]
	seq_1  = read[1]
	qual_1 = read[2]
	seq_2  = fq_2_list[name][0]
	qual_2 = fq_2_list[name][1]
	seq    = seq_1 + seq_2
	qual   = qual_1 + qual_2

	fq_merge.write(name + '\n')
	fq_merge.write(seq + '\n')
	fq_merge.write('+\n')
	fq_merge.write(qual + '\n')

fq_merge.close()