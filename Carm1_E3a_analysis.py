import sys, regex, pysam, mappy
from itertools import islice

def revComp(s):
        revCDict = {'G':'C','C':'G','T':'A','A':'T','N':'N','-':'-'}
        return ''.join([revCDict[x] for x in s])[::-1]

infile = sys.argv[1]
infile = mappy.fastx_read(infile, read_comment=True)
for read in infile:
	name = read[0]
	seq  = read[1]
	qual = read[2]

	pattern = "(GCTCCAAGTCACACAGCTCTTATGC){e<=3}"
	matches = regex.finditer(pattern, seq, overlapped = False)
	for it in matches:
		bar_name = 'E3a'


	seq = revComp(seq)
	qual = qual[::-1]
	pattern = "(GCTCCAAGTCACACAGCTCTTATGC){e<=3}"
	matches = regex.finditer(pattern, seq, overlapped = False)
	for it in matches:
		bar_name = 'E3a'
 
	out = [name, seq, qual, bar_name]
	out = '\t'.join([str(x) for x in out])
	print(out)

