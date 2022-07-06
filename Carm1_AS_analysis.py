import sys, regex, pysam, mappy
from itertools import islice

def revComp(s):
        revCDict = {'G':'C','C':'G','T':'A','A':'T','N':'N','-':'-'}
        return ''.join([revCDict[x] for x in s])[::-1]

barcode_dict = {}
barcodefile = sys.argv[2]
barcodefile = pysam.FastxFile(barcodefile)
for read in islice(barcodefile, None):
        barcode_dict[read.name] = read.sequence

AS_type = []
infile = sys.argv[1]
infile = mappy.fastx_read(infile, read_comment=True)
for read in infile:
	name = read[0]
	seq  = read[1]
	qual = read[2]
	bar_name = 'other'
	for barcode_name, barcode_seq in barcode_dict.items():
		pattern = "("+ barcode_seq + "){e<=3}"
		matches = regex.finditer(pattern, seq, overlapped = False)
		for it in matches:
			bar_name = barcode_name

	seq = revComp(seq)
	qual = qual[::-1]
	for barcode_name, barcode_seq in barcode_dict.items():
		pattern = "("+ barcode_seq + "){e<=3}"
		matches = regex.finditer(pattern, seq, overlapped = False)
		for it in matches:
			bar_name = barcode_name
 
	out = [name, seq, qual, bar_name]
	out = '\t'.join([str(x) for x in out])
	print(out)
