#!/usr/bin/env python
import re
import sys
import regex
import pysam
import parasail
from itertools import islice

def revComp(s):
        revCDict = {'G':'C','C':'G','T':'A','A':'T','N':'N','-':'-'}
        return ''.join([revCDict[x] for x in s])[::-1]

def splitCCSs(read_name, read_pass, barcode_dict, seq, qual):
        carm = []
        all_bc_loci = [(0, 0)]

        for barcode_name, barcode_seq in barcode_dict.items():
                bc_loci = []
                pattern = "(ATGTCCACCTCGCAC){e<=1}("+ barcode_seq + ")"
                matches = regex.finditer(pattern, seq, overlapped = False)
                for it in matches:
                        bc_loci.append((barcode_name, it.end()))
                all_bc_loci.extend(bc_loci)
        sort_loci = sorted(all_bc_loci, key=lambda all_bc_loci: all_bc_loci[1])

        for i in range(len(sort_loci)-1):
                barcode_name = sort_loci[i+1][0]
                barcode_seq = barcode_dict[barcode_name]
                target_seq = seq[sort_loci[i][1]:sort_loci[i+1][1]]
                target_qual = qual[sort_loci[i][1]:sort_loci[i+1][1]]
                read = (read_name, read_pass, barcode_name, target_seq, target_qual)
                carm.append(read)
        return carm

infile = sys.argv[1]
passfile = sys.argv[2]
barcodefile = sys.argv[3]

if len(sys.argv) != 4:
        msg = ['python', sys.argv[0], 'infile', 'barcodefile','passfile', '1>out.txt', '2>err.txt']
        msg = ' '.join([str(x) for x in msg])
        print(msg)
        print('\tinfile: Input raw CCS reads in fastq format')
        print('\tpassfile: CCS read pass file')
        print('\tbarcodefile: Barcode file in fasta format')
        print('\tPlease contact WJQ if you have any questions')
        quit()

tag_dict = {}
passfile = open(passfile,'r')
for line in passfile:
        (read_name, read_pass) = line.strip('\n').split('\t')
        pass_dict[read_name] = read_pass

carms = []
ccs_count = 0
infile = pysam.FastxFile(infile)
for ccs in islice(infile, None):
        ccs_count = ccs_count + 1
        ccs_name = ccs.name
        ccs_pass = pass_dict[ccs_name]

        seq = ccs.sequence
        qual = ccs.quality
        carm = splitCCSs(ccs_name, ccs_pass, barcode_dict, seq, qual)

        strand = '-'
        seq = revComp(ccs.sequence)
        qual = ccs.quality[::-1]
        carm.extend(splitCCSs(ccs_name, ccs_pass, barcode_dict, seq, qual))
        carms.extend(carm)

for read in carms:
        end = start = ''
        (read_name, read_pass, carm_barcode, seq, qual) = read
        primer_r = '(TCCCCACCATTGGTG){e<=1}'
        matches = regex.finditer(primer_r, seq, overlapped = False)
        for it in matches:
                end = it.end()
        if end != '':
                seq = seq[:(end + 15)]
                qual = qual[:(end + 15)]

                primer_f = '(ACAGCGTCCTCATCCAGTTTGC){e<=2}'
                matches = regex.finditer(primer_f, seq, overlapped = False)
                for it in matches:
                        start = it.start()
                if start != '':
                        seq = seq[start:]
                        qual = qual[start:]
                        out = [read_name, read_pass, carm_barcode, 'expected', seq, qual]
                        out = '\t'.join([str(x) for x in out])
                        print(out)
                else:
                        err= [read_name, read_pass, carm_barcode, 'Wrong5', seq, qual]
                        err = '\t'.join([str(x) for x in err])
                        sys.stderr.write(err + '\n')

        else:
                err = [read_name, read_pass, carm_barcode, 'Wrong3', seq, qual]
                err = '\t'.join([str(x) for x in err])
                sys.stderr.write(err + '\n')

logfile = sys.argv[1] + '.log.txt'
log_file = open(logfile, 'w')
log_file.write(str(ccs_count) + '\n')
log_file.close()
