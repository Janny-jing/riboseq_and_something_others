#!/usr/bin/env python
#-*- encoding=utf-8 -*-

import sys
from Bio import SeqIO

def format_refseq(fa_name, fmt_fa_name):

    dict_seqs = {}
    fa = SeqIO.parse(fa_name, 'fasta')
    for seq in fa:
        dict_seqs[seq.description] = str(seq.seq)
   
    with open(fmt_fa_name, 'w') as f:
        for dsp, seq in dict_seqs.items():
            seq = 'CCTAAGAGCAAGAAGAAGCCTGGN'+seq+'CCA'+'GGCTTCTTCTTGCTCTTAGGAAAAAAAAAA'
            f.write('>'+dsp+'\n')
            f.write(seq+'\n')
    
if __name__ == '__main__':
    format_refseq(sys.argv[1], sys.argv[2])