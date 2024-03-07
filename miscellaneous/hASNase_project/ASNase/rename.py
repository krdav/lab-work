from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import sys

filename = 'ASNase_human_like_dna_seq_id50_MAFFTalign.fasta'

name2number = dict()
for i, record in enumerate(SeqIO.parse(open(filename, 'r'), 'fasta')):
    name2number[record.id] = i
    print('>{}\n{}'.format(i, str(record.seq)))
    # print('{}, {}'.format(record.id, i))



