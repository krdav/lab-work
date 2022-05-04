from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import sys

aa_file = 'ASNase_human_like_aa_seq_id50_MAFFTalign.fasta'
dna_file = 'ASNase_human_like_nuc_seq.fasta'

dna_seqs = dict()
for record in SeqIO.parse(open(dna_file, 'r'), 'fasta'):
    dna_seqs[record.id] = str(record.seq)

for record in SeqIO.parse(open(aa_file, 'r'), 'fasta'):
    dna_seq = dna_seqs[record.id]
    dna_align = ''
    i = 0
    for aa in str(record.seq):
        if aa == '-':
            dna_align += '---'
        else:
            dna_align += dna_seq[i:(i+3)]
            i += 3

    assert(str(Seq(dna_align.replace('-', ''), generic_dna).translate()) == str(record.seq).replace('-', ''))
    # print(str(Seq(dna_align.replace('-', ''), generic_dna).translate()))
    # print(str(record.seq).replace('-', ''))
    print('>{}\n{}'.format(record.id, dna_align))



