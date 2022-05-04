from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import sys
import distance
min_similarity = 0.5

gb_file = 'ASNase_human_like_nuc_seq.gb'
sources2seqs = dict()
sources = dict()
seq_set = set()
comp_seqs = list()
for gb_record in SeqIO.parse(open(gb_file, 'r'), 'genbank'):
    # print(gb_record.name)
    # print(gb_record.annotations['source'])
    source = '_'.join(gb_record.annotations['source'].split())
    dup = 1
    while True:
        source_numb = source + '_' + str(dup)
        if source_numb in sources2seqs:
            dup += 1
        else:
            source = source_numb
            break
    # print(source)
    for feature in gb_record.features:
        if feature.type == 'CDS':
            # print(feature)
            cds_feature = feature
            break
    start = cds_feature.location.nofuzzy_start
    end = cds_feature.location.nofuzzy_end - 3
    # print(start, end)
    # print(cds_feature.location.strand)

    if cds_feature.location.strand == 1:
        seq = str(gb_record.seq[start:end])
    else:
        seq = str(gb_record.seq[start:end].reverse_complement())
    if seq not in seq_set and seq[0:3].upper() == 'ATG' and '*' not in str(Seq(seq, generic_dna).translate()):
        sources2seqs[source] = seq
        seq_set.add(seq)
        sources[source] = len(sources) + 1
        if 'human' in source.lower() or 'guinea_pig' in source.lower() or 'zebrafish' in source.lower() or 'fruit_fly' in source.lower():
            comp_seqs.append(seq)

    # print(str(gb_record.seq))
    # print(seq)
    # sys.exit()

# print(comp_seqs)
for name, seq in sources2seqs.items():
    # print('>{}\n{}'.format(sources[name], seq))
    # print('>{}\n{}'.format(name, seq))
    # print('>{}\n{}'.format(name, Seq(seq, generic_dna).translate()))
    dist_cut = False
    for comp_seq in comp_seqs:
        comp_seq_t = str(Seq(comp_seq, generic_dna).translate())
        seq_t = str(Seq(seq, generic_dna).translate())
        # print(comp_seq_t, seq_t)
        edits = distance.levenshtein(comp_seq_t, seq_t)
        # print(edits)
        seq_len = max((len(comp_seq_t), len(seq_t)))
        similarity = (seq_len - float(edits)) / seq_len
        # print(similarity)
        if similarity > min_similarity:
            dist_cut = True
            break
    if dist_cut:
        print('>{}\n{}'.format(name, Seq(seq, generic_dna).translate()))





