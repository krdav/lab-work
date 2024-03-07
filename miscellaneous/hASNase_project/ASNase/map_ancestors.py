from __future__ import division, print_function
import sys
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from ete3 import Tree, TreeNode, NodeStyle, TreeStyle, TextFace, add_face_to_node, CircleFace, faces, AttrFace
from Bio import AlignIO

def map_asr_to_tree(asr_seq, leaf_seq, tree, id_dict):
    '''Takes a IQ-TREE asr states and returns the matching ete3 tree node.'''

    DNA_order = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    # Parse input sequences:
    leafs = list(SeqIO.parse(leaf_seq, "phylip"))
    leafs = {r.id: str(r.seq).upper() for r in leafs}

    # Parse the ASR states from IQ-TREE:
    flag = False
    node_seqs = dict()
    with open(asr_seq) as fh:
        for l in fh:
            if not l.startswith('#') and not flag:
                assert(l.startswith('Node'))
                flag = True
            elif flag:
                state = l.strip().split()
                if state[0] not in node_seqs:
                    node_seqs[state[0]] = list()
                max_prop = max(map(float, state[3:7]))
                MAP_base = set([DNA_order[i] for i, p in enumerate(map(float, state[3:7])) if p == max_prop])
                if state[2] in MAP_base:
                    node_seqs[state[0]].append(state[2])
                elif state[2] == '-':
                    node_seqs[state[0]].append(MAP_base.pop())
                else:
                    raise Exception('Sequence reconstruction by IQTREE is inconsistent: {}'.format(l))
                assert(len(node_seqs[state[0]]) == int(state[1]))
    node_seqs = {k: ''.join(v) for k, v in node_seqs.items()}

    s1 = {len(s) for s in node_seqs.values()}
    s2 = {len(s) for s in leafs.values()}
    assert(s1 == s2)

    # Add info to tree:
    for node in tree.traverse():
        if node.name in leafs:
            node_seq = leafs[node.name]
        elif node.name in node_seqs:
            node_seq = node_seqs[node.name]
        else:
            raise Exception('Could not find node name in ancestral or input sequences:', node.name)

        node.add_feature('sequence', node_seq)
        if node.name in id_dict:
            node.name = id_dict[node.name]

    return tree




tree_file = 'ASNase_tree/ASNase_human_like_aa_seq_id50_MAFFTalign_rename.phylip.treefile'
leaf_seq = 'ASNase_tree/ASNase_human_like_aa_seq_id50_MAFFTalign_rename.phylip'
asr_seq = 'ASNase_tree/ASNase_human_like_aa_seq_id50_MAFFTalign_rename.phylip.state'
ids = 'ASNase_human_like_aa_seq_id50_MAFFTalign_rename.ids'

tree = Tree(tree_file, format=1)
id_dict = dict()
with open(ids) as fh:
    for l in fh.readlines():
        source, ID = l.split()
        source = source.replace('(', '@')
        source = source.replace(')', '@')
        id_dict[ID] = source

tree = map_asr_to_tree(asr_seq, leaf_seq, tree, id_dict)
print(tree)


huID = '@human@_1'
gpID = '@domestic_guinea_pig@_1'
for node in tree.traverse():
    if huID in node.name:
        huNODE = node
    elif gpID in node.name:
        gpNODE = node

common_ancestor = tree.get_common_ancestor([huNODE, gpNODE])
print(common_ancestor)

lineage = list()
while huNODE.up is not common_ancestor:
    # print(huNODE.name, huNODE.sequence)
    lineage.append([huNODE.name, huNODE.sequence])
    huNODE = huNODE.up

l_tmp = list()
while gpNODE.up is not common_ancestor:
    # print(gpNODE.name, gpNODE.sequence)
    l_tmp.append([gpNODE.name, gpNODE.sequence])
    gpNODE = gpNODE.up

lineage.extend(l_tmp[::-1])
mask_seq = lineage[-1][1]
for l in lineage:
    l[1] = ''.join([s for s, ms in zip(l[1], mask_seq) if ms != '-'])
    assert(len(l[1]) == len(mask_seq.replace('-', '')))
    l[1] = str(Seq(l[1], generic_dna).translate())


collapsed_lineage = [lineage[0]]
for l in lineage[1:]:
    if l[1] != collapsed_lineage[-1][1]:
        collapsed_lineage.append(l)
lineage = collapsed_lineage

for i, l in enumerate(lineage[1:-1]):
    l[0] = 'intermediate_' + str(i+1)


# print(lineage)
huAA = lineage[0][1]
for i, l in enumerate(lineage):
    similarity = sum(1 for a1, a2 in zip(huAA, l[1]) if a1 == a2) / float(len(huAA))
    if i != 0:
        sim2 = sum(1 for a1, a2 in zip(lineage[i-1][1], l[1]) if a1 == a2) / float(len(lineage[i-1][1]))
        print('>{} sim_to_human={} sim_to_previous={}\n{}'.format(l[0], similarity, sim2, l[1]))
    else:
        print('>{} similarity to human={}\n{}'.format(l[0], similarity, l[1]))








