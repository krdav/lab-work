from ete3 import Tree

tree = 'ASNase_tree/ASNase_human_like_aa_seq_id50_MAFFTalign_rename.phylip.treefile'
ids = 'ASNase_human_like_aa_seq_id50_MAFFTalign_rename.ids'

t = Tree(tree, format=1)

with open(ids) as fh:
    for l in fh.readlines():
        source, ID = l.split()
        source = source.replace('(', '@')
        source = source.replace(')', '@')
        for leaf in t:
            if leaf.name == ID:
                leaf.name = source

print(t)
t.write(format=1, outfile='ASNase_tree/ASNase_human_like_aa_seq_id50_MAFFTalign.phylip.treefile')

