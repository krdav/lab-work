# Stepwise procedure to make ancestral reconstruction of human asparaginase (ASNase)

### Download proteins with the typical human ASNase gene architecture (L-asparaginase, N-terminal domain followed by an Ankyrin repeat-containing domain):
`wget -O ASNase_human_like.fasta 'http://www.ebi.ac.uk/interpro/entry/IPR027474/proteins-matched?ida=495&export=fasta'`


### Extract protein IDs:
`grep '>' ASNase_human_like.fasta | tr -d '>' > ASNase_human_like.ids`


### Convert protein IDs to "RefSeq Nucleotide" IDs via UniProt:
https://www.uniprot.org/uploadlists/


### Download DNA sequences as Genbank file via NCBI:
https://www.ncbi.nlm.nih.gov/sites/batchentrez


### Extract the coding DNA sequence and translate to amino acids. Then filter out sequences that are less than 50% similar (according to Levenstein distance) to either 'human', 'guinea_pig', 'zebrafish' or 'fruit_fly':
`python ./extract_CDS_from_GBfile.py > ASNase_human_like_aa_seq_id50.fasta`


### Make a multiple sequence alignment on the selected protein sequences:
https://www.ebi.ac.uk/Tools/msa/mafft/
Download alignment as fasta and rename it XXX


### Rename sequence entries:
`python ./rename.py`


### Reverse translate amino acid alignment into DNA:
`python ./aa2dna_align.py >ASNase_human_like_dna_seq_id50_MAFFTalign.fasta`


### Convert to PHYLIP alignment format:
http://sequenceconversion.bugaco.com/converter/biology/sequences/fasta_to_phylip.php
Maybe this can be skipped by using fasta directly to IQ-TREE


### Perform phylogenetic and ancestral sequence reconstruction using IQ-TREE:
some command


### Create the maximum a posteriori sequences for each ancestor and map it back to the tree:
python3 ./map_ancestors.py






* Make ancestral sequence lineage between human and guinea pig

* Apply "red colored regions" on 3d structure as a mask i.e. discard all changes from human to guinea pig if outside of these regions.

* Choose three variant along this lineage.

* Compare the consensus sequence of the low vs. high Km ASNases to infer which residues are conserved in the low Km variant.

* Compare ^ conserved residues to the three chosen reconstructed variants.

* Run each variant through NetMHCIIpan, remove endogeneous peptides from human ASNase and estimate immunogenicity.



*** Variants to test ***
1. Loop reversion and conserved residues. Simply revert all ASN toutching loops and conserved regions to gpASNase.
2. Mutations from intermediate_3 -> gpASNase that overlap with the ASN toutching loops or conserved regions.
3. Take all 90% consensus positions and from "low_Km_variants_AA_MAFFTalign.fasta" and transplant these into huASNase.

4. As a positive control the variant suggested by Rigouin et al.

@NOTE Remember to revert the gaps of the C-term gpASNase to the full huASNase C-term.
***




### Low Km (active) ASNases
fruit fly, zebrafish, guinea pig, Aspergillus_terreus


### High Km:
Human, Aspergillus_nidulans, mouse




