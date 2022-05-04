# Container for quantitative work on 143B and H1299 cell lines
Various quantifications of metabolites/pathways related to aspartate metabolism in 143B and H1299 cell lines.
The folders contain data and analysis scripts for the following quantifications:


### `acid-hydrolysis`
Acid hydrolysis followed by LCMS to quantify total intracellular level of aspartate and nucleotides in 143B and H1299 cell lines.

### `cell_composition_prediction`
Contains all relevant information to predict cell composition to check the acid hydrolysis results.
Three subfolders store data/methods used to perform the prediction:

#### `intracellular_amino-acids`
LCMS based quantification of intracellular concentration of amino acids in 143B and H1299 cell lines.
Quantification is based on U-15N, U-13C internal standards and calibration curves with three replicates and either with/without 1 mM asparagine in the medai.

#### `BCA-prot-quant`
BCA based quantification of total protein in 143B and H1299 cell lines.
BCA is known to be a quick and robust way of quantifying whole cell protein; however, to absolute quantity is based on a BSA calibration curve which is not repressentative of all proteins in a cell.
Thus, the BCA can systematically over or under-estimate total protein.

#### `amino_acid_abundance_estimate`
A computational estimate of amino acid abundance based on global proteomics to determine the relative abundance of each protein in the proteome and the amino acid sequence of each protein.


### `AA_influx_143B-GOT-DKO`
Amino acid uptake flux into 143B cells with both GOT1 and GOT2 knocked out using LCMS based quantification.
The GOT double knockout makes the cells aspartate auxotrophic and expression of the aspartate transporter SLC1A3 enables them to sustain proliferation with 1 mM aspartate in the media.
In these cells the total aspartate uptake flux is the sum of all aspartate consuming reactions.

### `Asn-consumption-flux_143B`
Using labelled asparagine the asparagine consumption flux is determined using U-15N, U-13C internal standards and LCMS.




