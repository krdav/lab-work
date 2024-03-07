# Amino acid influx
Here, I determine the influx of amino acids in aspartate auxotroph 143B expressing SLC1A3 with GOT1 and GOT2 knocked out.
See the notebook for more details on the experimental setup and on how to calculate flux.

Carefull inspection of the output data will lead one to wonder where the cell disposes all these amino acids.
Typically, 143B cells have 3-400 pg intracellular protein per cell at 4500-5500 fL cell volume and with a proteome based amino acid distribution that gives a sum of ~0.7 M on all amino acids.
The influx data presented here show the total amino acid accumulation in a cell to be ~3.5 M.
The discrepancy could be due to the large extracellular matrix secreted by 143B cells, as this is not accounted for in the 3-400 pg protein per cell.

Measurement error, wrongly noted dilution or other trivial errors could also explain it; however, the glutamine (and serine) influx is within range of previously measured and published values:
1. Sullivan et al 2018 Nat. Cell Biol. report "anaplerotic glutamine consumption" (glutamine consumption rate – glutamate production rate) in 143B cells to be ~170 fmol/(cell\*h).
This translates into ~38 mM/h vs. 54 mM/h presented in this dataset (at 31 h sampling time).
2. Fan et. al 2013 Mol. Syst. Biol. report 40 mM/h glutamine influx for parental iBMK cells and 54 and 45 mM/h for RAS and Akt activate iBMK, respectively, vs. ~70 mM/h presented in this dataset (at 31 h sampling time).
3. Park et al. 2016 Nat. Chem Biol. report 53 mM/h glutamine influx for iBMK cells (supposedly same cells as Fan et al.^.) vs. ~70 mM/h presented in this dataset (at 31 h sampling time).
Park et al also report serine uptake to be 8 mM/h vs ~9 mM/h presented in this dataset (at 23 h sampling time).

According to Fan et al. iBMK have a doubling time of 24 h.
Doubling time for 143B is ~20 h.
This could further explain the difference in glutamine influx as faster proliferation leads to faster nutrient influx.
In summary, the values presented here appear high but are supported by previous published values.



