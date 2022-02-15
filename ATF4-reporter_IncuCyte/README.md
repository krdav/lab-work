# ATF4 reporter assay using IncuCyte
Here lives data, data processing and plots related to an IncuCyte based ATF4 reporter assay.


### Method
Reporter used in Gu et al. 2021 (https://www.addgene.org/141281/)

After infection one or two stable expressing clones were chosen.
Often designated low/high indicating their baseline expression before induction.
ATF4 signal is measured in the RFP channel either as an endpoint or continuous measurements.
First pass processing is preformed by the IncuCyte software to subtract bachground fluorescence, then the integrated RFP signal per well is exported and further processed for normalization and plotting.

The "sample_annotations.xlsx" contains information about the content of each plate and other relevant information.


