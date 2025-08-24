# 2. Background Theory: Olink Proteomics
This section provides a brief overview of the technology that generates the data analyzed by the app.

## 2.1 What is Proteomics?

*Proteomics* is the large-scale study of proteins. It's a fundamental field of study in biology that aims to understand the structure, function, and interactions of proteins, as well as how they change in response to different conditions. Proteins are the workhorses of the cell, and their levels and activity provide a dynamic snapshot of biological processes.

## 2.2 Proximity Extension Assay (PEA)

The Olink platform is based on the Proximity Extension Assay (PEA) technology. This highly sensitive and specific method enables the simultaneous measurement of hundreds of proteins from a very small sample volume (typically 1-3 microliters of plasma, serum, or other biological fluid). The core principle is as follows:

1. **Antibody Pairs:** Each protein of interest is targeted by a pair of antibodies, each labeled with a unique DNA oligonucleotide.
2. **Proximity:** If both antibodies successfully bind to the same target protein, their attached DNA oligos are brought into close proximity.
3. **Hybridization and Extension:** The close proximity allows the oligos to hybridize and serve as a template for a DNA polymerase. The polymerase extends the hybridized strands, creating a new, unique DNA sequence (a "reporter sequence") for that specific protein.
4. **Quantification:** The newly formed reporter sequences are then quantified using high-throughput real-time PCR (qPCR), providing a readout proportional to the amount of the target protein in the original sample. This value is reported as Normalized Protein Expression (NPX), which is on a log2 scale.

The PEA technology's high specificity, as it requires two antibodies to bind to the same protein, minimizes false positive signals and allows for robust, multiplexed protein analysis.

## References

Assarsson, E. et al. (2014). "A single-tube, quantitative technique for high-throughput protein analysis." Nature Methods, 11(6), 665–670.

Olink Proteomics Official Website: https://www.olink.com/