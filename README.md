# ENCODE TE project (ENCOTE for short)
Various code and files for "Regulatory Transposable Elements in the Encyclopedia of DNA Elements"

## ENCOTE public datasets
Files containing links to publicly available datasets/files used for analysis.

## ENCOTE supplementary code
Unix commands, python scripts, and R scripts to perform analyses and create visualizations/plots in "Regulatory Transposable Elements in the Encyclopedia of DNA Elements" are provided here.

Code is grouped into different directories based on common themes of analysis or by individual writing the code. The primary purpose/aim of each directory is listed below.

Commands are generally listed in a commands.sh file. The exception is for the "combined_MPRA_feature_comparison_and_TF_distance_turnover" directory, where commands are listed in the jc_code.txt file.

For some commands, input files are obtained through commands from other directories. As such, commands in each directory should be run according to the following order.

### Reference files
Create reference files from public files to be used to subsequent analyses.

### TE-derived cCREs in human
Quantify TE-derived cCREs in humans.

### Human-mouse cCRE comparison
Compare human and mouse cCREs to identify shared and lineage-specific cCREs, using the information to quantify TE contributions to cCREs after human-mouse divergence.

### TF TE origins
Quantify cCRE associated transcription factor (TF) origins in TEs at TE subfamily level.

### Relative distance TE-cCREs
Quantify genomic distance of TEs (cCRE-associated and non-associated) to non-TE cCREs.

### Combined MPRA feature comparison and TF distance turnover
Several analyses combined together:
1. Quantify genomic distance of TEs (TF-bound and non-bound) to non-TE TF binding sites.
2. Compare MPRA activity of (primarily) TE-derived sequences and non-TE sequences using ENCODE K562 lentiMPRA data.
3. Compare feature overlap of TE-derived cCREs and non-TE cCREs. Features are K562 lentiMPRA activity, K562 TF ChIP-seq peak, ATAC-seq peak, and phastCons score.
4. Quantify TF binding site turnover between analogous human (K562) and mouse (MEL) cell lines.

### Compare TE non-TE human population
Quantify common human population variants (>1% allele frequency) in TE-derived cCREs compared to flanking sequence and non-TE cCREs.

### Compare TE non-TE GWAS
Quantify GWAS SNPs in TE-derived cCREs and compare to all cCREs and non-TE cCREs.
