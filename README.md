# GHR-KO-proteomics-2022

Reproducibility repo for Shashikadze et al. "Structural and proteomic repercussions of growth hormone receptor deficiency on the pituitary gland – lessons from a translational pig model". 

Maxquant [proteingroups](https://github.com/ShashikadzeB/GHR-KO-proteomics-2022/blob/main/input%20files/proteinGroups.txt) file was first processed with perseus. Proteins detected in at least 70% of all replicates in at least one condition were kept for quantitative analysis. To handle missing values, data imputation from a normal distribution was performed using the Perseus default parameters (width - 0.3 and down-shift - 1.8). The volcano plots were generated with the Perseus using two-tailed Student’s t-test statistics and permutation-based FDR cut-off of 0.05, together with an s0-parameter of 0.1 to additionally consider fold changes. Figure 3 and figure 4 in the manuscript were generated with the following [code](https://github.com/ShashikadzeB/GHR-KO-proteomics-2022/blob/main/figure3_figure4.md).

2 way anova analysis was performed followed with Tukey's HSD, [R code to reproduce the results](https://github.com/ShashikadzeB/GHR-KO-proteomics-2022/blob/main/ANOVA%20analysis/ANOVA_THSD.md)


The dataset has been submitted to the ProteomeXchange Consortium via the PRIDE partner repository (PXD$$$$$). 
