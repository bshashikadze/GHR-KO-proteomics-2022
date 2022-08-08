# GHR-KO-proteomics-2022

Reproducibility repo for Shashikadze et al. "Structural and proteomic repercussions of growth hormone receptor deficiency on the pituitary gland – lessons from a translational pig model". 

Maxquant protein groups file was first processed with Perseus. Proteins with valid values in at least 70% of all replicates in at least one condition (genotype) were quantified. Missing values were imputed from a normal distribution using the Perseus default parameters (width - 0.3 and down-shift - 1.8). The volcano plot was generated with the Perseus, using Student’s t-test and permutation-based FDR cut-off of 0.05, alongside with an s0-parameter of 0.1. [The imputed data with Perseus-generated statistics](https://github.com/ShashikadzeB/GHR-KO-proteomics-2022/blob/main/input%20files/perseus_output.txt) was used for downstream analysis in R. Figure 3 and figure 4 in the manuscript were generated with the following [code](https://github.com/ShashikadzeB/GHR-KO-proteomics-2022/blob/main/figure3_figure4.md).

Next, 2-way anova analysis was performed followed with Tukey's HSD - [R code to reproduce the results](https://github.com/ShashikadzeB/GHR-KO-proteomics-2022/blob/main/ANOVA%20analysis/ANOVA_THSD.md)


