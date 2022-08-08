## 2-way ANOVA

[R script to reproduce the results](https://github.com/bshashikadze/GHR-KO-proteomics-2022/blob/main/ANOVA%20analysis/ANOVA_THSD.md)

[input](https://github.com/bshashikadze/GHR-KO-proteomics-2022/blob/main/ANOVA%20analysis/perseus_output.csv) should (at least) contain unique identifier (gene/protein name) and protein intensities for each sample. In this example samples are columns and rows are gene names.

[conditions file](https://github.com/bshashikadze/GHR-KO-proteomics-2022/blob/main/ANOVA%20analysis/conditions.txt) provides grouping information for each sample. 

Steps:
1. Perform 2-way ANOVA for each gene
2. Pool all resulting p-values and correct for multiple hypothesis testing
3. In this study, each factor contains two levels only, therefore Tukey's HSD was performed on significant interaction hits (adjusted p-value < 0.05)

[result](https://github.com/bshashikadze/GHR-KO-proteomics-2022/blob/main/ANOVA%20analysis/anova_results.csv)
