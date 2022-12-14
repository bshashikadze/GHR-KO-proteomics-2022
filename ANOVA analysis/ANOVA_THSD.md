2 way ANOVA with Tukey’s HSD and Cohen’s D
================
BS
07/08/2022

## load libraries

``` r
library(tidyverse)
```

## load data

``` r
perseus_data        <- read.csv("perseus_output.csv", check.names = F)
conditions          <- read.delim("conditions.txt")
conditions$Animal   <- as.character(conditions$Animal)
```

## 2 way anova

``` r
# perform two way ANOVA
anova  <- perseus_data %>% 
          pivot_longer(names_to = "Animal", 
          values_to = "Abundance", -Gene) %>% 
          left_join(conditions) %>% 
          group_by(Gene) %>% 
          summarise(anova_p_val = 
          summary(aov(Abundance ~ Genotype*Sex))[[1]][["Pr(>F)"]][1:3]) 

# correct all resulting p-values (pool) for multiple hypothesis testing
anova$adjusted_p_val  <- p.adjust(anova$anova_p_val, method = "BH")

# prepare empty data frame with proper comparisons
anova_results <- as.data.frame(rep(c("genotype", "sex", "genotype:sex"), length = length(anova$Gene)))

# rename column 
names(anova_results) <- "Comparison"
```

## data frame with anova statistics

``` r
# final anova results
anova_results  <- as.data.frame(cbind(anova, anova_results))

# long to wide table format
anova_results <- anova_results %>%
                 pivot_wider(names_from = Comparison, 
                             values_from = c(anova_p_val, adjusted_p_val), Gene) %>% 
                 left_join(perseus_data)
```

## Tukey’s HSD for significant interactions (genotype:sex adjusted p-val \<0.05)

``` r
# filter data for significant interactions
data_tukey <- anova_results %>% 
              filter(`adjusted_p_val_genotype:sex` < 0.05) %>% 
              select(Gene, 8:15) %>% 
              pivot_longer(names_to = "Animal", values_to = "Abundance", 
                           -Gene) %>% 
              left_join(conditions)



# for loop. for each gene which showed significant interaction from
# 2 way anova, THSD is calculated. significant pairs are extracted

# make an empty list, where to each gene significant interactions will be assigned
my_vec <- list()

# for which genes
int_sig_genes <- unique(data_tukey$Gene)

# for loop
for (i in int_sig_genes) {
    sub_df              <- data_tukey[data_tukey$Gene %in% i,]
    anova_model         <- aov(data= sub_df,  Abundance ~Genotype*Sex)
    anova_tukey         <- TukeyHSD(anova_model)
    tuk_interactions    <- anova_tukey[3][[1]] %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    filter(`p adj` < 0.05) %>% 
    select(rowname)
    tuk_interactions_t  <- t(tuk_interactions)
    sig_tuk_interaction <- matrix(apply(tuk_interactions_t,1,paste,collapse=";"),
                                nrow=1)
    my_vec[i] <- sig_tuk_interaction
  }


# final data with anova and THSD statistics
anova_with_TKHSD <- anova_results %>% 
               left_join(my_vec %>% 
               unlist() %>% 
               as.data.frame() %>% 
               rownames_to_column() %>% 
               rename(Gene=1) %>% 
               rename(THSD_pair = 2))
```

## effect size calculation (Cohens’d)

``` r
anova_effect1 <- perseus_data %>% 
                 pivot_longer(names_to = "Animal", values_to = "Abundance", 
                           -Gene) %>% 
                 left_join(conditions) %>% 
                 group_by(Gene, Genotype) %>% 
                 summarise(mean=mean(Abundance), sd = sd(Abundance)) %>% 
                 pivot_wider(names_from = Genotype, 
                             values_from = c(mean, sd), Gene) %>% 
                 mutate(Cohens.d_genotype = 
                (`mean_GHR-KO` - mean_Control)/sqrt((`sd_GHR-KO`^2 +   sd_Control^2)/2))


anova_effect2 <- perseus_data %>% 
                 pivot_longer(names_to = "Animal", values_to = "Abundance", 
                           -Gene) %>% 
                 left_join(conditions) %>% 
                 group_by(Gene, Sex) %>% 
                 summarise(mean=mean(Abundance), sd = sd(Abundance)) %>% 
                 pivot_wider(names_from = Sex, values_from = c(mean, sd), Gene) %>%
                 mutate(Cohens.d_sex = 
                (mean_Male - mean_Female)/sqrt((sd_Male^2 +    sd_Female^2)/2))
```

## final data with anova, THSD, and effect size

``` r
final_data <- anova_effect1 %>% 
  left_join(anova_effect2) %>% 
  left_join(anova_with_TKHSD)
```

``` r
write.table(final_data, "anova_results.csv", row.names = F, sep = ",")
```
