# Mendelianization: Concentrating Polygenic Signal into a Single Causal Locus

The [Mendelianization algorithm](https://onlinelibrary.wiley.com/doi/10.1002/gepi.70053) uses summary z-statistics to learn weighted combinations of outcome variables for complex conditions (e.g., symptom dimensions in depression) so that each composite outcome is associated with a single causal locus.

All code was tested in R version 4.3.1.

# Installation

> library(devtools)

> install_github("ericstrobl/Mendelianization")

> library(Mendelianization)

# Sample Run

Download the real sample data in the _Data_ folder and place it in your working directory. The dataset consists of [Pan-UK Biobank](https://pan.ukbb.broadinstitute.org/downloads) summary z-statistics for quality-controlled chromosome 7 variants, evaluated across 52 outcomes that capture diverse dimensions of lifetime depression and anxiety.

> load("DepAnx_zstats.RData") # load sample chromosome 7 z-statistics (also load chromosomes and positions dataframe)

Run the Mendelianization algorithm:

> out = Mendelianization(Z,SoM=T,chr_pos)

The **Score of Mendelianism** (SoM) quantifies the extent to which the outcome optimized for a single variant concentrates at one LD locus, ranging from 0 to 1. Larger values correspond to stronger Mendelian-like (single locus) behavior. For this dataset, the expected result is perfect Mendelianism (SoM = 1) for a lead variant, since we have optimized over 52 diverse measures of depression and anxiety:

> print(out$SoMs) # print SoMs for leads

In addition, the algorithm produces well-calibrated p-values. As a diagnostic, the histogram of p-values should approximate a uniform distribution:

> hist(out$pval) # should be close to uniform

To extract the lead variant and examine its p-value (which incorporates outcome learning):

> leads = as.numeric(names(out$SoMs))

> print(out$pval[leads[1]]) # p-value of first lead variant

A Manhattan plot can also be constructed, which should exhibit a single prominent signal (``tower''), consistent with the SoM of 1:

> plot_tower(out,Z,leads[1]) # Manhattan plot of fixed outcome of the first lead variant applied genome-wide

Finally, the canonical coefficients associated with the lead variant (or any variant) are directly interpretable across outcomes. In other words, their magnitudes are comparable across different outcomes, thereby facilitating interpretation:

> out$Alpha_p[,leads[1]] # interpretable canonical coefficients of the first lead variant

# Updates

**8/27/2026:** We identified that p-values can be poorly calibrated when the outcome panel includes traits with diffuse associations across many variants genome-wide, such as quantitative blood-cell or neuroimaging traits. In this setting, associated variants can make a nonvanishing aggregate contribution to the genome-wide covariance estimator, violating Condition 2 (vanishing aggregate contamination) of Theorem 2. Selecting a sufficiently nonredundant subset of outcomes can also be burdensome when the outcome panel is large. We therefore added two automated options:

- `screen_Gamma = TRUE` uses adaptive central-null estimation of the covariance matrix `Gamma` to improve null calibration. Candidate trimming proportions are controlled by `threshold_grid`, which defaults to `seq(0, 0.20, length.out = 20)` but can be customized.
- `screen_outcomes = TRUE` iteratively removes redundant outcomes until the null-correlation condition number is at most `max_condition`. The default is `max_condition = 100`, but this threshold can also be customized.

The new default behavior is `screen_Gamma = TRUE`, `screen_outcomes = TRUE`, `threshold_grid = seq(0, 0.20, length.out = 20)`, and `max_condition = 100`. To recover the fully legacy behavior, set both screening options to `FALSE`.
