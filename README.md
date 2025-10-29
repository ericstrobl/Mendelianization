# Mendelianization: Concentrating Polygenic Signal into a Single Causal Locus

The Mendelianization algorithm uses summary z-statistics to learn weighted combinations of outcome variables for complex conditions (e.g., symptom dimensions in depression) so that each composite outcome is associated with a single causal locus.

# Installation

> library(devtools)

> install_github("ericstrobl/Mendelianization")

> library(Mendelianization)

# Sample Run

Download the real sample data in the Data folder and place it in your working directory.

> load("DepAnx_zstats.RData") # load sample chromosome 7 z-statistics (also load chromosomes and positions dataframe)

> out = Mendelianization(Z,SoM=T,chr_pos) # with Score of Mendelianization

> print(out$SoMs) # print SoMs for leads

> leads = as.numeric(names(out$SoMs))

> print(out$pval[leads[1]]) # p-value of first lead variant

> plot_tower(out,Z,leads[1],chr_pos) # Manhattan plot of fixed outcome of the first lead variant applied genome-wide

> out$Alpha_p[,leads[1]] # interpretable canonical coefficients of the first lead variant
