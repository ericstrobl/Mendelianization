# Mendelianization: Concentrating Polygenic Signal into a Single Causal Locus

An algorithm that uses summary z-statistics to learn weighted combinations of outcome variables for complex conditions (e.g., symptom dimensions in depression) so that each composite outcome is associated with a single causal locus.

# Installation

> library(devtools)
> install_github("ericstrobl/Mendelianization")
> library(Mendelianization)

# Sample Run

> load("DepAnx_zstats.RData") # load z-statistics (also load chrosomes and positions dataframe for Score of Mendelianization)
> out = Mendelianization(Z,SoM=F) # without Score of Mendelianization
> out_SoM = Mendelianization(Z,SoM=T,chr_pos) # without Score of Mendelianization
> hist(out$pval) # p-values
> out$Alpha_p[,1] # interpretable canonical coefficients for first variant
