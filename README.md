# Mendelianization: Concentrating Polygenic Signal into a Single Causal Locus

The Mendelianization algorithm uses summary z-statistics to learn weighted combinations of outcome variables for complex conditions (e.g., symptom dimensions in depression) so that each composite outcome is associated with a single causal locus.

# Installation

> library(devtools)

> install_github("ericstrobl/Mendelianization")

> library(Mendelianization)

# Sample Run

Download the real sample data in the _Data_ folder and place it in your working directory:

> load("DepAnx_zstats.RData") # load sample chromosome 7 z-statistics (also load chromosomes and positions dataframe)

Run the Mendelianization algorithm:

> out = Mendelianization(Z,SoM=T,chr_pos)

The **Score of Mendelianism** (SoM) quantifies the extent to which the lead variant follows Mendelian inheritance, ranging from 0 to 1. Larger values correspond to stronger Mendelian behavior. For this dataset, the expected result is perfect Mendelianism (SoM = 1):

> print(out$SoMs) # print SoMs for leads

In addition, the algorithm produces well-calibrated p-values. As a diagnostic, the histogram of p-values should approximate a uniform distribution:

> hist(out$pval) # should be close to uniform

To extract the lead variant and examine its p-value (which incorporates outcome learning):

> leads = as.numeric(names(out$SoMs))

> print(out$pval[leads[1]]) # p-value of first lead variant

A Manhattan plot can also be constructed, which should exhibit a single prominent signal (``tower''), consistent with an SoM of 1:

> plot_tower(out,Z,leads[1],chr_pos) # Manhattan plot of fixed outcome of the first lead variant applied genome-wide

Finally, the canonical coefficients associated with the lead variant (or any variant) are directly interpretable across outcomes. In other words, their magnitudes are comparable across different outcomes, thereby facilitating interpretation:

> out$Alpha_p[,leads[1]] # interpretable canonical coefficients of the first lead variant
