# Updates

**8/27/2026:** We identified that p-values can be poorly calibrated when the outcome panel includes traits with diffuse associations across many variants genome-wide, such as quantitative blood-cell or neuroimaging traits. In this setting, associated variants can make a nonvanishing aggregate contribution to the genome-wide covariance estimator, violating Condition 2 (vanishing aggregate contamination) of Theorem 2. Selecting a sufficiently nonredundant subset of outcomes can also be burdensome when the outcome panel is large. We therefore added two automated options:

- `screen_Gamma = TRUE` uses adaptive central-null estimation of the covariance matrix `Gamma` to improve null calibration. Candidate trimming proportions are controlled by `threshold_grid`, which defaults to `seq(0, 0.20, length.out = 20)` but can be customized.
- `screen_outcomes = TRUE` iteratively removes redundant outcomes until the null-correlation condition number is at most `max_condition`. The default is `max_condition = 100`, but this threshold can also be customized.

The new default behavior is `screen_Gamma = TRUE`, `screen_outcomes = TRUE`, `threshold_grid = seq(0, 0.20, length.out = 20)`, and `max_condition = 100`. To recover the fully legacy behavior, set both screening options to `FALSE`.

**10/18/2026:** We added `Marginal_p` to the output of Mendelianization to facilitate interpretation of locus-specific phenotypes. This matrix contains the marginal variant–outcome associations, standardized by the estimated outcome standard deviations: `Zstat[j,a] / sqrt(Gamma[a,a])`. Unlike `Alpha_p,` which represents associations with residualized (unique) outcome components, `Marginal_p` retains the original marginal associations without adjusting for correlations among outcomes.
