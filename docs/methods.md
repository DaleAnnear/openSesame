# Methods and scientific limits

SeSAMe reads paired IDATs and applies the explicit configured preparation sequence. Detection and masking metrics are retained. Cohort matrices are constructed only from common probes and retain sample/probe order. DMPs use limma on M-values, with BH adjustment. For a two-group group contrast, reported beta-scale effects include both the covariate-adjusted beta-model coefficient (`delta_beta_adjusted`) and the observed difference of group mean beta values (`delta_beta_observed`). Linear-model coefficients are not constrained to the beta range; the observed delta beta is bounded to -1 through 1 and is used for the DMP effect-size threshold.

No biological-variable-associated probe removal is implemented. Sex discordance must be interpreted cautiously; no EPICv2-validated sex or fingerprint module is bundled yet. DMR execution is gated on the EPICv2 DMRcate annotation packages and is not allowed to fall back to EPICv1 coordinates.
