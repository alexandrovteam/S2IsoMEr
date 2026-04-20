# Plot Heatmap of MSEA Results Across Multiple Conditions

This function generates a heatmap of combined MSEA results. The heatmap
displays normalized enrichment scores (NES) across multiple conditions
and highlights significant results based on a specified alpha cutoff.

## Usage

``` r
plot_MSEA_Multi_cond(combined_MSEA_res, alpha_cutoff = 0.05)
```

## Arguments

- combined_MSEA_res:

  A data frame containing combined MSEA results for each pairwise
  comparison from either enrichment results in
  [`Run_bootstrap_MSEA`](https://alexandrovteam.github.io/S2IsoMEr/reference/Run_bootstrap_MSEA.md)
  or
  [`Run_simple_MSEA`](https://alexandrovteam.github.io/S2IsoMEr/reference/Run_simple_MSEA.md)

- alpha_cutoff:

  A numeric value specifying the cutoff for significance. Terms with
  p-values below this threshold will be marked as significant. Default
  is 0.05.

## Value

A heatmap plot showing the NES values for each term across different
condition pairs. Significant results are highlighted with an asterisk.

## Details

The function checks for the presence of required columns and ensures
there are enough conditions to generate a meaningful plot. The heatmap
includes color coding for NES values, with significant results indicated
by asterisks.

## Examples

``` r
if (FALSE) { # \dontrun{
data("example_MSEA_multicond", package = "S2IsoMErData")
p = plot_MSEA_Multi_cond(combined_MSEA_res = example_MSEA_multicond,
                         alpha_cutoff = 0.05)
} # }
```
