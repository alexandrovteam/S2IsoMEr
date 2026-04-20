# Perform Bootstrap Sampling on Metabolite Annotations

This function performs bootstrap sampling on metabolite annotations,
allowing for resampling of molecules with or without associated weights.
It generates multiple bootstrap samples for further analysis.

## Usage

``` r
metabo_bootstrap(annot_list, annot_weights = NULL, n_bootstraps = 50)
```

## Arguments

- annot_list:

  A list of character vectors where each vector contains metabolites or
  molecules for a particular annotation. Each element in the list
  corresponds to a specific sum formula.

- annot_weights:

  A list of numeric vectors where each vector contains weights
  associated with the metabolites in `annot_list`. Must have the same
  length as `annot_list`. If NULL, all metabolites are equally weighted.

- n_bootstraps:

  An integer specifying the number of bootstrap samples to generate.

## Value

A list of length `n_bootstraps`, where each element is a list of
bootstrap samples for each sum formula. Each bootstrap sample is a
character vector of sampled metabolites.

## Details

This function is used to create multiple bootstrap samples of metabolite
annotations by randomly sampling metabolites from each sum formula. If
weights are provided, the sampling is performed according to these
weights; otherwise, the sampling is uniform.

## Examples

``` r
if (FALSE) { # \dontrun{
annot_list <- list(
  Sf1 = c("Metabolite1", "Metabolite2", "Metabolite3"),
  Sf2 = c("MetaboliteA", "MetaboliteB")
)
annot_weights <- list(
  Sf1 = c(0.1, 0.2, 0.7),
  Sf2 = c(0.5, 0.5)
)
bootstrapped_samples <- metabo_bootstrap(annot_list, annot_weights, n_bootstraps = 100)
} # }
```
