# Wrapper function for Simple Over-Representation Analysis (ORA)

This function performs a simple over-representation analysis (ORA) on a
given list of metabolites and background.

## Usage

``` r
Run_simple_ORA(
  marker_list,
  background,
  custom_universe = NULL,
  alpha_cutoff = 0.05,
  min_intersection = 3
)
```

## Arguments

- marker_list:

  A vector or list of marker sets to analyze.

- background:

  A list of named terms where each term is a character vector of
  metabolites.

- custom_universe:

  An optional vector specifying a custom universe of terms. If not
  provided all metabolites in background will be used as universe.

- alpha_cutoff:

  A numeric value indicating the alpha cutoff for significance.

- min_intersection:

  An integer specifying the minimum intersection required between input
  list and a given term in background.

## Value

A dataframe containing the ORA results for each group of markers.

## References

Badia-i-Mompel, P., Nagai, J. S., & Saez-Rodriguez, J. (2022).
decoupleR: A flexible tool to handle various modes of biological network
analysis. *Bioinformatics Advances*, 2(1), vbac016.
[https://doi.org/10.1093/bioadv/vbac016](https://academic.oup.com/bioinformaticsadvances/article/2/1/vbac016/6544613)

## Examples

``` r
if (FALSE) { # \dontrun{
data("example_ORA_markers", package = "S2IsoMErData")
bg = Load_background(mol_type = "Metabo",bg_type = "main_class",feature_type = "sf")
enrich_res = Run_simple_ORA(marker_list = example_ORA_markers,background = bg)
} # }
```
