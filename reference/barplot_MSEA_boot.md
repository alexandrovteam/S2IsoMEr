# Plot bootstrap enrichment analysis

Plot bootstrap enrichment analysis

## Usage

``` r
barplot_MSEA_boot(
  object,
  min.annotations = 2,
  q.value.cutoff = 0.1,
  bootstrap.fraction.cutoff = 0.5,
  by.statistic = "ES"
)
```

## Arguments

- object:

  A S2IsoMEr object after enrichment analysis.

- min.annotations:

  An integer describing the minimal number of annotations each term
  should include

- q.value.cutoff:

  A numeric between 0 and 1. Only terms with q-values lower than this
  value will be displayed.

- bootstrap.fraction.cutoff:

  A numeric between 0 and 1 (default = 0.5), indicating the minimal
  fraction that the metabolite set is present in all bootstrap
  iterations.

- by.statistic:

  A character indicating how the x-axis will be arranged. Can be either
  'ES' (enrichment score) or 'q.value' (default).

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE) { # \dontrun{
data("example_MSEA_obj", package = "S2IsoMErData")
p = barplot_MSEA_boot(object = example_MSEA_obj,
      q.value.cutoff = 0.2, by.statistic = "ES")
} # }
```
