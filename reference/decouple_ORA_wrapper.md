# Run ORA Analysis with Optional Bootstrapping

This function performs over-representation analysis (ORA) on a given
list of markers using either a bootstrapped or a simple method. It is a
wrapper around the `run_ora` function from the `decoupleR` package. This
function is also used internally in the `Run_bootstrap_ORA` and
`Run_simple_ORA` functions.

## Usage

``` r
decouple_ORA_wrapper(
  marker_list,
  term_list,
  universe,
  pass_adjust = F,
  seed = 42,
  ORA_boot = T
)
```

## Arguments

- marker_list:

  A list of metabolites per condition. If there are no conditions,
  please provide it as list("condition" = c("metabo_1", "metabo_2",
  ...))

- term_list:

  A list of named terms where each term is a character vector of
  metabolites.

- universe:

  A character vector specifying the universe of metabolites.

- pass_adjust:

  A logical indicating whether to pass the adjustment parameter to the
  ORA function. If TRUE, the contingency table won't be adjusted. Only
  required for bootstrap-ORA.

- seed:

  An integer specifying the seed for reproducibility.

- ORA_boot:

  A logical indicating whether to use bootstrapped ORA. If TRUE, the
  function will use the bootstrapped ORA method; otherwise, it will use
  a simple ORA method.

## Value

A list containing two elements:

- ORA_res:

  A data frame with the ORA results, including columns for statistic,
  source, condition, score, and p-value.

- ORA_conting:

  A data frame with the ORA contingency table, including columns for
  source, condition, TP, FP, FN, and TN.

## Details

This function is a wrapper around the `run_ora` function from the
`decoupleR` package. For more details, please refer to the [decoupleR
documentation](https://saezlab.github.io/decoupleR/reference/run_ora.html).

## References

Badia-i-Mompel, P., Nagai, J. S., & Saez-Rodriguez, J. (2022).
decoupleR: A flexible tool to handle various modes of biological network
analysis. *Bioinformatics Advances*, 2(1), vbac016.
[https://doi.org/10.1093/bioadv/vbac016](https://academic.oup.com/bioinformaticsadvances/article/2/1/vbac016/6544613)

## Examples

``` r
if (FALSE) { # \dontrun{
marker_list <- c("B2", "B5", "B8")
term_list <- list(
  "Term1" = c("B1", "B2", "B3", "B4"),
  "Term2" = c("B4", "B5", "B6", "B9")
)
universe <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9")
results <- decouple_ORA_wrapper(
  marker_list, term_list, universe,
  pass_adjust = TRUE, seed = 42, ORA_boot = FALSE
)
} # }
```
