# Run Metabolite/Lipids enrichment analysis for single cell metabolomics

Run_enrichment() A wrapper to run either ORA or MSEA based on the
initialized enrichment object

## Usage

``` r
Run_enrichment(
  object,
  Run_DE = FALSE,
  DE_pval_cutoff = 0.05,
  DE_LFC_cutoff = 1,
  min.pct.diff = 0.1,
  ...
)
```

## Arguments

- object:

  A numeric matrix of n metabolites (rows) and m cells or measurments
  (columns).

- Run_DE:

  A logical indicating whether to run differential analysis using
  limma's rank Sum Test With Correlation. Ignored if enrichment type is
  'MSEA'.

- DE_pval_cutoff:

  A numeric indicating p-value cutoff. Default is 0.05.

- DE_LFC_cutoff:

  A numeric indicating the minimum log2 fold change for differential
  analysis

- min.pct.diff:

  A numeric indicating the minimum percentage difference between
  samples/cells in both conditions for a marker to be considered
  differentially abundant.

- ...:

  Additional arguments passed to functions used internally, such as
  [`Run_bootstrap_ORA`](https://alexandrovteam.github.io/S2IsoMEr/reference/Run_bootstrap_ORA.md),
  [`Run_simple_ORA`](https://alexandrovteam.github.io/S2IsoMEr/reference/Run_simple_ORA.md),
  [`Run_bootstrap_MSEA`](https://alexandrovteam.github.io/S2IsoMEr/reference/Run_bootstrap_MSEA.md),
  and
  [`Run_simple_MSEA`](https://alexandrovteam.github.io/S2IsoMEr/reference/Run_simple_MSEA.md).
  Consult the specific function documentation for details on these
  arguments.

## Value

Data.frame with enrichment results

## References

Badia-i-Mompel, P., Nagai, J. S., & Saez-Rodriguez, J. (2022).
decoupleR: A flexible tool to handle various modes of biological network
analysis. *Bioinformatics Advances*, 2(1), vbac016.
[https://doi.org/10.1093/bioadv/vbac016](https://academic.oup.com/bioinformaticsadvances/article/2/1/vbac016/6544613)

## Examples

``` r
if (FALSE) { # \dontrun{
data("example_ORA_obj", package = "S2IsoMErData")
data("example_ORA_custom_universe", package = "S2IsoMErData")
input_scm = example_ORA_obj$scmatrix
conds = example_ORA_obj$conditions
cond_x = "U"
cond_y = "F"
ORA_obj <- initEnrichment(
  scmatrix = input_scm,
  conditions = conds,
  enrichment_type = "ORA",
  annot_db = "HMDB",
  consider_isomers = TRUE,
  consider_isobars = TRUE,
  polarization_mode = "positive",
  background_type = "sub_class",
  molecule_type = "Metabo",
  condition.x = cond_x,
  condition.y = cond_y
)
ORA_res <- Run_enrichment(
  object = ORA_obj,
  custom_universe = example_ORA_custom_universe,
  report_ambiguity_scores = TRUE,
  DE_LFC_cutoff = 0,
  min.pct.diff = 0
)
} # }

```
