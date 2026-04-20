# Dot Plot for Over-Representation Analysis (ORA) Results

This function creates a dot plot for ORA results, with options to color
the dots by enrichment score or significance and to facet by a specified
variable.

## Usage

``` r
dotplot_ORA(
  ORA_res,
  alpha_cutoff = 0.05,
  color_by = "Significance",
  facet_by = NULL
)
```

## Arguments

- ORA_res:

  A data frame containing the ORA results from either filtered results
  in output of
  [`Run_bootstrap_ORA`](https://alexandrovteam.github.io/S2IsoMEr/reference/Run_bootstrap_ORA.md)
  or direct output of
  [`Run_simple_ORA`](https://alexandrovteam.github.io/S2IsoMEr/reference/Run_simple_ORA.md).

- alpha_cutoff:

  A numeric value specifying the alpha cutoff for significance. Default
  is 0.05.

- color_by:

  A character string indicating whether to color the dots by
  "Enrichment_score" or "Significance"

- facet_by:

  An optional character string specifying a column name by which to
  facet the plot.

## Value

A ggplot object displaying the dot plot of ORA results.

## Details

The function creates a dot plot where the size of the dots represents
the size of the term (e.g., number of true positives), and the color
represents either the enrichment score or the significance (-log10 of
the p-value). If multiple conditions are present, the plot is adjusted
accordingly. The plot can also be faceted by a specified variable.

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
multi_cond_res = collapse_ORA_boot_multi_cond(ORA_boot_res_list = ORA_res)
p = dotplot_ORA(ORA_res = multi_cond_res$clean_enrich_res)
} # }
```
