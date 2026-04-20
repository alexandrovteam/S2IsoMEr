# Ridge Plot of Bootstrap Enrichment Results for Terms of Interest

This function creates a ridge plot visualizing the distribution of
intersection sizes (TP) for specified terms of interest from bootstrap
enrichment results.

## Usage

``` r
ridge_bootstraps(enrich_res, terms_of_interest, condition = NULL)
```

## Arguments

- enrich_res:

  A data frame containing the enrichment results per bootstrap,
  including columns for Term, n (TP), and fraction.

- terms_of_interest:

  A character vector specifying the terms of interest to be plotted.

- condition:

  Optional. A character string specifying a particular condition to
  filter the enrichment results. If `NULL` (the default), results from
  all conditions will be included.

## Value

A ggplot object displaying the ridge plot of intersection sizes for the
specified terms of interest.

## Details

The function filters the enrichment results to include only the
specified terms of interest and creates a ridge plot showing the
distribution of intersection sizes (TP). The terms are labeled with
their names and the fraction of bootstraps in which they appear, and the
plot includes quantile lines to indicate distribution quartiles.

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
p4 = ridge_bootstraps(enrich_res = multi_cond_res$unfiltered_enrich_res,
                     terms_of_interest = c(multi_cond_res$clean_enrich_res$Term),
                     condition = "upregulated")
} # }
```
