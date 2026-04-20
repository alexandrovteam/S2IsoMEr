# Bar Plot for ORA Bootstrapped Results

This function creates a bar plot for Over-Representation Analysis (ORA)
bootstrapped results, including error bars to show the range of
enrichment scores (OR) and indicating the significance based on
q-values.

## Usage

``` r
barplot_ORA_boot(ORA_boot_res, collapse_multi_cond = F)
```

## Arguments

- ORA_boot_res:

  A list containing the ORA bootstrapped results. The list should
  include `unfiltered_enrich_res` and `clean_enrich_res` data frames.

- collapse_multi_cond:

  A logical value indicating whether to collapse results across multiple
  conditions. Default is `FALSE`.

## Value

A ggplot object displaying the bar plot of ORA bootstrapped results.

## Details

The function creates a bar plot where the terms are reordered by their
median enrichment scores. Bars are filled with colors representing the
combined q-values, and error bars show the minimum and maximum
enrichment scores (ES_min and ES_max) across bootstraps. The number of
true positives (n) is included in the term labels.

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
p = barplot_ORA_boot(ORA_boot_res = ORA_res,collapse_multi_cond = T)
} # }
```
