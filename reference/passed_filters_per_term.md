# Report filters applied to bootstrap-based enrichment results

This function filters enrichment results based on various criteria, such
as minimum intersection, significance thresholds, and bootstrapping
fractions, reporting which terms passed / didn't pass which filter.

## Usage

``` r
passed_filters_per_term(
  unfiltered_df,
  enrich_type = c("ORA", "MSEA"),
  min_intersection = 3,
  alpha_cutoff = 0.05,
  q.val_cutoff = 0.2,
  boot_fract_cutoff = 0.5
)
```

## Arguments

- unfiltered_df:

  A data frame containing unfiltered enrichment results from
  [`Run_bootstrap_ORA`](https://alexandrovteam.github.io/S2IsoMEr/reference/Run_bootstrap_ORA.md)
  or
  [`Run_bootstrap_MSEA`](https://alexandrovteam.github.io/S2IsoMEr/reference/Run_bootstrap_MSEA.md).

- enrich_type:

  A character string specifying the type of enrichment analysis. Must be
  either "ORA" (Over-Representation Analysis) or "MSEA" (Metabolite Set
  Enrichment Analysis).

- min_intersection:

  An integer specifying the minimum number of true positives (TP)
  required for a term to pass the filter.

- alpha_cutoff:

  A numeric value specifying the p-value cutoff for significance.
  Default is 0.05.

- q.val_cutoff:

  A numeric value specifying the q-value cutoff for significance.
  Default is 0.2.

- boot_fract_cutoff:

  A numeric value specifying the minimum fraction of bootstraps in which
  a term must be present to be considered.

## Value

A data frame with binary columns indicating which terms pass the
specified criteria. The data frame columns include:

- Term

- min_TP (minimum true positives)

- significant_adj_boot (significant adjusted bootstrap p-value)

- significant_adj_terms (significant adjusted term p-value)

- pass_boot_fraction (pass bootstrapping fraction)

- pass_all_filts (Whether term passes all filters)

## Details

The function adjusts p-values using the False Discovery Rate (FDR)
method and calculates combined p-values using the `metap_sumlog_pvals`
function. It then applies several filters to determine which terms pass
all criteria.

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
enrich_ORA_summary <- passed_filters_per_term(
  unfiltered_df = ORA_res$upregulated$unfiltered_enrich_res,
  enrich_type = "ORA"
)
} # }
```
