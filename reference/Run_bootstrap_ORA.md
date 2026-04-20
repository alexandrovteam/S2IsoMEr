# Wrapper function for Bootstrap Over-Representation Analysis (ORA)

This function performs a bootstrap-based over-representation analysis
(ORA) on a given list of metabolites and background.

## Usage

``` r
Run_bootstrap_ORA(
  marker_list,
  background,
  custom_universe = NULL,
  alpha_cutoff = 0.05,
  min_intersection = 3,
  consider_isobars = T,
  polarization_mode = NA,
  mass_range_ppm = 3,
  annot_db = "HMDB",
  annot_custom_db = NULL,
  use_LION = F,
  endogenous_only = T,
  pathway_assoc_only = F,
  remove_expected_predicted = T,
  annot_list = NULL,
  annot_weights = NULL,
  n_bootstraps = 50,
  boot_fract_cutoff = 0.5,
  q.val_cutoff = 0.2,
  selected_terms = NULL,
  adjust_contingency = T,
  report_ambiguity_scores = F
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

- consider_isobars:

  A logical indicating whether to consider isobars.

- polarization_mode:

  A parameter for polarization mode. If provided, it will be used to
  consider isomers and isobars.

- mass_range_ppm:

  A numeric value for mass range in parts per million (ppm).

- annot_db:

  A character string specifying the annotation database(s) used for
  annotation. Check
  [`build_iso_bg`](https://alexandrovteam.github.io/S2IsoMEr/reference/build_iso_bg.md)
  for more information.

- annot_custom_db:

  An optional custom annotation database. If provided, it will be used
  alongside or instead of the specified annotation database.

- use_LION:

  A logical indicating whether to use LION ontology for select
  isomers/isobars for background.

- endogenous_only:

  A logical indicating whether to consider only endogenous compounds.
  Applies only for HMDB and CoreMetabolome annotation databases.

- pathway_assoc_only:

  A logical indicating whether to consider only pathway-associated
  compounds. Applies only for HMDB and CoreMetabolome annotation
  databases.

- remove_expected_predicted:

  A logical indicating whether to remove expected and predicted
  annotations. Applies only for HMDB and CoreMetabolome annotation
  databases.

- annot_list:

  An optional list of annotations. If provided, it will be used instead
  of generating isomers from built-in formula to molecule mapping.

- annot_weights:

  An optional list of annotation weights. If provided, it will be used
  in the ambiguity score calculation and bootstrap analysis.

- n_bootstraps:

  An integer specifying the number of bootstrap iterations.

- boot_fract_cutoff:

  A numeric value specifying the minimum bootstrap fraction cutoff
  required. Default to 0.5 .

- q.val_cutoff:

  A numeric value specifying the q-value cutoff for significance.

- selected_terms:

  An optional vector of selected terms to focus on. If provided, the
  analysis will be restricted to these terms.

- adjust_contingency:

  A logical indicating whether to adjust the contingency table to
  account for isomeric ambiguity.

- report_ambiguity_scores:

  A logical indicating whether to report ambiguity scores. If TRUE,
  ambiguity scores will be included in the results.

## Value

A list containing the ORA results for each group of markers. Both
filtered and per-bootstrap results are provided.

## References

Badia-i-Mompel, P., Nagai, J. S., & Saez-Rodriguez, J. (2022).
decoupleR: A flexible tool to handle various modes of biological network
analysis. *Bioinformatics Advances*, 2(1), vbac016.
[https://doi.org/10.1093/bioadv/vbac016](https://academic.oup.com/bioinformaticsadvances/article/2/1/vbac016/6544613)

## Examples

``` r
if (FALSE) { # \dontrun{
data("example_ORA_markers", package = "S2IsoMErData")
data("example_ORA_obj", package = "S2IsoMErData")
object = example_ORA_obj
enrich_res <- Run_bootstrap_ORA(
  marker_list = example_ORA_markers,
  background = object$pathway_list,
  polarization_mode = object$polarization_mode,
  mass_range_ppm = object$mass_range_ppm,
  annot_db = object$Annotation_database,
  annot_custom_db = object$Custom_database,
  use_LION = ifelse(stringr::str_detect(object$background_name, "LION"), TRUE, FALSE),
  endogenous_only = object$endogenous_only,
  pathway_assoc_only = object$pathway_assoc_only,
  remove_expected_predicted = object$remove_expected_predicted,
  annot_weights = object$annotation.weights,
  consider_isobars = object$consider_isobars,
  annot_list = object$annotations
)
} # }
```
