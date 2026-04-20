# Generate S2IsoMEr enrichment object

initEnrichment() creates object to perform bootstrapping metabolite set
enrichment analysis

## Usage

``` r
initEnrichment(
  scmatrix,
  conditions,
  enrichment_type = c("ORA", "MSEA"),
  annot_db = "HMDB",
  annot_custom_db = NULL,
  endogenous_only = T,
  pathway_assoc_only = F,
  remove_expected_predicted = T,
  annotations = NULL,
  annotation.weights = NULL,
  consider_isomers = TRUE,
  consider_isobars = FALSE,
  mass_range_ppm = 3,
  polarization_mode = "positive",
  include = NULL,
  molecule_type = c("Lipid", "Metabo"),
  background_type = c("LION", "main_class", "super_class", "sub_class", "pathways"),
  custom_bg = NULL,
  termsOfInterest = "all",
  condition.x = NULL,
  condition.y = NULL,
  ranking.by = c("wilcox.test", "t.test", "BWS", "logFC"),
  gsea.method = c("fgsea", "ks_signed")
)
```

## Arguments

- scmatrix:

  A numeric matrix of n metabolites (rows) and m cells or measurments
  (columns).

- conditions:

  A vector of length m with condition identifiers.

- enrichment_type:

  A character specifying whether Overrepresentation analysis (ORA) or
  Metabolite Set Enrichment Analysis (MSEA) is performed

- annot_db:

  A character or character vector specifying which annotation
  database(s) were used. Current databases include ("CoreMetabolome",
  "HMDB","SwissLipids","LipidMaps"). Multiple databases can be selected

- annot_custom_db:

  An optional dataframe from which the isomers will be sampled. It
  should have 2 columns : formula and molecule name. If provided,
  `annot_db` will be ignored.

- endogenous_only:

  A logical indicating whether to only consider endogenous metabolites
  (default = TRUE).

- pathway_assoc_only:

  A logical indicating whether to only consider metabolites associated
  with a biological pathway (default = FALSE)

- remove_expected_predicted:

  A logical indicating whether to remove expected and predicted isomers
  based on HMDB status (default = TRUE)

- annotations:

  An optional custom list of length n, with each element contains a
  vector of isomer names. If not specified, S2IsoMEr uses the
  CoreMetabolome, LIPIDMAPS, SwissLipids, and HMDB databases from
  METASPACE (https://metaspace2020.eu/) to generate an annotation list
  automatically.

- annotation.weights:

  An optional list of length n, each element contains a vector of isomer
  weights. Only when annotations is provided as list.

- consider_isomers:

  A logical indicating whether to include isomers (default = TRUE)

- consider_isobars:

  A logical indicating whether to include isobars (default = FALSE)

- mass_range_ppm:

  A numeric indicating the mass range in ppm (default: mass_range_ppm =
  3). Molecular formulas + adducts within this range will be treated as
  isobars. Only required when isobars = TRUE.

- polarization_mode:

  A character with either 'positive' (default) or 'negative'. Only
  required when isobars = TRUE. When set to 'positive', included adducts
  are '+H', '+Na', and '+K'. When set to 'negative', included adducts
  are '-H', '+Cl'.

- include:

  An optional logical vector of length n indicating whether to include
  the annotations in the analysis.

- molecule_type:

  A character specifying the feature type for the enrichment background,
  either Metabolites or Lipids. Valid choices are c("Metabo", "Lipid").

- background_type:

  A character specifying the background type for enrichment, choose one
  from c("LION","main_class","super_class","sub_class","pathways").

- custom_bg:

  A named list with character vectors of metabolite names. Default is
  NULL

- termsOfInterest:

  A character containing 'selection' (for default LION-term selection),
  'all', or a vector of term names (see 'pathway').

- condition.x:

  first condition identifier for pairwise comparison.

- condition.y:

  second condition identifier for pairwise comparison.

- ranking.by:

  A character of either 'wilcox.test' (default), 't.test', 'logFC' or
  'BWS', to rank metabolites using the respective statistic. Check
  [`rankScore`](https://alexandrovteam.github.io/S2IsoMEr/reference/rankScore.md)
  for more details. Ignored if `enrichment_type` is 'ORA'.

- gsea.method:

  A character of either 'ks_signed' or 'fgsea'. Ignored if
  `enrichment_type` is 'ORA'.

## Value

An object of class S2IsoMEr.

## References

Korotkevich G, Sukhov V, Sergushichev A (2019). “Fast gene set
enrichment analysis.” bioRxiv. doi:10.1101/060012,
<http://biorxiv.org/content/early/2016/06/20/060012>.

Napoli F (2017). “signed-ks-test.”
<https://github.com/franapoli/signed-ks-test/blob/master/signed-ks-test.R>.

Murakami, Hidetoshi. (2012). Modified Baumgartner Statistics for the
Two-Sample and Multisample Problems: A Numerical Comparison. *Journal of
Statistical Computation and Simulation* 82 (5): 711–28.
[https://doi.org/10.1080/00949655.2010.551516](https://www.tandfonline.com/doi/pdf/10.1080/00949655.2010.551516)

## Examples

``` r
if (FALSE) { # \dontrun{
data("example_ORA_obj", package = "S2IsoMErData")
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
} # }

```
