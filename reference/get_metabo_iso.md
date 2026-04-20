# Get Metabolite Isomers and Isobars

This function retrieves metabolite isomers and isobars for a given list
of sum formula, considering various parameters such as adducts,
polarization mode, and mass range.

## Usage

``` r
get_metabo_iso(
  sf_vec,
  consider_isobars = T,
  polarization_mode = NA,
  mass_range_ppm = 3,
  annot_db = "HMDB",
  annot_custom_db = NULL,
  use_LION = F,
  endogenous_only = T,
  pathway_assoc_only = F,
  remove_expected_predicted = T
)
```

## Arguments

- sf_vec:

  A character vector of metabolite formulas with optional adducts. The
  formulas should be formatted as `formula[+adduct]` or
  `formula[-adduct]`.

- consider_isobars:

  A logical indicating whether to consider isobars. If TRUE, the
  function will identify and include isobars based on the provided mass
  range.

- polarization_mode:

  A character string specifying the ionization mode. Options are
  "positive" or "negative". Only relevant if `consider_isobars` is TRUE.

- mass_range_ppm:

  A numeric value specifying the mass range in parts per million (ppm)
  for identifying isobars.

- annot_db:

  A character vector specifying the annotation databases to use. Options
  include "HMDB", "LipidMaps", "SwissLipids", and "CoreMetabolome".

- annot_custom_db:

  A custom annotation database. If not NULL, this custom database will
  be used instead of the specified annotation databases.

- use_LION:

  A logical indicating whether to use the LION version of the LipidMaps
  database if "LipidMaps" is included in `annot_db`.

- endogenous_only:

  A logical indicating whether to include only endogenous metabolites.
  Applies only to "HMDB" and "CoreMetabolome".

- pathway_assoc_only:

  A logical indicating whether to include only metabolites associated
  with pathways. Applies only to "HMDB" and "CoreMetabolome".

- remove_expected_predicted:

  A logical indicating whether to remove metabolites with a status of
  "expected" or "predicted" in HMDB and CoreMetabolome databases.

## Value

A named list where each element corresponds to a formula from `sf_vec`.
Each element is a character vector containing metabolite names matching
the formula and its isobars, if `consider_isobars` is TRUE.

## Details

This function processes a vector of metabolite formulas and retrieves
the corresponding metabolite names from the specified annotation
databases. If `consider_isobars` is TRUE, it will also find isobars
based on the specified mass range and polarization mode.

## Examples

``` r
if (FALSE) { # \dontrun{
sf_vec <- c("C6H12O6+H", "C6H12O6-H")
results <- get_metabo_iso(
  sf_vec,
  consider_isobars = TRUE,
  polarization_mode = "positive",
  mass_range_ppm = 5
)
} # }
```
