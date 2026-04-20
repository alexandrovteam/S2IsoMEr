# Build Isomer/Isobar Background

This function constructs an isomeric/isobaric background database based
on specified annotation databases and additional filtering criteria.

## Usage

``` r
build_iso_bg(
  annot_db = "HMDB",
  annot_custom_db = NULL,
  use_LION = F,
  endogenous_only = T,
  pathway_assoc_only = F,
  remove_expected_predicted = T
)
```

## Arguments

- annot_db:

  A character vector specifying the annotation databases to use. Options
  include "HMDB", "LipidMaps", "SwissLipids", and "CoreMetabolome".

- annot_custom_db:

  A custom annotation database. If not NULL, this custom database will
  be used instead of the specified annotation databases.

- use_LION:

  A logical indicating whether to use the LION ontology or internal
  LipidMaps backgrounds if "LipidMaps" is included in `annot_db`.

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

A data frame containing the isomers and isobers background based on the
specified criteria.

## Details

This function builds an isomer/isobar background by either using the
specified annotation databases or a provided custom annotation database.
Various filtering options are available to refine the background.

## Examples

``` r
if (FALSE) { # \dontrun{
iso_bg <- build_iso_bg(annot_db = "HMDB", use_LION = FALSE, endogenous_only = TRUE,
                       pathway_assoc_only = FALSE, remove_expected_predicted = TRUE)
custom_db <- data.frame(
  ID = 1:3,
  Name = c("Met1", "Met2", "Met3"),
  Endogenous = "Yes",
  Pathway_assoc = "Yes",
  HMDB_status = "known"
)
iso_bg_custom <- build_iso_bg(annot_custom_db = custom_db)
} # }
```
