# Map True Positive Markers to Input Ions from single cell matrix.

This function maps true positive markers to their corresponding ions by
modifying and matching the marker and ion names.

## Usage

``` r
map_TP_markers_to_ions(markers, scm_ions)
```

## Arguments

- markers:

  A character vector containing true positive markers obtained from
  [`get_TP_markers_per_Term`](https://alexandrovteam.github.io/S2IsoMEr/reference/get_TP_markers_per_Term.md)

- scm_ions:

  A character vector containing input ion names.

## Value

A character vector containing the ions that match the given true
positive markers.

## Examples

``` r
if (FALSE) { # \dontrun{
markers <- c("sf.Na", "sf.H")
scm_ions <- c("sf+Na", "sf+H")
map_TP_markers_to_ions(markers, scm_ions)
} # }
```
