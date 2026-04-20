# Get intersection between term and query for a Specific Term

This function extracts unique true positive markers for a specified term
from the ORA bootstrap result.

## Usage

``` r
get_TP_markers_per_Term(ORA_boot_df, term_of_interest)
```

## Arguments

- ORA_boot_df:

  A data frame containing the ORA bootstrap results.

- term_of_interest:

  A string specifying the term for which true positive markers should be
  extracted.

## Value

A character vector containing the unique true positive markers
associated with the specified term.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Example usage with hypothetical data
  ORA_boot_df <- data.frame(
    Term = c("Term1", "Term2", "Term1"),
    TP_markers = c("Marker1;Marker2", "Marker3;Marker4", "Marker2;Marker5")
  )
  markers <- get_TP_markers_per_Term(ORA_boot_df, "Term1")
  print(markers)
} # }
```
