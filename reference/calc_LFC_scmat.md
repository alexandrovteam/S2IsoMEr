# Calculate Log Fold Change (LFC) from Single-Cell Metabolomics Matrix

This method calculates the Log Fold Change (LFC) between two conditions
in a single-cell metabolomics matrix. It compares the mean metabolite
expression levels between the specified conditions and returns the LFC
values for each metabolite.

## Usage

``` r
calc_LFC_scmat(object)
```

## Arguments

- object:

  A `S2IsoMEr` object that contains the single-cell metabolomics matrix
  (`scmatrix`) and metadata about conditions (`conditions`,
  `condition.x`, and `condition.y`).

## Value

A numeric vector of Log Fold Change (LFC) values for each metabolite.
The length of the vector matches the number of metabolites in the
matrix.

## Details

This method computes the mean expression levels of metabolites for two
conditions specified in the `S2IsoMEr` object. The LFC is calculated as
the difference in log2-transformed mean expression levels between the
two conditions.

## Examples

``` r
if (FALSE) { # \dontrun{
data("example_ORA_obj")
lfc_values <- calc_LFC_scmat(example_ORA_obj)
} # }
```
