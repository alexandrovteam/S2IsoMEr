# Bar Plot for Simple ORA Results

This function creates a bar plot for simple (no bootstraps)
Over-Representation Analysis (ORA) results, indicating the significance
of terms based on a q-value cutoff.

## Usage

``` r
barplot_ORA_simple(ORA_simple_res, q_val_cutoff = 0.05)
```

## Arguments

- ORA_simple_res:

  A data frame containing the ORA results.

- q_val_cutoff:

  A numeric value specifying the q-value cutoff for significance.
  Default is 0.05.

## Value

A ggplot object displaying the bar plot of ORA results.

## Details

The function creates a bar plot where the terms are reordered by their
scores in descending order. Bars are colored based on whether the term's
score is above or below the negative log10 of the q-value cutoff. A
dashed line indicates the cutoff for significance.

## Examples

``` r
if (FALSE) { # \dontrun{
data("example_ORA_markers", package = "S2IsoMErData")
bg = Load_background(mol_type = "Metabo",bg_type = "main_class",feature_type = "sf")
enrich_res = Run_simple_ORA(marker_list = example_ORA_markers,background = bg)
p = barplot_ORA_simple(enrich_res, q_val_cutoff = 0.2)
} # }

```
