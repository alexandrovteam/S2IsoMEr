# Rank metabolites for S2IsoMEr enrichment object

rankScore() ranks metabolites of S2IsoMEr object to perform
bootstrapping metabolite set enrichment analysis

## Usage

``` r
rankScore(object, ranking.by = "wilcox.test", alternative = "greater")
```

## Arguments

- object:

  A S2IsoMEr object.

- ranking.by:

  A character string specifying the ranking method. Options include:

  - `"wilcox.test"` (default) - Ranks metabolites based on the results
    of a Wilcoxon test.

  - `"t.test"` - Ranks metabolites based on the results of a t-test.

  - `"logFC"` - Ranks metabolites based on log fold change.

  - `"BWS"` - Ranks metabolites based on BWS
    (Baumgartner-Weiss-Schindler ).

- alternative:

  A character string specifying the alternative hypothesis for
  statistical tests. Options are:

  - `"two.sided"` - Test for differences in both directions.

  - `"less"` - Test if the second condition is less than the first

  - `"greater"` - Test if the second condition is greater than the first

  Default is `"greater"`.

## Value

An object of class S2IsoMEr.

## References

Murakami, Hidetoshi. (2012). Modified Baumgartner Statistics for the
Two-Sample and Multisample Problems: A Numerical Comparison. *Journal of
Statistical Computation and Simulation* 82 (5): 711–28.
[https://doi.org/10.1080/00949655.2010.551516](https://www.tandfonline.com/doi/pdf/10.1080/00949655.2010.551516)

## Examples

``` r
if (FALSE) { # \dontrun{
data("example_ORA_obj")
example_ORA_obj <-
rankScore(object = example_ORA_obj, ranking.by = 't.test')
} # }

```
