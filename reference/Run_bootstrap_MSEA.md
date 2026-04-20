# Run Bootstrapped Metabolite Set Enrichment Analysis (MSEA)

This function performs bootstrapped metabolite set enrichment analysis
on a given dataset. It ensures that the ranking conditions are set
properly and conducts the enrichment analysis using either the
[KS-signed
method](https://github.com/franapoli/signed-ks-test/blob/master/signed-ks-test.R)
or
[fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html).

## Usage

``` r
Run_bootstrap_MSEA(
  object,
  n_bootstraps = 50,
  min_pathway_size = 3,
  report_ambiguity_scores = F,
  boot_fract_cutoff = 0.5
)
```

## Arguments

- object:

  A S2IsoMEr object initialized by
  [`initEnrichment`](https://alexandrovteam.github.io/S2IsoMEr/reference/initEnrichment.md)
  containing the necessary data and parameters, including annotations,
  annotation weights, rankings, pathway list, and additional settings.

- n_bootstraps:

  An integer specifying the number of bootstrap samples to generate.

- min_pathway_size:

  An integer specifying the minimum number of metabolites that must be
  present in a given term for it to be considered.

- report_ambiguity_scores:

  A logical value indicating whether to calculate and report ambiguity
  scores for the annotations.

- boot_fract_cutoff:

  A numeric value specifying the minimum fraction of bootstraps in which
  a pathway must appear to be considered in the final results.

## Value

A list containing the results of the enrichment analysis, including:

- summarized enrichment results

- per-bootstrap enrichment results

- the number of bootstraps

- the fraction matched to the pathway

- the comparison conditions

- annotation ambiguity scores if calculated

## Details

This function performs a bootstrapped metabolite set enrichment analysis
by resampling the annotations and calculating enrichment scores for each
bootstrap sample. It handles both [KS-signed
method](https://github.com/franapoli/signed-ks-test/blob/master/signed-ks-test.R)
or
[fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html)
methods for enrichment analysis and returns a comprehensive summary of
the results.

## References

Korotkevich G, Sukhov V, Sergushichev A (2019). “Fast gene set
enrichment analysis.” bioRxiv. doi:10.1101/060012,
<http://biorxiv.org/content/early/2016/06/20/060012>.

Napoli F (2017). “signed-ks-test.”
<https://github.com/franapoli/signed-ks-test/blob/master/signed-ks-test.R>.

## Examples

``` r
if (FALSE) { # \dontrun{
data("example_MSEA_obj", package = "S2IsoMErData")
result <- Run_bootstrap_MSEA(object = example_MSEA_obj, n_bootstraps = 50)
} # }
```
