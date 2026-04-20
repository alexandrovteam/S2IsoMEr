# Run simple fGSEA enrichment

simple_fgsea calls a simple fGSEA from the fgsea package

## Usage

``` r
simple_fgsea(
  pathways,
  stats,
  minSize = 1,
  maxSize = length(stats) - 1,
  eps = 1e-50,
  scoreType = c("std", "pos", "neg"),
  nPermSimple = 1000
)
```

## Arguments

- pathways:

  List of metabolite sets to check.

- stats:

  Named vector of metabolite-level stats. Names should be the same as in
  'pathways'

- minSize:

  Minimal size of a metabolite set to test. All pathways below the
  threshold are excluded.

- maxSize:

  Maximal size of a metabolite set to test. All pathways above the
  threshold are excluded.

- eps:

  This parameter sets the boundary for calculating the p value.

- scoreType:

  This parameter defines the GSEA score type. Possible options are
  ("std", "pos", "neg")

- nPermSimple:

  Number of permutations in the simple fgsea implementation for
  preliminary estimation of P-values

## Value

A table with GSEA results. Each row corresponds to a tested pathway. The
columns are the following:

- pathway – name of the pathway as in `names(pathway)`;

- pval – an enrichment p-value;

- padj – a BH-adjusted p-value;

- ES – enrichment score, same as in Broad GSEA implementation;

- NES – enrichment score normalized to mean enrichment of random samples
  of the same size;

- nMoreExtreme`-- a number of times a random gene set had a more extreme enrichment score value; \item size -- size of the pathway after removing genes not present in`names(stats)\`.

- leadingEdge – vector with indexes of leading edge genes that drive the
  enrichment, see
  <https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm#_Running_a_Leading>.

## Examples

``` r
bg = list("Term_A" = paste0("mol_", sample(100,20)),"Term_B" = paste0("mol_", sample(100,20)))
mol_ranks = runif(10, -3, 3)
names(mol_ranks) = sample(unique(unlist(bg)), 10)
mol_ranks = mol_ranks[order(mol_ranks)]
simple_fgsea(pathways = bg, stats = mol_ranks, minSize = 3, scoreType = "std")
#>    pathway      pval      padj         ES       NES nMoreExtreme  size
#>     <char>     <num>     <num>      <num>     <num>        <num> <int>
#> 1:  Term_A 0.1153846 0.1236443  0.8000000  1.484242           65     5
#> 2:  Term_B 0.1236443 0.1236443 -0.7664702 -1.448042           56     5
#>                           leadingEdge
#>                                <list>
#> 1: mol_45,mol_55,mol_79,mol_84,mol_73
#> 2:         mol_53,mol_7,mol_96,mol_70
```
