# Combine p-values Using the Sum of Logs Method

This function combines multiple p-values using Fisher's method (sum of
logarithms of p-values), which is often used in meta-analysis. It
computes a chi-squared statistic and the corresponding p-value, based on
the assumption that the input p-values are independent. This function is
adapted from the `sumlog` function in the `metap` package.

## Usage

``` r
metap_sumlog_pvals(p, log.p = FALSE)
```

## Arguments

- p:

  A numeric vector of p-values. Each p-value must lie in the range (0,
  1\]. Invalid values (e.g., \<= 0 or \> 1) are excluded.

- log.p:

  Logical. If `TRUE`, the returned p-value is given as `log(p)`.

## Value

A list with the following components:

- `chisq`:

  The chi-squared statistic computed as `-2 * sum(log(p))`.

- `df`:

  The degrees of freedom, equal to `2 * length(validp)`.

- `p`:

  The combined p-value corresponding to the chi-squared statistic.

- `validp`:

  A vector of valid p-values after filtering out invalid ones.

## Details

This function applies Fisher's method for combining p-values. It first
filters out invalid p-values (values \<= 0 or \> 1). If fewer than two
valid p-values remain, a warning is issued, and `NA` values are
returned. If some p-values are filtered out, the function issues a
warning indicating that some studies were omitted.

Fisher's method is based on the chi-squared distribution, with degrees
of freedom proportional to the number of valid p-values (2 \* number of
valid p-values).

For further details on the original implementation, please refer to the
`sumlog` function in the `metap` package:
<https://cran.r-project.org/web/packages/metap/metap.pdf>

## References

`sumlog` function from the `metap` package:
<https://cran.r-project.org/web/packages/metap/metap.pdf>

## Examples

``` r
# Example with valid p-values
pvals <- c(0.01, 0.03, 0.05, 0.2)
result <- metap_sumlog_pvals(pvals)
print(result)
#> $chisq
#> [1] 25.4338
#> 
#> $df
#> [1] 8
#> 
#> $p
#> [1] 0.001312015
#> 
#> $validp
#> [1] 0.01 0.03 0.05 0.20
#> 

# Example with some invalid p-values
pvals <- c(0.01, 0, 1.2, 0.05)
result <- metap_sumlog_pvals(pvals)
#> Warning: Some studies omitted
print(result)
#> $chisq
#> [1] 15.2018
#> 
#> $df
#> [1] 4
#> 
#> $p
#> [1] 0.004300451
#> 
#> $validp
#> [1] 0.01 0.05
#> 
```
