# Set conditions for enrichment analysis

Set conditions for enrichment analysis

## Usage

``` r
setConditions(object, condition.x = NULL, condition.y = NULL)
```

## Arguments

- object:

  A S2IsoMEr object.

- condition.x:

  A optional character describing the reference condition.

- condition.y:

  A optional character describing condition to interest.

## Value

An object of class `S2IsoMEr`.

## Examples

``` r
if (FALSE) { # \dontrun{
data("example_MSEA_obj", package = "S2IsoMErData")
example_MSEA_obj = setConditions(example_MSEA_obj, condition.x = 'U', condition.y = "FI")
} # }
```
