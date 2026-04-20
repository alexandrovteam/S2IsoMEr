# Load enrichment background

Load_background() retrieves a list of terms and associated metabolites
for enrichment based on requested molecule type, background type and
input feature type.

## Usage

``` r
Load_background(
  mol_type = c("Lipid", "Metabo"),
  bg_type = c("LION", "main_class", "super_class", "sub_class", "pathway"),
  feature_type = c("sf", "name")
)
```

## Arguments

- mol_type:

  Character indicating molecule type. Can be either `Lipid` for Lipids
  and `Metabo` for small molecules

- bg_type:

  Character indicating background type :

  `"LION"`

  :   LION ontology for lipids only.

  `"main_class"`

  :   Level 1 classification - finer classification compared to
      super_class.

  `"super_class"`

  :   Level 0 classification.

  `"sub_class"`

  :   Level 2 classification - finer classification compared to
      main_class.

  `"pathways"`

  :   Biological pathways based on KEGG, Reactome, and SMPDB.

- feature_type:

  Character indicating input feature type. Can be either `sf` for sum
  formula or ion and `name` for molecule name.

## Value

List where each element name correspond to a given term and values are
either metabolite names or sum formulas.

## Examples

``` r
if (FALSE) { # \dontrun{
bg <-
Load_background(mol_type = "Metabo",
                   bg_type = "main_class",
                   feature_type = "name")
} # }

```
