# bmetenrichr

## Overview

`Bmetenrichr` is an R package for metabolite enrichment analysis in imaging MS. `Bmetenrichr` can be applied to spatial and single-cell metabolomics datasets and addresses the challenge of metabolite identification ambiguity. The key idea to handle molecular isomers and/or isobars is to propagate the molecular ambiguity to the enrichment results as follows: We apply iterative random sampling (bootstrapping), perform enrichment analysis for each iteration, and report summarized results.

## Installation

```{r, eval=FALSE}
# If you don't have devtools installed yet, you can install it via:
install.packages("devtools")

# Install the development version from GitHub
devtools::install_github("Bisho2122/bmetenrichr")
```

## Features

`Bmetenrichr` supports both overrepresentation analysis (ORA) and metabolite set enrichment analysis (MSEA) for metabolite and lipid-based backgrounds. Metabolite backgrounds are curated from [RAMP-DB 2.0](https://academic.oup.com/bioinformatics/article/39/1/btac726/6827287), encompassing biological pathways and metabolic classes, which are further categorized into super, main, and sub-classes. The pathways from RAMP-DB integrate multiple resources, including SMPDB, Reactome, KEGG, and WikiPathways. For lipids, similar background types are provided, with the addition of the [LION](http://lipidontology.com/) lipid ontology. Except for the LION ontology, each term is mapped to either molecule names or sum formulas, allowing enrichment analysis to be performed with or without consideration of isomeric/isobaric ambiguity.

## Licence

`Bmetenrichr` is licensed under the MIT License. See the [LICENSE](https://github.com/bisho2122/bmetenrichr?tab=MIT-2-ov-file) file for more details.

## Citation

bioRxiv paper coming soon.
