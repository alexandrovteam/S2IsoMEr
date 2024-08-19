#' Example S2IsoMEr ORA Object
#'
#' The `example_ORA_obj` dataset is a pre-built S2IsoMEr enrichment object. This object was created using the `initEnrichment()` function and contains all the necessary parameters for performing bootstrapping-based ORA.
#'
#' @format An object of class `S2IsoMEr` with the following components:
#' \describe{
#'   \item{scmatrix}{A numeric matrix of metabolite measurements, with rows representing metabolites and columns representing cells or measurements.}
#'   \item{polarization_mode}{Character indicating the polarization mode, either "positive" or "negative".}
#'   \item{enrichment_type}{Character specifying the enrichment analysis type, either "ORA" or "MSEA".}
#'   \item{Annotation_database}{Character or character vector specifying the annotation database(s) used, e.g., "CoreMetabolome", "HMDB", "SwissLipids", or "LipidMaps".}
#'   \item{Custom_database}{Optional custom annotation database provided as a dataframe with two columns: "formula" and "molecule name".}
#'   \item{annotations}{A list of length `n`, where each element contains a vector of isomer names associated with each metabolite.}
#'   \item{annotation.weights}{An optional list of length `n`, with each element containing a vector of isomer weights.}
#'   \item{isobars_list}{A list of isobaric species, only used when `consider_isobars = TRUE`.}
#'   \item{conditions}{A vector of condition identifiers, with length equal to the number of columns in `scmatrix`.}
#'   \item{include}{A logical vector indicating which metabolites to include in the analysis.}
#'   \item{consider_isomers}{Logical indicating whether isomers are considered in the analysis.}
#'   \item{consider_isobars}{Logical indicating whether isobars are included in the analysis.}
#'   \item{mass_range_ppm}{Numeric indicating the mass range (in ppm) for isobar identification.}
#'   \item{background_name}{Character specifying the enrichment background type, either "custom" or a combination of molecule and background types.}
#'   \item{endogenous_only}{Logical indicating whether to only consider endogenous metabolites.}
#'   \item{pathway_assoc_only}{Logical indicating whether to include only metabolites associated with a biological pathway.}
#'   \item{remove_expected_predicted}{Logical indicating whether to remove expected or predicted isomers based on HMDB status.}
#'   \item{pathway_list}{A list of pathways associated with the analyzed metabolites.}
#'   \item{LUT}{A lookup table used for annotating metabolites.}
#'   \item{termsSelection}{Character specifying the terms of interest for the enrichment analysis.}
#'   \item{condition.x}{Character specifying the first condition identifier for pairwise comparison.}
#'   \item{condition.y}{Character specifying the second condition identifier for pairwise comparison.}
#'   \item{ranking.by}{Character specifying the ranking method used in MSEA, e.g., "wilcox.test", "t.test", "logFC", or "BWS".}
#'   \item{gsea.method}{Character specifying the GSEA method, either "ks_signed" or "fgsea".}
#' }
#'
#' @usage data(example_ORA_obj)
#'
#' @examples
#' # Load the example dataset
#' data(example_ORA_obj)
#'
#'
#' @keywords datasets
"example_ORA_obj"

#' Example S2IsoMEr MSEA Object
#'
#' The `example_MSEA_obj` dataset is a pre-built S2IsoMEr enrichment object. This object was created using the `initEnrichment()` function and contains all the necessary parameters for performing a bootstrapping-based MSEA.
#'
#' @format An object of class `S2IsoMEr` with the following components:
#' \describe{
#'   \item{scmatrix}{A numeric matrix of metabolite measurements, with rows representing metabolites and columns representing cells or measurements.}
#'   \item{polarization_mode}{Character indicating the polarization mode, either "positive" or "negative".}
#'   \item{enrichment_type}{Character specifying the enrichment analysis type, either "ORA" or "MSEA".}
#'   \item{Annotation_database}{Character or character vector specifying the annotation database(s) used, e.g., "CoreMetabolome", "HMDB", "SwissLipids", or "LipidMaps".}
#'   \item{Custom_database}{Optional custom annotation database provided as a dataframe with two columns: "formula" and "molecule name".}
#'   \item{annotations}{A list of length `n`, where each element contains a vector of isomer names associated with each metabolite.}
#'   \item{annotation.weights}{An optional list of length `n`, with each element containing a vector of isomer weights.}
#'   \item{isobars_list}{A list of isobaric species, only used when `consider_isobars = TRUE`.}
#'   \item{conditions}{A vector of condition identifiers, with length equal to the number of columns in `scmatrix`.}
#'   \item{include}{A logical vector indicating which metabolites to include in the analysis.}
#'   \item{consider_isomers}{Logical indicating whether isomers are considered in the analysis.}
#'   \item{consider_isobars}{Logical indicating whether isobars are included in the analysis.}
#'   \item{mass_range_ppm}{Numeric indicating the mass range (in ppm) for isobar identification.}
#'   \item{background_name}{Character specifying the enrichment background type, either "custom" or a combination of molecule and background types.}
#'   \item{endogenous_only}{Logical indicating whether to only consider endogenous metabolites.}
#'   \item{pathway_assoc_only}{Logical indicating whether to include only metabolites associated with a biological pathway.}
#'   \item{remove_expected_predicted}{Logical indicating whether to remove expected or predicted isomers based on HMDB status.}
#'   \item{pathway_list}{A list of pathways associated with the analyzed metabolites.}
#'   \item{LUT}{A lookup table used for annotating metabolites.}
#'   \item{termsSelection}{Character specifying the terms of interest for the enrichment analysis.}
#'   \item{condition.x}{Character specifying the first condition identifier for pairwise comparison.}
#'   \item{condition.y}{Character specifying the second condition identifier for pairwise comparison.}
#'   \item{ranking.by}{Character specifying the ranking method used in MSEA, e.g., "wilcox.test", "t.test", "logFC", or "BWS".}
#'   \item{gsea.method}{Character specifying the GSEA method, either "ks_signed" or "fgsea".}
#' }
#'
#' @usage data(example_MSEA_obj)
#'
#' @examples
#' # Load the example dataset
#' data(example_MSEA_obj)
#'
#'
#' @keywords datasets
"example_MSEA_obj"

#' Example custom universe for ORA
#'
#' `example_ORA_custom_universe` is a character vector of ions representing a custom universe for ORA
#'
#' @format A character vector with 1618 elements, each representing an ion (sum formula + adduct)
#'
#' @usage data(example_ORA_custom_universe)
#'
#' @examples
#' # Load the example vector
#' data(example_ORA_custom_universe)
#'
#' # View the first few elements
#' head(example_ORA_custom_universe)
"example_ORA_custom_universe"

#' Example Multi condition MSEA results
#'
#' The `example_MSEA_multicond` dataset is a data frame that contains results of running \code{\link{Run_bootstrap_MSEA}} on multiple pairwise conditions compared to same reference.
#'
#' @format A data frame with 100 rows and 9 variables:
#' \describe{
#'   \item{Term}{Character: Term name of metabolite set (e.g. Triacylglycerols)}
#'   \item{n}{Numeric: Median term/query overlap size over bootstraps.}
#'   \item{ES_median}{Numeric: Median of NES (normalized enrichment score) for fgsea or ES for ks_signed}
#'   \item{ES_sd}{Numeric: stadndard deviation of NES (normalized enrichment score) for fgsea or ES for ks_signed}
#'   \item{p.value_combined}{Numeric: Combined p-value using `metap::sumlog`}
#'   \item{q.value_combined}{Numeric: Combined adjusted p-value using `metap::sumlog`}
#'   \item{fraction.bootstrap.presence}{Numeric: Proportion of bootstraps the given term was tested}
#'   \item{condition.x}{Character: Reference condition}
#'   \item{condition.y}{Character: Query condition}
#' }
#'
#' @usage data(example_MSEA_multicond)
#'
#' @examples
#' # Load the example data frame
#' data(example_MSEA_multicond)
#'
#' # View the first few rows
#' head(example_MSEA_multicond)
#'
#'
#' @keywords datasets
"example_MSEA_multicond"

#' Example List of Markers for ORA
#'
#' `example_ORA_markers` is a list that provides markers for use in ORA
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{upregulated}{A character vector of upregulated ions}
#'   \item{downregulated}{A character vector of downregulated ions}
#'   \item{all}{A character vector of all markers, including both upregulated and downregulated}
#' }
#'
#' @usage data(example_ORA_markers)
#'
#' @examples
#' # Load the marker list
#' data(example_ORA_markers)
#'
#' @keywords datasets
"example_ORA_markers"
