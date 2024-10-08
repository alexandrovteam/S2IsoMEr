% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{example_ORA_obj}
\alias{example_ORA_obj}
\title{Example S2IsoMEr ORA Object}
\format{
An object of class \code{S2IsoMEr} with the following components:
\describe{
\item{scmatrix}{A numeric matrix of metabolite measurements, with rows representing metabolites and columns representing cells or measurements.}
\item{polarization_mode}{Character indicating the polarization mode, either "positive" or "negative".}
\item{enrichment_type}{Character specifying the enrichment analysis type, either "ORA" or "MSEA".}
\item{Annotation_database}{Character or character vector specifying the annotation database(s) used, e.g., "CoreMetabolome", "HMDB", "SwissLipids", or "LipidMaps".}
\item{Custom_database}{Optional custom annotation database provided as a dataframe with two columns: "formula" and "molecule name".}
\item{annotations}{A list of length \code{n}, where each element contains a vector of isomer names associated with each metabolite.}
\item{annotation.weights}{An optional list of length \code{n}, with each element containing a vector of isomer weights.}
\item{isobars_list}{A list of isobaric species, only used when \code{consider_isobars = TRUE}.}
\item{conditions}{A vector of condition identifiers, with length equal to the number of columns in \code{scmatrix}.}
\item{include}{A logical vector indicating which metabolites to include in the analysis.}
\item{consider_isomers}{Logical indicating whether isomers are considered in the analysis.}
\item{consider_isobars}{Logical indicating whether isobars are included in the analysis.}
\item{mass_range_ppm}{Numeric indicating the mass range (in ppm) for isobar identification.}
\item{background_name}{Character specifying the enrichment background type, either "custom" or a combination of molecule and background types.}
\item{endogenous_only}{Logical indicating whether to only consider endogenous metabolites.}
\item{pathway_assoc_only}{Logical indicating whether to include only metabolites associated with a biological pathway.}
\item{remove_expected_predicted}{Logical indicating whether to remove expected or predicted isomers based on HMDB status.}
\item{pathway_list}{A list of pathways associated with the analyzed metabolites.}
\item{LUT}{A lookup table used for annotating metabolites.}
\item{termsSelection}{Character specifying the terms of interest for the enrichment analysis.}
\item{condition.x}{Character specifying the first condition identifier for pairwise comparison.}
\item{condition.y}{Character specifying the second condition identifier for pairwise comparison.}
\item{ranking.by}{Character specifying the ranking method used in MSEA, e.g., "wilcox.test", "t.test", "logFC", or "BWS".}
\item{gsea.method}{Character specifying the GSEA method, either "ks_signed" or "fgsea".}
}
}
\usage{
data(example_ORA_obj)
}
\description{
The \code{example_ORA_obj} dataset is a pre-built S2IsoMEr enrichment object. This object was created using the \code{initEnrichment()} function and contains all the necessary parameters for performing bootstrapping-based ORA.
}
\examples{
# Load the example dataset
data(example_ORA_obj)


}
\keyword{datasets}
