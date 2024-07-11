#' Run simple fGSEA enrichment
#'
#' simple_fgsea calls a simple fGSEA from the fgsea package
#'
#' @param pathways List of metabolite sets to check.
#' @param stats Named vector of metabolite-level stats. Names should be the same as in 'pathways'
#' @param minSize Minimal size of a metabolite set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a metabolite set to test. All pathways above the threshold are excluded.
#' @param eps This parameter sets the boundary for calculating the p value.
#' @param scoreType This parameter defines the GSEA score type. Possible options are ("std", "pos", "neg")
#' @param nPermSimple Number of permutations in the simple fgsea implementation for preliminary estimation of P-values
#'
#' @return A table with GSEA results. Each row corresponds to a tested pathway.
#' The columns are the following:
#' \itemize{
#'  \item pathway -- name of the pathway as in `names(pathway)`;
#'  \item pval -- an enrichment p-value;
#'  \item padj -- a BH-adjusted p-value;
#'  \item ES -- enrichment score, same as in Broad GSEA implementation;
#'  \item NES -- enrichment score normalized to mean enrichment of random samples of the same size;
#'  \item nMoreExtreme` -- a number of times a random gene set had a more
#'      extreme enrichment score value;
#'  \item size -- size of the pathway after removing genes not present in `names(stats)`.
#'  \item leadingEdge -- vector with indexes of leading edge genes that drive the enrichment, see \url{http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Running_a_Leading}.
#' }
#'
#' @export
#'
#' @examples
#' bg = list("Term_A" = paste0("mol_", sample(100,20)),"Term_B" = paste0("mol_", sample(100,20)))
#' mol_ranks = runif(10, -3, 3)
#' names(mol_ranks) = sample(unique(unlist(bg)), 10)
#' mol_ranks = mol_ranks[order(mol_ranks)]
#' simple_fgsea(pathways = bg, stats = mol_ranks, minSize = 3, scoreType = "std")
simple_fgsea = function(pathways,
                        stats,
                        minSize = 1,
                        maxSize = length(stats)-1,
                        eps = 1e-50,
                        scoreType   = c("std", "pos", "neg"),
                        nPermSimple = 1000){
  scoreType <- match.arg(scoreType)
  pp <- fgsea:::preparePathwaysAndStats(pathways, stats, minSize, maxSize, 1, scoreType)
  pathwaysFiltered <- pp$filtered
  pathwaysSizes <- pp$sizes
  stats <- pp$stats
  m <- length(pathwaysFiltered)
  minSize <- max(minSize, 1)
  eps <- max(0, min(1, eps))
  gseaStatRes <- do.call(rbind,
                         lapply(pathwaysFiltered,
                                fgsea:::calcGseaStat,
                                stats             = stats,
                                returnLeadingEdge = TRUE,
                                scoreType         = scoreType))
  leadingEdges <- mapply("[", list(names(stats)), gseaStatRes[, "leadingEdge"], SIMPLIFY = FALSE)
  pathwayScores <- unlist(gseaStatRes[, "res"])
  seeds <- sample.int(10^9, 1)
  simpleFgseaRes <- fgsea:::fgseaSimpleImpl(pathwayScores=pathwayScores, pathwaysSizes=pathwaysSizes,
                                            pathwaysFiltered=pathwaysFiltered, leadingEdges=leadingEdges,
                                            permPerProc=nPermSimple, seeds=seeds, toKeepLength=m,
                                            stats=stats, BPPARAM= BiocParallel::SerialParam(), scoreType=scoreType)

  return(simpleFgseaRes)
}
