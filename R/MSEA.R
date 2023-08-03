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
#' @return
#' @export
#'
#' @examples
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

#' @export
Run_simple_MSEA = function(object, ranking_by = c("t.test", "wilcox.test",
                                                  "BWS", "logFC"),
                           min_pathway_size = 3){
  object = rankScore.bmetenrich(object = object, ranking.by = match.arg(ranking_by))
  scmat_data = object$scmatrix[object$rankings$rank,]
  bg = object$pathway_list

  if (object$gsea.method == "ks.signed"){
    enrichment_analysis = sapply(names(bg), function(term){

      members_logi <- rownames(scmat_data) %in% bg[[term]]
      if(sum(members_logi)==0){
        data.frame(Term = term,
                   n = 0,
                   ES = NA,
                   p.value = NA)
      } else if (all(members_logi)){
        data.frame(Term = term,
                   n = sum(members_logi),
                   ES = NA,
                   p.value = NA)
      } else {
        ks_results <- ks.test.signed(which(members_logi), which(!members_logi))

        data.frame(Term = term,
                   n = sum(members_logi),
                   ES = ks_results$ES,
                   p.value = ks_results$p.value)
      }

    }, simplify = F) %>% bind_rows()
  }
  else{
    mols = rownames(scmat)
    mols_ranks = enrich_obj$rankings$statistic
    names(mols_ranks) = mols
    mols_ranks = mols_ranks[order(mols_ranks)]

    fgsea_res = simple_fgsea(pathways = bg,
                               stats    = mols_ranks,
                               minSize  = min_pathway_size, scoreType = "std")
    enrichment_analysis = fgsea_res
  }
  return(enrichment_analysis)
}
#' @export
Run_bootstrap_MSEA = function(object,min_pathway_size = 3){
  if(!all(c(object$condition.x, object$condition.y) ==   object$rankings$comparison)){
    message("condition comparison of the ranking is not the same as set conditions")
    message("run rankScore before enrichment analysis")
    stop(paste0("condition comparison of the ranking is not the same as set conditions\n",
                "run rankScore before enrichment analysis"))
  }

  cat("\n")
  cat("Bootstrapping...")
  cat("\n")

  bootstrapped_sublist <- pbapply::pbsapply(seq(n),       ## bootstrapping
                                            function(n_i) {
                                              sapply(seq(dim(object$scmatrix)[1]), function(row_number_i) {

                                                molecules_to_sample <- object$annotations[[row_number_i]]

                                                if(length(molecules_to_sample)==0){
                                                  molecules_to_sample <- NA
                                                }

                                                if(!is.null(object$annotation.weights)){
                                                  weights_to_sample <- object$annotation.weights[[row_number_i]]
                                                } else {
                                                  weights_to_sample <- rep(x = 1, times = length(molecules_to_sample))
                                                }


                                                if(length(molecules_to_sample)!=1){
                                                  sample(
                                                    x = molecules_to_sample,
                                                    size = 1,    ## take one
                                                    prob = weights_to_sample)


                                                } else {          ## length to_sample == 1, R sample() fails with n=1
                                                  molecules_to_sample
                                                }

                                              })
                                            }) %>% data.frame


  ## rank
  ## >> test whether ranking is done for set comparison

  bootstrapped_sublist <- bootstrapped_sublist[object$rankings$rank,]

  if(!is.null(object$include)){
    bootstrapped_sublist <- bootstrapped_sublist[object$include,]
  }

  cat("\n")
  cat("Match to pathway...")
  cat("\n")

  fraction_matched_to_pathway <-
    rowMeans(
      pbapply::pbsapply(bootstrapped_sublist, function(i) {
        i %in% object$pathway_list$all
      })
    )


  message(paste0(
    format(
      weighted.mean(fraction_matched_to_pathway) * 100,
      digits = 4,
      nsmall = 1
    ),
    "% of annotations were matched to pathway"
  ))

  ## prune pathway_list to speed up
  molecules_in_dataset <- unique(unlist(bootstrapped_sublist))

  pathway_list_slim <- sapply(object$pathway_list, function(i){
    i[i %in% molecules_in_dataset]
  }, simplify = F)

  pathway_list_slim <- pathway_list_slim[sapply(pathway_list_slim, length) > 0]
  object$pathway_list <- pathway_list_slim

  ## perform enrichment
  cat("\n")
  cat("Perform enrichment analysis...")
  cat("\n")
  if (object$gsea.method == "ks_signed"){
    enrichment_analysis <-
      pbapply::pbsapply(seq(n), function(bootstrap_i){
        sapply(names(pathway_list_slim), function(term){

          members_logi <- bootstrapped_sublist[[bootstrap_i]] %in% pathway_list_slim[[term]]
          if(sum(members_logi)==0){
            data.frame(Term = term,
                       bootstrap = bootstrap_i,
                       n = 0,
                       ES = NA,
                       p.value = NA)
          } else if (all(members_logi)){
            data.frame(Term = term,
                       bootstrap = bootstrap_i,
                       n = sum(members_logi),
                       ES = NA,
                       p.value = NA)
          } else {
            ks_results <- ks.test.signed( which(members_logi), which(!members_logi))

            data.frame(Term = term,
                       bootstrap = bootstrap_i,
                       n = sum(members_logi),
                       ES = ks_results$ES,
                       p.value = ks_results$p.value)
          }

        }, simplify = F) %>% bind_rows()
      },simplify = F) %>% bind_rows()
  }
  else if (object$gsea.method == "fgsea"){
    enrichment_analysis <-
      pbapply::pblapply(seq(n), function(bootstrap_i){
        boot_mols = bootstrap_list[,bootstrap_i]
        boot_mols[is.na(boot_mols)] = names(object$rankings$statistic)[sapply(object$annotations,
                                                                              length) == 0]
        boot_ranks = object$rankings$statistic
        names(boot_ranks) = boot_mols
        boot_ranks = boot_ranks[order(boot_ranks)]

        fgsea_res = simple_fgsea(pathways = pathway_list_slim,
                                 stats    = boot_ranks,
                                 minSize  = min_pathway_size,
                                 scoreType = ifelse(any(boot_ranks < 0), "std", "pos"))
        fgsea_res$bootstrap = boot
      }) %>% dplyr::bind_rows()
    colnames(enrichment_analysis)[which(colnames(enrichment_analysis) == "pval")] = "p.value"
    colnames(enrichment_analysis)[which(colnames(enrichment_analysis) == "size")] = "n"
    colnames(enrichment_analysis)[which(colnames(enrichment_analysis) == "pathway")] = "Term"
  }


  enrichment_analysis$Term <-                    ## match LION name to LION ID
    object$LUT$name[match(enrichment_analysis$Term, object$LUT$ID)]

  object$enrichment_analysis <- list(table = enrichment_analysis,
                                     n_bootstraps = n,
                                     fraction_matched_to_pathway = fraction_matched_to_pathway,
                                     comparison = object$rankings$comparison)
  return(object)
}


ks.test.signed <- function (x, y, ..., alternative = c("two.sided", "less", "greater"), exact = NULL, maxCombSize=10000)
{
  ### https://github.com/franapoli/signed-ks-test/blob/master/signed-ks-test.R

  alternative <- match.arg(alternative)
  DNAME <- deparse(substitute(x))
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 1L)
    stop("not enough 'x' data")
  PVAL <- NULL
  if (is.numeric(y)) {
    DNAME <- paste(DNAME, "and", deparse(substitute(y)))
    y <- y[!is.na(y)]
    n.x <- as.double(n)
    n.y <- length(y)
    if (n.y < 1L)
      stop("not enough 'y' data")
    if (is.null(exact)) {
      exact <- (n.x * n.y < maxCombSize)
      if(!exact)
        warning(paste("P-value not computed exactly because",
                      "of combined sample size"))
    }
    METHOD <- "Two-sample Kolmogorov-Smirnov test"
    TIES <- FALSE
    n <- n.x * n.y/(n.x + n.y)
    w <- c(x, y)
    z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
    if (length(unique(w)) < (n.x + n.y)) {
      if (exact) {
        warning("cannot compute exact p-value with ties")
        exact <- FALSE
      }
      else warning("p-value will be approximate in the presence of ties")
      z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]
      TIES <- TRUE
    }
    STATISTIC <- switch(alternative, two.sided = max(abs(z)),
                        greater = max(z), less = -min(z))

    edge <- which.max(abs(z))
    ES <- z[edge]

    nm_alternative <- switch(alternative, two.sided = "two-sided",
                             less = "the CDF of x lies below that of y", greater = "the CDF of x lies above that of y")
    if (exact && (alternative == "two.sided") && !TIES)
      PVAL <- 1 - .Call(stats:::C_pSmirnov2x, STATISTIC, n.x, n.y)
  }
  else {
    if (is.character(y))
      y <- get(y, mode = "function", envir = parent.frame())
    if (!is.function(y))
      stop("'y' must be numeric or a function or a string naming a valid function")
    METHOD <- "One-sample Kolmogorov-Smirnov test"
    TIES <- FALSE
    if (length(unique(x)) < n) {
      warning("ties should not be present for the Kolmogorov-Smirnov test")
      TIES <- TRUE
    }
    if (is.null(exact))
      exact <- (n < 100) && !TIES
    x <- y(sort(x), ...) - (0:(n - 1))/n
    STATISTIC <- switch(alternative, two.sided = max(c(x,
                                                       1/n - x)), greater = max(1/n - x), less = max(x))
    if (exact) {
      PVAL <- 1 - if (alternative == "two.sided")
        result = tryCatch({
          .C(C_pkolmogorov2x, p = as.double(STATISTIC),
             as.integer(n), PACKAGE = "stats")$p
        }, warning = function(w) {
          warning(w)
        }, error = function(e) {
          .Call(C_pKolmogorov2x, STATISTIC, n)
        }, finally = {
        })

      else {
        pkolmogorov1x <- function(x, n) {
          if (x <= 0)
            return(0)
          if (x >= 1)
            return(1)
          j <- seq.int(from = 0, to = floor(n * (1 -
                                                   x)))
          1 - x * sum(exp(lchoose(n, j) + (n - j) * log(1 -
                                                          x - j/n) + (j - 1) * log(x + j/n)))
        }
        pkolmogorov1x(STATISTIC, n)
      }
    }
    nm_alternative <- switch(alternative, two.sided = "two-sided",
                             less = "the CDF of x lies below the null hypothesis",
                             greater = "the CDF of x lies above the null hypothesis")
  }
  names(STATISTIC) <- switch(alternative, two.sided = "D",
                             greater = "D^+", less = "D^-")
  if (is.null(PVAL)) {
    pkstwo <- function(x, tol = 1e-06) {
      if (is.numeric(x))
        x <- as.double(x)
      else stop("argument 'x' must be numeric")
      p <- rep(0, length(x))
      p[is.na(x)] <- NA
      IND <- which(!is.na(x) & (x > 0))
      if (length(IND))
        p[IND] <- tryCatch({
          tryRes <- .C(stats:::C_pkstwo, length(x[IND]), p = x[IND],
                       as.double(tol), PACKAGE = "stats")$p
        }, warning = function(w) {
          warning(w)
        }, error = function(e) {
          tryRes <- .Call(stats:::C_pKS2, p = x[IND], tol)
        }, finally = {
        })
      p
    }
    PVAL <- ifelse(alternative == "two.sided", 1 - pkstwo(sqrt(n) *
                                                            STATISTIC), exp(-2 * n * STATISTIC^2))
  }
  PVAL <- min(1, max(0, PVAL))
  RVAL <- list(statistic = STATISTIC, p.value = PVAL, alternative = nm_alternative,
               method = METHOD, data.name = DNAME, ES = ES, edge = edge)
  class(RVAL) <- "htest"
  return(RVAL)
}

#' Rank metabolites for bmetenrichr enrichment object
#'
#' rankScore() ranks metabolites of bmetenrichr object to perform bootstrapping metabolite set enrichment analysis
#'
#' @param object A bmetenrichr object.
#' @param ranking.by A character of either 't.test' (default) or 'wilcox.test', to rank metabolites using the respective statistic.
#'
#' @return An object of class bmetenrich.
#'
#' @examples
#' myTestRun <-
#' rankScore(object = myTestRun, ranking.by = 't.test')
#'
#'
#' @export
rankScore <- function (object, ...) {
  UseMethod("rankScore", object)
}

#' @export
rankScore.bmetenrich <- function(object,
                                 ranking.by = NULL,
                                 alternative = c("two.sided", "less", "greater")){

  #IDEA add limma and use p-values as ranking


  if (is.null(ranking.by) & is.null(object$ranking.by)){
    stop("no valid ranking algorithm selected")
  }

  if(!is.null(ranking.by)){
    object$ranking.by <- ranking.by
  }

  if (is.null(object$condition.x) | is.null(object$condition.y)){
    stop("No valid conditions given. Use setConditions().")
  }

  if(object$ranking.by == "t.test"){
    rank_score <-
      apply(object$scmatrix, 1, function(i){
        t.test(x = i[object$conditions == object$condition.x],
               y = i[object$conditions == object$condition.y], alternative = alternative)$statistic
      })


  } else if(object$ranking.by == "wilcox.test"){
    rank_score <-
      apply(object$scmatrix, 1, function(i){
        wilcox.test(x = i[object$conditions == object$condition.x],
                    y = i[object$conditions == object$condition.y], alternative = alternative)$statistic
      })
  } else if (object$ranking.by == "BWS") {
    rank_score <-
      apply(object$scmatrix, 1, function(i){
        BWStest::bws_test(x = i[object$conditions == object$condition.x],
                          y = i[object$conditions == object$condition.y], alternative = alternative)$statistic
      })
  }
  else if (object$ranking.by == "logFC"){
    rank_score <-
      apply(object$scmatrix, 1, function(i){
        x = log2(mean(i[object$conditions == object$condition.x], na.rm = T))
        y = log2(mean(i[object$conditions == object$condition.y], na.rm = T))

        y - x
      })
  }
  else{
    stop("no valid ranking algorithm selected")
  }


  ties <- sum(duplicated(rank_score[object$include]))
  if(ties > 0){
    message(paste0('number of ties: ', ties,
                   " (",
                   format(x = ties / dim(object$scmatrix[object$include,])[1] * 100, digits = 3),
                   "%)"))
  }


  object$rankings <- list(rank = order(rank_score,                   ## first by rankscore
                                       sample(seq_along(rank_score)) ## then by random
  ),
  comparison = c(object$condition.x, object$condition.y),
  statistic = rank_score,
  ranking.by = object$ranking.by)


  return(object)
}

