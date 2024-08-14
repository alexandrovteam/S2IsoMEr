#' Run Simple Metabolite Set Enrichment Analysis (MSEA)
#'
#' @description This function performs simple metabolite set enrichment analysis on a given dataset.
#' It ensures that the ranking conditions are set properly and conducts the enrichment analysis using either the [KS-signed method](https://github.com/franapoli/signed-ks-test/blob/master/signed-ks-test.R) or [fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html).
#'
#' @param object A S2IsoMEr object initialized by \code{\link{initEnrichment}} containing the necessary data and parameters, including annotations, annotation weights, rankings, pathway list, and additional settings.
#' @param min_pathway_size An integer specifying the minimum number of metabolites that must be present in a given term for it to be considered.
#' @return A list containing the results of the enrichment analysis, including:
#'   \itemize{
#'     \item enrichment results
#'     \item the comparison conditions
#'   }
#' @details This function performs a bootstrapped metabolite set enrichment analysis by resampling the annotations and calculating enrichment scores for each bootstrap sample.
#' It handles both KS-signed method](https://github.com/franapoli/signed-ks-test/blob/master/signed-ks-test.R) or [fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html) methods for enrichment analysis and returns a comprehensive summary of the results.
#'
#' @references
#' Korotkevich G, Sukhov V, Sergushichev A (2019). “Fast gene set enrichment analysis.” bioRxiv. doi:10.1101/060012, \url{http://biorxiv.org/content/early/2016/06/20/060012}.
#'
#' Napoli F (2017). “signed-ks-test.” \url{https://github.com/franapoli/signed-ks-test/blob/master/signed-ks-test.R}.
#' @examples
#' \dontrun{
#' # Assuming `my_data` is a properly formatted object of class `S2IsoMEr` initialized by \code{\link{initEnrichment}}
#' result <- Run_bootstrap_MSEA(object = my_data, n_bootstraps = 100)
#' }
#' @export
Run_simple_MSEA = function(object,min_pathway_size = 3){
  object = rankScore.S2IsoMEr(object = object, ranking.by = NULL,
                                alternative = "greater")
  scmat_data = object$scmatrix[object$rankings$rank,]
  bg = object$pathway_list

  if (object$gsea.method == "ks.signed"){
    enrichment_analysis = sapply(names(bg), function(term){

      members_logi <- rownames(scmat_data) %in% bg[[term]]
      if(sum(members_logi)==0){
        data.frame(Term = term,
                   n = 0,
                   ES = NA,
                   p_value = NA)
      } else if (all(members_logi)){
        data.frame(Term = term,
                   n = sum(members_logi),
                   ES = NA,
                   p_value = NA)
      } else {
        ks_results <- ks.test.signed(which(members_logi), which(!members_logi))

        data.frame(Term = term,
                   n = sum(members_logi),
                   ES = ks_results$ES,
                   p_value = ks_results$p.value)
      }

    }, simplify = F) %>% bind_rows()
  }
  else{
    mols = names(object$annotations)
    mols_ranks = object$rankings$statistic
    if (any(is.infinite(mols_ranks))){
      min_rank = min(mols_ranks[!is.infinite(mols_ranks)])
      max_rank = max(mols_ranks[!is.infinite(mols_ranks)])
      mols_ranks[is.infinite(mols_ranks) & mols_ranks < 0] = min_rank - 1
      mols_ranks[is.infinite(mols_ranks) & mols_ranks > 0] = max_rank + 1
    }
    names(mols_ranks) = mols
    mols_ranks = mols_ranks[order(abs(mols_ranks), decreasing = T)]
    mols_ranks = mols_ranks[!duplicated(names(mols_ranks))]
    mols_ranks = mols_ranks[order(mols_ranks)]

    fgsea_res = simple_fgsea(pathways = bg,
                               stats    = mols_ranks,
                               minSize  = min_pathway_size, scoreType = "std")
    fgsea_res = fgsea_res %>% dplyr::select(pathway, size, ES, pval, padj, NES)
    enrichment_analysis = fgsea_res
  }

  colnames(enrichment_analysis)[which(colnames(enrichment_analysis) == "pval")] = "p_value"
  colnames(enrichment_analysis)[which(colnames(enrichment_analysis) == "size")] = "n"
  colnames(enrichment_analysis)[which(colnames(enrichment_analysis) == "pathway")] = "Term"

  enrichment_analysis = enrichment_analysis %>%
    dplyr::filter(Term != "all")

  enrichment_analysis$Term <-
    object$LUT$name[match(enrichment_analysis$Term, object$LUT$ID)]

  object$enrichment_analysis <- list(enrichment_results = enrichment_analysis,
                                     comparison = object$rankings$comparison)

  return(object)
}

#' Run Bootstrapped Metabolite Set Enrichment Analysis (MSEA)
#'
#' @description This function performs bootstrapped metabolite set enrichment analysis on a given dataset. It ensures that the ranking conditions are set properly and conducts the enrichment analysis using either the [KS-signed method](https://github.com/franapoli/signed-ks-test/blob/master/signed-ks-test.R) or [fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html).
#'
#' @param object A S2IsoMEr object initialized by \code{\link{initEnrichment}} containing the necessary data and parameters, including annotations, annotation weights, rankings, pathway list, and additional settings.
#' @param n_bootstraps An integer specifying the number of bootstrap samples to generate.
#' @param min_pathway_size An integer specifying the minimum number of metabolites that must be present in a given term for it to be considered.
#' @param report_ambiguity_scores A logical value indicating whether to calculate and report ambiguity scores for the annotations.
#' @param boot_fract_cutoff A numeric value specifying the minimum fraction of bootstraps in which a pathway must appear to be considered in the final results.
#' @return A list containing the results of the enrichment analysis, including:
#'   \itemize{
#'     \item summarized enrichment results
#'     \item per-bootstrap enrichment results
#'     \item the number of bootstraps
#'     \item the fraction matched to the pathway
#'     \item the comparison conditions
#'     \item annotation ambiguity scores if calculated
#'   }
#' @details This function performs a bootstrapped metabolite set enrichment analysis by resampling the annotations and calculating enrichment scores for each bootstrap sample.
#' It handles both [KS-signed method](https://github.com/franapoli/signed-ks-test/blob/master/signed-ks-test.R) or [fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html) methods for enrichment analysis and returns a comprehensive summary of the results.
#'
#' @references
#' Korotkevich G, Sukhov V, Sergushichev A (2019). “Fast gene set enrichment analysis.” bioRxiv. doi:10.1101/060012, \url{http://biorxiv.org/content/early/2016/06/20/060012}.
#'
#' Napoli F (2017). “signed-ks-test.” \url{https://github.com/franapoli/signed-ks-test/blob/master/signed-ks-test.R}.
#'
#' @examples
#' \dontrun{
#' # Assuming `my_data` is a properly formatted object of class `S2IsoMEr` initialized by \code{\link{initEnrichment}}
#' result <- Run_bootstrap_MSEA(object = my_data, n_bootstraps = 100)
#' }
#' @export
Run_bootstrap_MSEA = function(object,n_bootstraps = 50,
                              min_pathway_size = 3,
                              report_ambiguity_scores = F,
                              boot_fract_cutoff = 0.5){

  if(!all(c(object$condition.x, object$condition.y) ==   object$rankings$comparison)){
    message("condition comparison of the ranking is not the same as set conditions")
    message("run rankScore before enrichment analysis")
    stop(paste0("condition comparison of the ranking is not the same as set conditions\n",
                "run rankScore before enrichment analysis"))
  }

  message("Calculating ranks ... ")

  object = rankScore.S2IsoMEr(object = object, ranking.by = NULL,
                                alternative = "greater")

  if(report_ambiguity_scores){
    ambig_scores = calc_ambiguity(input_iso_list = object$annotations,
                                  weights = object$annotation.weights)
  }
  else{
    ambig_scores = NULL
  }

  message("\nBootstrapping ... \n")

  bootstrapped_sublist <- pbapply::pbsapply(seq(n_bootstraps),       ## bootstrapping
                                            function(n_i) {
                                              sapply(seq(dim(object$scmatrix)[1]), function(row_number_i) {

                                                molecules_to_sample <- object$annotations[[row_number_i]]

                                                if(length(molecules_to_sample)==0){
                                                  mols <- names(object$annotations)[row_number_i]
                                                }

                                                else{
                                                  if(!is.null(object$annotation.weights)){
                                                    weights_to_sample <- object$annotation.weights[[row_number_i]]
                                                  } else {
                                                    weights_to_sample <- rep(x = 1, times = length(molecules_to_sample))
                                                  }


                                                  if(length(molecules_to_sample)!=1){
                                                    mols <- sample(x = molecules_to_sample, size = 1, prob = weights_to_sample)
                                                  } else {          ## length to_sample == 1, R sample() fails with n=1
                                                    mols <- molecules_to_sample
                                                  }
                                                }

                                                names(mols) <- names(object$annotations)[row_number_i]
                                                mols
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

  all_path_metabo = unique(unlist(object$pathway_list))

  fraction_matched_to_pathway <-
    rowMeans(
      pbapply::pbsapply(bootstrapped_sublist, function(i) {
        i %in% all_path_metabo
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
      pbapply::pbsapply(seq(n_bootstraps), function(bootstrap_i){
        sapply(names(pathway_list_slim), function(term){

          members_logi <- bootstrapped_sublist[[bootstrap_i]] %in% pathway_list_slim[[term]]
          if(sum(members_logi)==0){
            data.frame(Term = term,
                       bootstrap = bootstrap_i,
                       n = 0,
                       NES = NA,
                       p_value = NA)
          } else if (all(members_logi)){
            data.frame(Term = term,
                       bootstrap = bootstrap_i,
                       n = sum(members_logi),
                       NES = NA,
                       p_value = NA)
          } else {
            ks_results <- ks.test.signed( which(members_logi), which(!members_logi))

            data.frame(Term = term,
                       bootstrap = bootstrap_i,
                       n = sum(members_logi),
                       NES = ks_results$ES,
                       p_value = ks_results$p.value)
          }

        }, simplify = F) %>% bind_rows()
      },simplify = F) %>% bind_rows()
  }
  else if (object$gsea.method == "fgsea"){
    enrichment_analysis <-
      pbapply::pblapply(seq(n_bootstraps), function(bootstrap_i){

        boot_ranks = object$rankings$statistic
        names(boot_ranks) = gsub("\\+|\\-",".",names(boot_ranks))
        boot_ranks = boot_ranks[rownames(bootstrapped_sublist)]

        boot_mols = bootstrapped_sublist[,bootstrap_i]
        # boot_mols[is.na(boot_mols)] = names(object$rankings$statistic)[sapply(object$annotations,
        #                                                                       length) == 0]

        if (any(is.infinite(boot_ranks))){
          min_rank = min(boot_ranks[!is.infinite(boot_ranks)])
          max_rank = max(boot_ranks[!is.infinite(boot_ranks)])
          boot_ranks[is.infinite(boot_ranks) & boot_ranks < 0] = min_rank - 1
          boot_ranks[is.infinite(boot_ranks) & boot_ranks > 0] = max_rank + 1
        }
        names(boot_ranks) = boot_mols
        boot_ranks = boot_ranks[order(abs(boot_ranks), decreasing = F)]
        boot_ranks = boot_ranks[!duplicated(names(boot_ranks))]
        boot_ranks = boot_ranks[order(boot_ranks)]

        fgsea_res = simple_fgsea(pathways = pathway_list_slim,
                                 stats    = boot_ranks,
                                 minSize  = min_pathway_size,
                                 scoreType = ifelse(any(boot_ranks < 0), "std", "pos"))
        fgsea_res$bootstrap = bootstrap_i
        fgsea_res = fgsea_res %>% dplyr::select(pathway, bootstrap, size, ES, pval, NES)
        fgsea_res
      }) %>% dplyr::bind_rows()
    colnames(enrichment_analysis)[which(colnames(enrichment_analysis) == "pval")] = "p_value"
    colnames(enrichment_analysis)[which(colnames(enrichment_analysis) == "size")] = "n"
    colnames(enrichment_analysis)[which(colnames(enrichment_analysis) == "pathway")] = "Term"
  }

  enrichment_analysis = enrichment_analysis %>%
    dplyr::group_by(Term) %>%
    dplyr::mutate(fraction = length(Term) / length(bootstrapped_sublist)) %>%
    dplyr::ungroup()

  enrichment_analysis$Term <-
    object$LUT$name[match(enrichment_analysis$Term, object$LUT$ID)]

  summarized_enrichment_results <- enrichment_analysis %>%
    dplyr::filter(fraction > boot_fract_cutoff, !is.na(NES)) %>%
    dplyr::group_by(bootstrap) %>%
    dplyr::mutate(q.value = p.adjust(p_value, method = "fdr"))  %>%
    dplyr::group_by(Term) %>%
    dplyr::summarise(n = median(n, na.rm = T),
              ES_median = median(NES, na.rm = T),
              ES_sd = sd(NES, na.rm = T),
              p.value_combined = metap::sumlog(p_value)[["p"]],
              q.value_combined = metap::sumlog(q.value)[["p"]],
              fraction.bootstrap.presence = median(fraction, na.rm = T)) %>%
    dplyr::arrange(q.value_combined)

  object$enrichment_analysis <- list(enrichment_results = summarized_enrichment_results,
                                     per_bootstrap_enrich_results = enrichment_analysis,
                                     n_bootstraps = n_bootstraps,
                                     fraction_matched_to_pathway = fraction_matched_to_pathway,
                                     comparison = object$rankings$comparison,
                                     annotation_ambiguity_scores = ambig_scores)
  return(object)
}
