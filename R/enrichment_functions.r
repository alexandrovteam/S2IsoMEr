
###

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

#' Generate bmetenrichr enrichment object
#'
#' initEnrichment() creates object to perform bootstrapping metabolite set enrichment analysis
#'
#' @param scmatrix A numeric matrix of n metabolites (rows) and m cells or measurments (columns).
#' @param annotations Either (i) a list of length n, with each element contains a vector of isomer names,
#' or (ii) a vector of length n containing molecular formulas with ("C44H84NO8P.H") or without adduct ("C44H84NO8P").
#' In the second case, bmetenrichr uses the CoreMetabolome, LIPIDMAPS, SwissLipids, and HMDB databases from METASPACE (https://metaspace2020.eu/) to generate an annotation list automatically.
#' @param annotation.weights An optional list of length n, each element contains a vector of isomer weights. Only when annotations is provided as list.
#' @param isobars A logical indicating whether to include isobars and isomers (default = FALSE), FALSE will only include isomers.
#' Will be neglected when the annotations are provided as list.
#' @param mass_range_ppm A numeric indicating the mass range in ppm (default: mass_range_ppm = 3). Molecular formulas + adducts within this range will be treated as isobars. Only required when isobars = TRUE.
#' @param polarization_mode A character with either 'positive' (default) or 'negative'. Only required when isobars = TRUE. When set to 'positive', included adducts are '+H', '+Na', and '+K'.
#' When set to 'negative', included adducts are '-H', '+Cl'.
#' @param conditions A vector of length m with condition identifiers.
#' @param include An optional logical vector of length n indicating whether to include the annotations in the analysis.
#' @param condition.x first condition identifier for pairwise comparison.
#' @param condition.y second condition identifier for pairwise comparison.
#' @param pathway A named list with character vectors of metabolite names. When default 'LION' is used, bmetenrichr uses the preset LION metabolite set.
#' @param termsOfInterest A character containing 'selection' (for default LION-term selection), 'all', or a vector of term names (see 'pathway').
#' @param ranking.by A character of either 't.test' or 'wilcox.test', to rank metabolites for the respective statistic.
#'
#' @return An object of class bmetenrich.
#' @examples
#' myTestRun <-
#' initEnrichment(scmatrix = scMatrix,
#'                    annotations = my_annotations,
#'                    annotation.weights = my_weights,
#'                    conditions = my_conditions,
#'                    condition.x = "A",
#'                    condition.y = "B" )
#'
#'
#' @export
initEnrichment <- function(scmatrix,
                               annotations,
                               annotation.weights = NULL,
                               isobars = FALSE,
                               mass_range_ppm = 3,
                               polarization_mode = "positive",
                               conditions,
                               include = NULL,
                               pathway = "LION",
                               termsOfInterest = "selection",
                               condition.x = NULL,
                               condition.y = NULL,
                               ranking.by = "t.test"){


  if (dim(scmatrix)[1] !=  length(annotations)){
    stop("single-cell matrix and annotations do not have the same length")
  }
  if (!is.null(include) & length(annotations) != length(include)){
    stop("annotations and include do not have the same length")
  }
  if (dim(scmatrix)[2] !=  length(conditions)){
    stop("single-cell matrix and conditions do not have the same length")
  }
  if (isobars & is.list(annotations)){
    stop("isobars = TRUE in combination with a custom list of annotations is not supported")
  }
  if(isobars){
    if(!polarization_mode %in% c('positive','negative')){
      stop("polarization_mode should be either 'positive' or 'negative'")
    } else {
      message(paste0("polarization_mode is: ", polarization_mode))
    }
  }


  if (!is.null(annotation.weights)) {
    ## are the weights in the right format?

    if (length(annotation.weights) != length(annotations)) {
      stop("annotations and annotation.weights do not have the same length")
    } else {
      ## extra check
      if (!all(sapply(annotation.weights, length) == sapply(annotations, length))) {
        stop("annotations and annotation.weights do not have the same length ")
      }
    }
  }

  if (!is.null(condition.x) ){
    if(!any(condition.x == unique(conditions))){
      stop("submitted condition.x not found in dataset")
    }
  }

  if (!is.null(condition.y) ){
    if(!any(condition.y == unique(conditions))){
      stop("submitted condition.y not found in dataset")
    }
  }


  if (length(pathway) == 1) {
    if (pathway == "LION") {
      pathway_list <- pathway_list_LION

    } else {
      stop("pathway input is not in the right format")
      }
  } else {
    ## pathway is longer than length 1
    if (is.list(pathway)) {
      ## pathway is a list
      pathway_list <- pathway
      pathway_list$all <- unique(unlist(pathway_list))
      ## make an *ad-hoc* LUT
      LUT <- data.frame(ID = names(pathway_list),
                        name = names(pathway_list))
    } else {
        stop("pathway input is not in the right format")
      }

  }


  if (is.list(annotations)){
    ## use provided list when used as input
    annotation_list <- annotations
  } else if (is.character(annotations)){

    ## when vector of moceular formulas is provided, generate list with molecular names

    ## change + or - adduct into .
    annotation_formulas_adduct <- gsub("\\+|\\-",".",annotations)

    annotation_formulas <- gsub("\\..+$","",annotation_formulas_adduct)   ## remove adduct

    annotation_adduct <- gsub("^.+\\.","",annotation_formulas_adduct)

    cat("\nParsing isomers...\n")

    annotation_list <-
    lapply(annotation_formulas, function(annotation_formula_i){
      metaspace_databases$name[metaspace_databases$formula == annotation_formula_i]
    })


  } else {
    stop("annotations not in the right format")
  }

   if(isobars){
     #TODO change the IF condition above to accomodate metabolites not just lion

     cat("Parsing potential isobars...\n\n")

     switch(polarization_mode,
            'positive' = {
              col_name <- paste0("pos",annotation_adduct)

              exact_masses_slim <- exact_masses[,c(1,which(grepl("pos",colnames(exact_masses))))]
              colnames(exact_masses_slim) <- gsub("^pos","",colnames(exact_masses_slim))
              exact_masses_slim <- exact_masses_slim %>% tidyr::pivot_longer(cols = -1, values_to = "mass", names_to = "adduct")
              exact_masses_slim$formula_adduct <- paste0(exact_masses_slim$formula,".",exact_masses_slim$adduct)
              },
            'negative' = {
              col_name <- paste0("neg",annotation_adduct)

              exact_masses_slim <- exact_masses[,c(1,which(grepl("neg",colnames(exact_masses))))]
              colnames(exact_masses_slim) <- gsub("^neg","",colnames(exact_masses_slim))
              exact_masses_slim <- exact_masses_slim %>% tidyr::pivot_longer(cols = -1, values_to = "mass", names_to = "adduct")
              exact_masses_slim$formula_adduct <- paste0(exact_masses_slim$formula,".",exact_masses_slim$adduct)
              })

     exact_masses_annotations <-
       mapply(annotation_formulas_i = annotation_formulas,
              col_name_i = col_name,
              function(annotation_formulas_i, col_name_i) {
                mass <-
                  exact_masses[exact_masses$formula == annotation_formulas_i, col_name_i]
                if (length(mass) == 0) {
                  mass = NA
                }
                mass

              })

     names(exact_masses_annotations) <- annotation_formulas_adduct

     isobars_list <-
     sapply(exact_masses_annotations, function(mass_i){
       if(is.na(mass_i)){
         character(0)
       } else {
         exact_masses_slim$formula_adduct[
           between(exact_masses_slim$mass,
                   left = mass_i - (mass_range_ppm * mass_i / 1e6),
                   right = mass_i + (mass_range_ppm * mass_i / 1e6))]
       }
     }, simplify = F)

     isobars_list[sapply(isobars_list, length) > 1]

     ## remove self isobars
     isobars_list <-
     sapply(names(isobars_list), function(i){
       isobars_list[[i]][!isobars_list[[i]] %in% i]
     }, simplify = F)


     annotation_list <-
       lapply(seq_along(annotation_formulas), function(i){
         c(isomer = annotation_list[[i]],
           isobar = metaspace_databases$name[metaspace_databases$formula %in% gsub("\\..+$","",isobars_list[[i]])])

       })

   } else {   ## isbars == FALSE
     isobars_list <- NULL
   }


  if(termsOfInterest == "selection" & pathway[[1]][1] == "LION"){

    ## filter pathway_list by terms of interest
    # termsSelection <- read.csv(file = 'data-raw/LION_selection.csv')
    pathway_list <- pathway_list[names(pathway_list) %in% c("all", termsSelection$LION_ID)]

  } else if(termsOfInterest == "all"){
    termsSelection <- names(pathway_list)
    ## no filtering required

  } else if(termsOfInterest == "selection" & pathway[[1]][1] != "LION"){

    stop("Setting termsOfInterest = 'selection' is only valid for pathway = 'LION'\nUse termsOfInterest = 'all' instead")
  } else {

    pathway_list <- pathway_list[termsOfInterest]
    termsSelection <- names(pathway_list)

  }


  object <-
  structure(
  list(scmatrix = scmatrix,
       annotations = annotation_list,
       annotation.weights = annotation.weights,
       isobars_list = isobars_list,
       conditions = conditions,
       include = include,
       pathway = if(length(pathway) == 1){"LION"}else{"custom"},
       pathway_list = pathway_list,
       LUT = LUT,
       termsSelection = termsSelection,
       condition.x = condition.x,
       condition.y = condition.y,
       ranking.by = ranking.by),
  class = "bmetenrich")

  print(object)
  return(object)
}

#' @export
print.bmetenrich <- function(object){
  cat("single-cell metabolomics matrix of", dim(object$scmatrix)[1], "metabolites and",
      dim(object$scmatrix)[2], "cells\n")
  cat("active pathway:", object$pathway,"\n\n")

  cat("conditions:", paste(unique(object$conditions), collapse = ", "),"\n\n")

  cat("condition.x:", object$condition.x,"\n")
  cat("condition.y:", object$condition.y,"\n")
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
                      ranking.by = NULL){


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
               y = i[object$conditions == object$condition.y], alternative = "greater")$statistic
      })


  } else if(object$ranking.by == "wilcox.test"){
    rank_score <-
      apply(object$scmatrix, 1, function(i){
        wilcox.test(x = i[object$conditions == object$condition.x],
                    y = i[object$conditions == object$condition.y], alternative = "greater")$statistic
      })
  } else {
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



#' Perform metabolite set enrichment for bmetenrichr enrichment objects
#'
#' calcEnrichment() ranks metabolites of bmetenrichr object to perform bootstrapping metabolite set enrichment analysis
#'
#' @param object A bmetenrichr object.
#' @param n A integeter describing the number of bootstraps (default = 50).
#'
#' @return An object of class bmetenrich.
#'
#' @examples
#' myTestRun <- calcEnrichment(object = myTestRun, n = 100)
#'
#'
#' @export

#' @export
calcEnrichment <- function (object, ...) {
  UseMethod("calcEnrichment", object)
}

#' @export
calcEnrichment.bmetenrich <- function(object, n = 50){

  options(dplyr.summarise.inform = FALSE)


  if(!all(c(object$condition.x,  object$condition.y) ==   object$rankings$comparison)){
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
  enrichment_analysis <-
      pbapply::pbsapply(seq(n), function(bootstrap_i){
      sapply(names(pathway_list_slim), function(term){

        members_logi <- bootstrapped_sublist[[bootstrap_i]] %in% pathway_list_slim[[term]]
        if(sum(members_logi)==0){
          data.frame(LION_ID = term,
                     bootstrap = bootstrap_i,
                     n = 0,
                     ES = NA,
                     p.value = NA)
        } else if (all(members_logi)){
          data.frame(LION_ID = term,
                     bootstrap = bootstrap_i,
                     n = sum(members_logi),
                     ES = NA,
                     p.value = NA)
        } else {
          ks_results <- ks.test.signed( which(members_logi), which(!members_logi))

          data.frame(LION_ID = term,
                     bootstrap = bootstrap_i,
                     n = sum(members_logi),
                     ES = ks_results$ES,
                     p.value = ks_results$p.value)
        }

      }, simplify = F) %>% bind_rows()
    },simplify = F) %>% bind_rows()


  enrichment_analysis$LION_name <-                    ## match LION name to LION ID
    object$LUT$name[match(enrichment_analysis$LION_ID, object$LUT$ID)]

  object$enrichment_analysis <- list(table = enrichment_analysis,
                                     n = n,
                                     fraction_matched_to_pathway = fraction_matched_to_pathway,
                                     comparison = object$rankings$comparison)

  return(object)

}

#' Plot bootstrap enrichment analysis
#'
#' @param object A bmetenrichr object after enrichment analysis.
#' @param min.annotations An integer describing the minimal number of annotations each term should include
#' @param q.value.cutoff A numeric between 0 and 1. Only terms with q-values lower than this value will be displayed.
#' @param bootstrap.fraction.cutoff A numeric between 0 and 1 (default = 0.5), indicating the minimal fraction that the metabolite set is present in all bootstrap iterations.
#' @param plotIDs A logical indicating whether term IDs should be displayed.
#' @param by.statistic A character indicating how the x-axis will be arranged.
#' Can be either 'ES' (enrichment score) or 'q.value' (default).
#'
#' @return A ggplot2 object.
#' @examples
#'
#' plotEnrichment(myTestRun)
#'
#' @export
plotEnrichment <- function (object, ...) {
  UseMethod("plotEnrichment", object)
}


#' @export
plotEnrichment.bmetenrich <- function(object, min.annotations = 2, q.value.cutoff = 0.1, bootstrap.fraction.cutoff = .5,plotIDs = FALSE, by.statistic = 'q.value'){
  options(dplyr.summarise.inform = FALSE)

  enrichment_analysis <- object$enrichment_analysis$table

  ## remove NAs
  enrichment_analysis <- enrichment_analysis %>% filter(!is.na(p.value), !is.na(ES))


  ## filter terms with lower presence than bootstrap.fraction.cutoff
  enrichment_analysis <-
    enrichment_analysis %>%
    group_by(LION_ID) %>%
    mutate(total_bootstrap = length(bootstrap),
           fraction = total_bootstrap / object$enrichment_analysis$n) %>%
    ungroup() %>%
    filter(fraction > bootstrap.fraction.cutoff)

  enrichment_analysis <-
    enrichment_analysis %>% group_by(LION_ID) %>%
    mutate(n = median(n, na.rm = T))

  if (plotIDs) {
    enrichment_analysis$LION_name <-paste0(enrichment_analysis$LION_name, " (", enrichment_analysis$LION_ID, ")")
  }


  enrichment_analysis <-
    enrichment_analysis %>% group_by(LION_ID) %>%
    mutate(p.value_median = median(p.value, na.rm = T),
           q.value_median = median(p.adjust(p.value, method = "fdr"), na.rm = T))


  switch(by.statistic, 'ES' = {

    enrichment_analysis <-
    enrichment_analysis %>%
      ### here now LION-terms are not filtered by grepl and the names, that's already done by the terms_of_interest step
      filter(n > min.annotations,                                ## only show LION-term with 2 or more molecules, this is still important
             q.value_median < q.value.cutoff,
             LION_ID != "all")                 ## remove LION term 'all'

    if(dim(enrichment_analysis)[1] < 1){
      stop("Not enough enriched terms to visualize")
    }

    enrichment_plot <-
      enrichment_analysis %>%
      {ggplot(data = .,
              aes(x = reorder(LION_name, ES),
                  y = ES,
                  fill = sapply(q.value_median, function(i){  min(10,-log(i, base = 10))})
              ))+

          coord_flip()+
          geom_bar(data = .  %>% group_by(LION_name) %>%
                     summarise(ES = median(ES,na.rm = T),
                               q.value_median = q.value_median[1]),
                   color = NA,                       stat = "identity")+
          geom_jitter(size = .1, width = .1, color = "gray30")+

          scale_fill_gradient2(low = "gray", mid = "gray",high = "red",          ## scale from gray to red, with 10 as max
                               midpoint = -log(0.05, base = 10),limits = c(0,10))+
          geom_hline(yintercept = 0, linetype = 3)+

          labs(x = "", y = "enrichment score", fill = expression(-LOG[10]~italic(q)~value),
               subtitle =  bquote(.(object$enrichment_analysis$comparison[2])~italic(vs.)~.(object$enrichment_analysis$comparison[1]))

          )+
          theme_minimal()+
          theme(plot.title = element_text(face = "bold", hjust = 1), axis.title.x = element_text(face = "bold")
          )
      }
  }, 'q.value' = {

    enrichment_analysis <-
      enrichment_analysis %>% ungroup() %>% group_by(bootstrap) %>%
      mutate(q.value = p.adjust(p = p.value, method = "fdr"))

    enrichment_analysis <-
    enrichment_analysis %>%
      ### here now LION-terms are not filtered by grepl and the names, that's already done by the terms_of_interest step
      filter(n > min.annotations,                                ## only show LION-term with 2 or more molecules, this is still important
             q.value_median < q.value.cutoff,
             LION_ID != "all") %>%                 ## remove LION term 'all'
      group_by(LION_ID) %>%
      mutate(ES = median(ES, na.rm = T),
             up_down = factor(ifelse(sign(ES)>0, "UP","DOWN"), levels = c("UP","DOWN")))

    if(dim(enrichment_analysis)[1] < 1){
      stop("Not enough enriched terms to visualize")
    }

    enrichment_plot <-
      enrichment_analysis  %>%
      {ggplot(data = .,
              aes(x = reorder(LION_name, desc(-log(`q.value`, base = 10))),
                  y = -log(`q.value`, base = 10)
              ))+

          coord_flip()+
           geom_bar(data = . %>% group_by(LION_name) %>%
                     summarise(up_down = up_down[1],
                               q.value = median(-log(`q.value`, base = 10)),
                               q.value_clipped = ifelse(q.value > 10, 10, q.value)),
                   aes(y = q.value, fill = q.value_clipped), stat = "identity"

                   )+
          geom_jitter(size = .1, width = .1, color = "gray30")+
          facet_grid(up_down~.,  space = "free", scales = "free")+

          scale_fill_gradient2(low = "gray", mid = "gray",high = "red",          ## scale from gray to red, with 10 as max
                               midpoint = -log(0.05, base = 10),limits = c(0,10))+
          geom_hline(yintercept = c(0,-log(0.05,base = 10)), linetype = 3)+
          labs(x = "", y = expression(-LOG[10]~italic(q)~value), fill = expression(-LOG[10]~italic(q)~value),
               subtitle =  bquote(.(object$condition.y)~italic(vs.)~.(object$condition.x))
          )+
          theme_minimal()+
          theme(plot.title = element_text(face = "bold", hjust = 1), axis.title.x = element_text(face = "bold")
          )
      }
  },{     ## else
    stop("No valid plot-mode selected. Use  either `ES` or `q.value`")
  })

  return(enrichment_plot)

}

#' Set conditions for enrichment analysis
#'
#' @param object A bmetenrichr object.
#' @param condition.x A optional character describing the reference condition.
#' @param condition.y A optional character describing condition to interest.
#'
#' @return An object of class bmetenrich.
#'
#' @return An object of class bmetenrich.
#' @examples
#'
#' setConditions(myTestRun, condition.x = 'CON', condition.y = "TREATMENT")
#'
#' @export
setConditions <- function (object, ...) {
  UseMethod("setConditions", object)
}

#' @export
setConditions.bmetenrich <- function(object, condition.x = NULL, condition.y = NULL){
  if (is.null(condition.x) & is.null(condition.y)){
    stop("no condition identifiers submitted, with no defaults")
  }

  if (!is.null(condition.x) ){
    if(!any(condition.x == unique(object$conditions))){
      stop("submitted condition not found in dataset")
    } else {
      object$condition.x <- condition.x
    }

  }

  if (!is.null(condition.y) ){
    if(!any(condition.y == unique(object$conditions))){
      stop("submitted condition not found in dataset")
    } else {
      object$condition.y <- condition.y
    }

  }

  return(object)


}

#' Export bootstrap enrichment analysis
#'
#' @param object A bmetenrichr object after enrichment analysis.
#' @param min.annotations An integer describing the minimal number of annotations each term should include (default = 2).
#' @param q.value.cutoff A numeric between 0 and 1. Only terms with q-values lower than this value will be displayed (default = 0.5).
#' @param bootstrap.fraction.cutoff A numeric between 0 and 1 (default = 0.5), indicating the minimal fraction that the metabolite set is present in all bootstrap iterations.
#'
#' @return A data.frame
#' @examples
#'
#' enrichmentTable(myTestRun)
#'
#' @export
enrichmentTable <- function (object, ...) {
  UseMethod("enrichmentTable", object)
}


#' @export
enrichmentTable.bmetenrich <- function(object, min.annotations = 2, q.value.cutoff = 0.5, bootstrap.fraction.cutoff = 0.5){
  options(dplyr.summarise.inform = FALSE)

  enrichment_analysis <- object$enrichment_analysis$table


  ## remove NAs
  enrichment_analysis <- enrichment_analysis %>% filter(!is.na(p.value), !is.na(ES))



  ## filter terms with lower presence than bootstrap.fraction.cutoff
  enrichment_analysis <-
    enrichment_analysis %>%
    group_by(LION_ID) %>%
    mutate(fraction = length(LION_ID) / object$enrichment_analysis$n) %>%
    filter(fraction > bootstrap.fraction.cutoff)


  enrichment_analysis <-
    enrichment_analysis %>% group_by(bootstrap) %>%
    mutate(q.value = p.adjust(p.value, method = "fdr"))  %>%
    group_by(LION_ID, LION_name) %>%
    summarise(n = median(n, na.rm = T),
              ES_median = median(ES, na.rm = T),
              ES_sd = sd(ES, na.rm = T),
              p.value_median = median(p.value, na.rm = T),
              p.value_sd = sd(p.value, na.rm = T),
              q.value_median = median(q.value, na.rm = T),
              q.value_sd = sd(q.value, na.rm = T),
              fraction.bootstrap.presence = median(fraction, na.rm = T)) %>%
    arrange(q.value_median)


  enrichment_analysis <-
    enrichment_analysis %>%
    ### here now LION-terms are not filtered by grepl and the names, that's already done by the terms_of_interest step
    filter(n > min.annotations,                  ## only show LION-term with 2 or more molecules, this is still important
           q.value_median < q.value.cutoff,
           LION_ID != "all"  )   %>%             ## remove LION term 'all')
    ungroup() %>% as.data.frame

  attr(enrichment_analysis, "comparison") <- paste0(object$enrichment_analysis$comparison[2], " vs. ",object$enrichment_analysis$comparison[1])

  return(enrichment_analysis)

}


