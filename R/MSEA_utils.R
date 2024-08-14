#' Rank metabolites for S2IsoMEr enrichment object
#'
#' rankScore() ranks metabolites of S2IsoMEr object to perform bootstrapping metabolite set enrichment analysis
#'
#' @param object A S2IsoMEr object.
#' @param ranking.by A character string specifying the ranking method. Options include:
#'   \itemize{
#'     \item \code{"wilcox.test"} (default) - Ranks metabolites based on the results of a Wilcoxon test.
#'     \item \code{"t.test"} - Ranks metabolites based on the results of a t-test.
#'     \item \code{"logFC"} - Ranks metabolites based on log fold change.
#'     \item \code{"BWS"} - Ranks metabolites based on BWS (Baumgartner-Weiss-Schindler ).
#'   }
#' @param alternative A character string specifying the alternative hypothesis for statistical tests. Options are:
#'   \itemize{
#'     \item \code{"two.sided"} - Test for differences in both directions.
#'     \item \code{"less"} - Test if the second condition is less than the first
#'     \item \code{"greater"} - Test if the second condition is greater than the first
#'   }
#'   Default is \code{"greater"}.
#'
#' @return An object of class S2IsoMEr.
#'
#' @examples
#' myTestRun <-
#' rankScore(object = myTestRun, ranking.by = 't.test')
#'
#'
#' @export
rankScore <- function (object,
                       ranking.by = "wilcox.test",
                       alternative = "greater") {
  UseMethod("rankScore")
}

#' @export
rankScore.S2IsoMEr <- function(object,
                                 ranking.by = NULL,
                                 alternative = "greater"){


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
    #Conditions are switched so that ranking matches the interpretation of logFC
    rank_score <-
      apply(object$scmatrix, 1, function(i){
        t.test(x = i[object$conditions == object$condition.y],
               y = i[object$conditions == object$condition.x], alternative = alternative)$statistic
      })


  } else if(object$ranking.by == "wilcox.test"){

    LFC = calc_LFC_scmat(object)

    rank_score <-
      apply(object$scmatrix, 1, function(i){
        #Conditions are switched so that ranking matches the interpretation of logFC
        calculate_wilcox_statistic(x = i[object$conditions == object$condition.y],
                    y = i[object$conditions == object$condition.x])
        # wilcox.test(x = i[object$conditions == object$condition.y],
        #             y = i[object$conditions == object$condition.x], alternative = alternative)$statistic
      })

    rank_score = rank_score * sign(LFC)

  } else if (object$ranking.by == "BWS") {
    LFC = calc_LFC_scmat(object)
    #Conditions are switched so that ranking matches the interpretation of logFC and BWS
    rank_score <-
      apply(object$scmatrix, 1, function(i){
        BWStest::bws_stat(x = i[object$conditions == object$condition.y],
                          y = i[object$conditions == object$condition.x])
      })
    rank_score = rank_score * sign(LFC)
  }
  else if (object$ranking.by == "logFC"){
    rank_score <- calc_LFC_scmat(object)

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


calculate_wilcox_statistic <- function(x, y = NULL, paired = FALSE) {
  if (!is.null(y)) {
    if (paired) {
      if (length(x) != length(y))
        stop("'x' and 'y' must have the same length")
      x <- x - y
      y <- NULL
    } else {
      y <- y[!is.na(y)]
    }
  }

  x <- x[!is.na(x)]

  if (is.null(y)) {
    # Wilcoxon signed rank test
    r <- rank(abs(x))
    STATISTIC <- sum(r[x > 0])
  }
  else {
    # Wilcoxon rank sum test
    r <- rank(c(x, y))
    STATISTIC <- sum(r[seq_along(x)]) - length(x) * (length(x) + 1) / 2
  }

  return(STATISTIC)
}

calc_LFC_scmat <- function (object) {
  UseMethod("calc_LFC_scmat")
}
calc_LFC_scmat.S2IsoMEr = function(object){
  cond_x_cells = colnames(object$scmatrix)[which(object$conditions == object$condition.x)]
  cond_y_cells = colnames(object$scmatrix)[which(object$conditions == object$condition.y)]
  # Calculate means for both conditions
  mean_cond_y <- rowMeans(object$scmatrix[, cond_y_cells])
  mean_cond_x <- rowMeans(object$scmatrix[, cond_x_cells])

  # Calculate Log Fold Change (LFC)
  LFC <- log2(mean_cond_y) - log2(mean_cond_x)
  return(LFC)
}

