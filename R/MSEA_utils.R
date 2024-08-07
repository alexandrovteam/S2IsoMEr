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
        #Conditions are switched so that ranking matches the interpretation of logFC and BWS
        calculate_wilcox_statistic(x = i[object$conditions == object$condition.y],
                    y = i[object$conditions == object$condition.x])
        # wilcox.test(x = i[object$conditions == object$condition.y],
        #             y = i[object$conditions == object$condition.x], alternative = alternative)$statistic
      })
  } else if (object$ranking.by == "BWS") {
    rank_score <-
      apply(object$scmatrix, 1, function(i){
        BWStest::bws_test(x = i[object$conditions == object$condition.x],
                          y = i[object$conditions == object$condition.y], alternative = alternative)$statistic
      })
  }
  else if (object$ranking.by == "logFC"){
    cond_x_cells = colnames(object$scmatrix)[which(object$conditions == object$condition.x)]
    cond_y_cells = colnames(object$scmatrix)[which(object$conditions == object$condition.y)]
    # Calculate means for both conditions
    mean_cond_y <- rowMeans(object$scmatrix[, cond_y_cells])
    mean_cond_x <- rowMeans(object$scmatrix[, cond_x_cells])

    # Calculate Log Fold Change (LFC)
    rank_score <- log2(mean_cond_y) - log2(mean_cond_x)
    # rank_score <-
    #
    #   apply(object$scmatrix, 1, function(i){
    #     x = log2(mean(i[object$conditions == object$condition.x], na.rm = T))
    #     y = log2(mean(i[object$conditions == object$condition.y], na.rm = T))
    #
    #     y - x})
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
