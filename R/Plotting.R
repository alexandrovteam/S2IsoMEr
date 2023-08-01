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

  #TODO Double check the p-value adjustment whether to also consider grouping per bootstrap or per term
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
