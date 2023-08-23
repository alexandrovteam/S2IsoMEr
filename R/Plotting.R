#' Plot bootstrap enrichment analysis
#'
#' @param object A bmetenrichr object after enrichment analysis.
#' @param min.annotations An integer describing the minimal number of annotations each term should include
#' @param q.value.cutoff A numeric between 0 and 1. Only terms with q-values lower than this value will be displayed.
#' @param bootstrap.fraction.cutoff A numeric between 0 and 1 (default = 0.5), indicating the minimal fraction that the metabolite set is present in all bootstrap iterations.
#' @param by.statistic A character indicating how the x-axis will be arranged.
#' Can be either 'ES' (enrichment score) or 'q.value' (default).
#'
#' @return A ggplot2 object.
#' @examples
#'
#' barplot_bootstrap(myTestRun)
#'
#' @export
barplot_bootstrap <- function (object, ...) {
  UseMethod("plotEnrichment", object)
}


#' @export
barplot_bootstrap.bmetenrich <- function(object, min.annotations = 2, q.value.cutoff = 0.1,
                                      bootstrap.fraction.cutoff = .5, by.statistic = 'q.value'){
  options(dplyr.summarise.inform = FALSE)

  enrichment_analysis <- object$enrichment_analysis$table

  if (object$gsea.method != "fgsea"){
    colnames(enrichment_analysis)[which(colnames(enrichment_analysis) == "ES")] = "NES"
  }

  ## remove NAs
  enrichment_analysis <- enrichment_analysis %>% dplyr::filter(!is.na(p.value), !is.na(NES))


  ## filter terms with lower presence than bootstrap.fraction.cutoff
  if (any(c(object$consider_isomers, object$consider_isobars))){
    if("fraction" %in% colnames(enrichment_analysis)){
      enrichment_analysis = enrichment_analysis %>% dplyr::filter(fraction > bootstrap.fraction.cutoff)
    }
    else{
      enrichment_analysis <-
        enrichment_analysis %>%
        group_by(Term) %>%
        mutate(total_bootstrap = length(bootstrap),
               fraction = total_bootstrap / object$enrichment_analysis$n) %>%
        ungroup() %>%
        filter(fraction > bootstrap.fraction.cutoff)
    }
    if ("n" %nin% colnames(enrichment_analysis)){
      enrichment_analysis <-
        enrichment_analysis %>% group_by(Term) %>%
        mutate(n = median(n, na.rm = T))
    }
  }

  #TODO Double check the p-value adjustment whether to also consider grouping per bootstrap or per term
  enrichment_analysis <-
    enrichment_analysis %>% group_by(Term) %>%
    mutate(p.value_median = median(p.value, na.rm = T),
           q.value_median = median(p.adjust(p.value, method = "fdr"), na.rm = T))


  switch(by.statistic, 'ES' = {

    enrichment_analysis <-
      enrichment_analysis %>%
      dplyr::filter(n > min.annotations,
             q.value_median < q.value.cutoff,
             Term != "all")

    if(dim(enrichment_analysis)[1] < 1){
      stop("Not enough enriched terms to visualize")
    }

    enrichment_plot <-
      enrichment_analysis %>%
      {ggplot(data = .,
              aes(x = reorder(Term, NES),
                  y = NES,
                  fill = sapply(q.value_median, function(i){  min(10,-log(i, base = 10))})
              ))+

          coord_flip()+
          geom_bar(data = .  %>% dplyr::group_by(Term) %>%
                     dplyr::summarise(NES = median(NES,na.rm = T),q.value_median = q.value_median[1]),
                   color = NA,stat = "identity" )+
          geom_jitter(size = .1, width = .1, color = "gray30")+

          scale_fill_gradient2(low = "gray", mid = "gray",high = "red",          ## scale from gray to red, with 10 as max
                               midpoint = -log(0.05, base = 10),limits = c(0,10))+
          geom_hline(yintercept = 0, linetype = 3)+

          labs(x = "", y = "Enrichment Score", fill = expression(-LOG[10]~italic(q)~value),
               subtitle =  bquote(.(object$enrichment_analysis$comparison[2])~italic(vs.)~.(object$enrichment_analysis$comparison[1]))

          )+
          ggpubr::theme_pubr()+
          theme(plot.title = element_text(face = "bold", hjust = 1), axis.title.x = element_text(face = "bold")
          )
      }
  }, 'q.value' = {

    enrichment_analysis <-
      enrichment_analysis %>% ungroup() %>% group_by(bootstrap) %>%
      mutate(q.value = p.adjust(p = p.value, method = "fdr"))

    enrichment_analysis <-
      enrichment_analysis %>%
      dplyr::filter(n > min.annotations,
             q.value_median < q.value.cutoff,
             Term != "all") %>%
      dplyr::group_by(Term) %>%
      mutate(NES = median(NES, na.rm = T),
             up_down = factor(ifelse(sign(NES)>0, "UP","DOWN"), levels = c("UP","DOWN")))

    if(dim(enrichment_analysis)[1] < 1){
      stop("Not enough enriched terms to visualize")
    }

    enrichment_plot <-
      enrichment_analysis  %>%
      {ggplot(data = .,
              aes(x = reorder(Term, desc(-log(`q.value`, base = 10))),
                  y = -log(`q.value`, base = 10)
              ))+

          coord_flip()+
          geom_bar(data = . %>% group_by(Term) %>%
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
          ggpubr::theme_pubr()()+
          theme(plot.title = element_text(face = "bold", hjust = 1), axis.title.x = element_text(face = "bold")
          )
      }
  },{     ## else
    stop("No valid plot-mode selected. Use  either `ES` or `q.value`")
  })

  return(enrichment_plot)

}


#TODO Write a function for barplot simple

dotplot_enrich = function(plot_data, alpha_cutoff = 0.05){
  cols = c("blue", "red")

  colnames(plot_data)[which(colnames(plot_data) == "NES")] = "Enrichment_Score"
  colnames(plot_data)[which(colnames(plot_data) == "pval")] = "sig_score"

  plot_data[["sig_score"]] = -log10(plot_data[["sig_score"]])

  p = ggplot(data = plot_data, aes(x = Enrichment_Score,
                               y = fct_reorder(term, Enrichment_Score))) +
    geom_point(aes(size = size,color = sig_score)) +
    scale_size(range = c(3, 6)) +
    scale_color_gradient(low = cols[1], high = cols[2]) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(size = guide_legend(title = 'Term/query overlap'),
           color = guide_colorbar(title = '-Log10(pvalue)')) +
    labs(
      x = 'Enrichment Score',
      y = ''
    ) +
    ggpubr::theme_pubr() +
    theme(axis.text.x = element_text(angle = 0,hjust = 1))
    # facet_grid(~section, scales = "free_x", space = "free_x", switch = "y") +
    # theme(
    #   panel.spacing = unit(x = 1, units = "lines"),
    #   strip.background = element_blank()) +
  return(p)
}

Multi_cond_heatmap = function(plot_data){
  if ("condition" %in% colnames(plot_data)){
    if (length(unique(plot_data$condition)) < 2){
      stop("Not enough conditions to plot heatmap")
    }
  }
  else{
    stop("Data is missing condition info")
  }

  if ("TP" %in% colnames(plot_data)){
    colnames(plot_data)[which(colnames(plot_data) == "TP")] = "size"
  }

  colnames(plot_data)[which(colnames(plot_data) == "p_value")] = "pval"

  plot_data$sig = ifelse(plot_data$pval < alpha_cutoff, "*", "")

  plot_data$hm_label = paste0(plot_data$size, "\n", plot_data$sig)

  long_df = plot_data %>% dplyr::select(condition, term, score)
  long_df = long_df[!duplicated(long_df),]

  long_df_label = plot_data %>% dplyr::select(condition, term, hm_label)
  long_df_label = long_df_label[!duplicated(long_df_label),]

  wide_df = long_df %>% tidyr::pivot_wider(names_from = condition,
                                    values_from = score)
  rNames = wide_df$term
  wide_df = wide_df[,-1] %>% as.matrix()
  rownames(wide_df) = rNames
  wide_df[is.na(wide_df)] = 0

  wide_df_label = long_df_label %>% tidyr::pivot_wider(names_from = condition,
                                                values_from = hm_label)
  rNames = wide_df_label$term
  wide_df_label = wide_df_label[,-1] %>% as.matrix()
  rownames(wide_df_label) = rNames
  wide_df_label[is.na(wide_df_label)] = ""


  pheatmap::pheatmap(mat = wide_df, cluster_rows = T, cluster_cols = F,
           color = grDevices::colorRampPalette(rev(
             RColorBrewer::brewer.pal(n = 10, name ="RdBu"))
             )(100),
           treeheight_row = 5,show_colnames = T,
           display_numbers = wide_df_label, fontsize = 14, angle_col = 45,
           fontsize_number = 14, number_color = "lightblue",
           legend = T)
}

#' @importFrom EnhancedVolcano EnhancedVolcano
DE_volcano = function(plot_data){
  #TODO prepare data for enhancedvolcano plot
}


ridge_bootstraps = function(plot_data){

  colnames(plot_data)[which(colnames(plot_data) == "n")] = "TP"
  colnames(plot_data)[which(colnames(plot_data) == "term")] = "Term"

  plot_data$Term = paste0(plot_data$Term, " (", plot_data$fraction * 100,
                          ")")

  ggplot2::ggplot(plot_data, aes(x = TP, y = Term,
                        fill = -log10(p_value))) +
    ggridges::stat_density_ridges(geom = "density_ridges_gradient") +
    ggplot2::scale_fill_viridis_c(name = "-Log10(pval)", direction = 1)
}

treeplot_enrich = function(plot_data){
  #TODO Write function and move to utils for pair_term_semsim for tree plot
  #TODO Write script to prepare data for tree plot
  #TODO Rewrite the tree plot or import it from enrichplot R package
  #https://rdrr.io/github/GuangchuangYu/enrichplot/src/R/treeplot.R
}


