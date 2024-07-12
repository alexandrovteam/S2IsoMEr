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
  UseMethod("barplot_bootstrap", object)
}


#' @export
barplot_bootstrap.bmetenrich <- function(object, min.annotations = 2, q.value.cutoff = 0.1,
                                         bootstrap.fraction.cutoff = .5, by.statistic = 'q.value'){
  options(dplyr.summarise.inform = FALSE)

  enrichment_analysis <- object$enrichment_analysis$per_bootstrap_enrich_results

  if (object$gsea.method != "fgsea"){
    colnames(enrichment_analysis)[which(colnames(enrichment_analysis) == "ES")] = "NES"
  }

  ## remove NAs
  enrichment_analysis <- enrichment_analysis %>%
    dplyr::filter(!is.na(p_value), !is.na(NES))


  ## filter terms with lower presence than bootstrap.fraction.cutoff
  if (any(c(object$consider_isomers, object$consider_isobars))){
    if("fraction" %in% colnames(enrichment_analysis)){
      enrichment_analysis = enrichment_analysis %>%
        dplyr::filter(fraction > bootstrap.fraction.cutoff)
    }
    else{
      enrichment_analysis <-
        enrichment_analysis %>%
        dplyr::group_by(Term) %>%
        dplyr::mutate(total_bootstrap = length(bootstrap),
                      fraction = total_bootstrap / object$enrichment_analysis$n) %>%
        dplyr::ungroup() %>%
        dplyr::filter(fraction > bootstrap.fraction.cutoff)
    }
    if ("n" %nin% colnames(enrichment_analysis)){
      enrichment_analysis <-
        enrichment_analysis %>%
        dplyr::group_by(Term) %>%
        dplyr::mutate(n = median(n, na.rm = T))
    }
  }

  enrichment_analysis <-
    enrichment_analysis %>%
    dplyr::group_by(bootstrap) %>%
    dplyr::mutate(q.value = p.adjust(p_value, method = "fdr")) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Term) %>%
    dplyr::mutate(p.value_combined = metap::sumlog(p_value)[["p"]],
                  q.value_combined = metap::sumlog(q.value)[["p"]]) %>%
    dplyr::arrange(q.value_combined) %>%
    dplyr::ungroup()%>%
    dplyr::filter(n > min.annotations,
                  q.value_combined < q.value.cutoff,
                  Term != "all")


  switch(by.statistic, 'ES' = {

    if (dim(enrichment_analysis)[1] < 1) {
      stop("Not enough enriched terms to visualize")
    }

    enrichment_plot <-
      enrichment_analysis %>%
      {
        ggplot(data = .,
               aes(
                 x = reorder(Term, NES),
                 y = NES,
                 fill = sapply(q.value_combined, function(i) {
                   min(10, -log(i, base = 10))
                 })
               )) +

          coord_flip() +
          geom_bar(
            data = .  %>%
              dplyr::group_by(Term) %>%
              dplyr::summarise(NES = median(NES, na.rm = T), q.value_combined = q.value_combined[1]),
            color = NA,
            stat = "identity"
          ) +
          geom_jitter(size = 1,
                      width = .1,
                      color = "gray30") +

          scale_fill_gradient2(
            low = "gray",
            mid = "gray",
            high = "red",
            ## scale from gray to red, with 10 as max
            midpoint = -log(0.05, base = 10),
            limits = c(0, 10)
          ) +
          geom_hline(yintercept = 0, linetype = 3) +

          labs(
            x = "",
            y = "Enrichment Score",
            fill = expression(-LOG[10] ~ italic(q) ~ value),
            subtitle =  bquote(
              .(object$enrichment_analysis$comparison[2]) ~ italic(vs.) ~ .(object$enrichment_analysis$comparison[1])
            )

          ) +
          ggpubr::theme_pubr() +
          theme(
            plot.title = element_text(face = "bold", hjust = 1),
            axis.title.x = element_text(face = "bold")
          )
      }
  },
  'q.value' = {

    enrichment_analysis <-
      enrichment_analysis %>%
      dplyr::group_by(Term) %>%
      mutate(NES = median(NES, na.rm = T),
             up_down = factor(ifelse(sign(NES)>0, "UP","DOWN"),
                              levels = c("UP","DOWN")))

    if(dim(enrichment_analysis)[1] < 1){
      stop("Not enough enriched terms to visualize")
    }

    enrichment_plot <-
      enrichment_analysis  %>%
      ggplot(data = .,aes(x = reorder(Term, desc(-log(`q.value`, base = 10))),
                          y = -log(`q.value`, base = 10))) +
      coord_flip() +
      geom_bar(data = . %>%
                 group_by(Term) %>%
                 summarise(up_down = up_down[1],
                           q.value = median(-log(`q.value`, base = 10)),
                           q.value_clipped = ifelse(q.value > 10, 10, q.value)),
               aes(y = q.value, fill = q.value_clipped), stat = "identity") +
      geom_jitter(size = 1, width = .1, color = "gray30") +
      facet_grid(up_down~.,  space = "free", scales = "free") +
      scale_fill_gradient2(low = "gray", mid = "gray",
                           high = "red",midpoint = -log(0.05, base = 10),
                           limits = c(0,10)) +
      geom_hline(yintercept = c(0,-log(0.05,base = 10)), linetype = 3) +
      labs(x = "", y = expression(-LOG[10]~italic(q)~value),
           fill = expression(-LOG[10]~italic(q)~value),
           subtitle =  bquote(.(object$condition.y)~italic(vs.)~.(object$condition.x))) +
      ggpubr::theme_pubr(legend = "right") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title.x = element_text(face = "bold"))
  },
  {     ## else
    stop("No valid plot-mode selected. Use  either `ES` or `q.value`")
  })

  return(enrichment_plot)

}

barplot_ORA_simple = function(ORA_simple_res, q_val_cutoff = 0.05){
  plot_data = ORA_simple_res
  plot_data$Significance = ifelse(plot_data$score < -log10(q_val_cutoff), "No", "Yes")

  plot_data %>%
    ggplot(aes(x = reorder(Term, score, max), y = score, fill = Significance)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = as.character(TP),hjust = 1), position = position_dodge(width = 0.5),
              color = "black", size = 6) +
    coord_flip() +
    ggpubr::theme_pubr() +
    ylab("-Log10 (q-value)") +
    xlab("") +
    geom_hline(yintercept = -log10(q_val_cutoff), linetype = "dashed", lwd = 1) +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14))
}

barplot_ORA_boot = function(ORA_boot_res){
  boot_summary = ORA_boot_res[["unfiltered_enrich_res"]] %>%
    dplyr::filter(Term %fin% ORA_boot_res[["clean_enrich_res"]]$Term) %>%
    dplyr::group_by(Term) %>%
    dplyr::summarise(ES_max = max(OR),
                     ES_min = min(OR))


  data_viz = ORA_boot_res[["clean_enrich_res"]] %>%
    dplyr::left_join(boot_summary) %>%
    dplyr::mutate(n = ceiling(n),
                  Term_label = paste0(Term, " (", n, ")")) %>%
    dplyr::select(Term_label, ES_median, ES_max, ES_min,
                  p.value_combined, q.value_combined)

  min_qval = min(data_viz$q.value_combine)

  plot <- ggplot(data_viz, aes(x = reorder(Term_label, ES_median,),
                               y = ES_median, fill = q.value_combined)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
    scale_y_continuous(n.breaks = 10) +
    geom_errorbar(aes_string(ymin = "ES_min", ymax = "ES_max"),
                  width = 0.2, position = position_dodge(0.7)) +
    scale_fill_gradient(low = "#B3CDE3", high = "#FBB4AE", name = "q.value_combined") +
    ggpubr::theme_pubr() +
    coord_flip() +
    xlab("") +
    ylab("Fold Enrichment Score")
    # theme(plot.title = element_text(hjust = 0.5))

    return(plot)
}

dotplot_enrich = function(plot_data, alpha_cutoff = 0.05){
  cols = c("blue", "red")

  colnames(plot_data)[which(colnames(plot_data) == "TP")] = "size"
  colnames(plot_data)[which(colnames(plot_data) == "n")] = "size"
  colnames(plot_data)[which(colnames(plot_data) == "ES_median")] = "Enrichment_Score"
  colnames(plot_data)[which(colnames(plot_data) == "score")] = "Enrichment_Score"
  colnames(plot_data)[which(colnames(plot_data) == "p_value")] = "sig_score"
  colnames(plot_data)[which(colnames(plot_data) == "p.value_combined")] = "sig_score"

  plot_data[["sig_score"]] = -log10(plot_data[["sig_score"]])

  p = ggplot(data = plot_data, aes(x = Enrichment_Score,
                               y = reorder(Term, Enrichment_Score))) +
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


ridge_bootstraps = function(enrich_res, terms_of_interest){

  plot_data = enrich_res[enrich_res$Term %in% terms_of_interest,]
  colnames(plot_data)[which(colnames(plot_data) == "n")] = "TP"

  plot_data$Term = paste0(plot_data$Term, " (", plot_data$fraction * 100,
                          ")")

  ggplot2::ggplot(plot_data, aes(x = TP, y = Term,
                        fill = factor(stat(quantile)))) +
    ggridges::stat_density_ridges(geom = "density_ridges_gradient",calc_ecdf = T,
                                  quantiles = 4, quantile_lines = T) +
    ggplot2::scale_fill_viridis_d(name = "Quartiles") +
    ggplot2::scale_x_continuous(n.breaks = 15, limits = c(0,max(plot_data$TP) +2)) +
    ylab("") +
    xlab("Intersection size") +
    ggpubr::theme_pubr()
}

#' Plot intensity distribution of metabolite across conditions
#'
#' @param object A bmetenrichr object after enrichment analysis.
#' @param metabolite A character indicating metabolite name
#'
#' @return A ggplot2 object.
#' @examples
#'
#' compare_metabo_distr(obj, metabolite="C6H8N2O3-H")
#'
#' @export
compare_metabo_distr = function(obj, metabolite){
  data = obj$scmatrix[metabolite,] %>%
    as.data.frame()
  colnames(data)[1] = "intens"
  data$condition = obj$conditions
  # data_viz = data %>%
  #   tidyr::gather(key = "metabo", value = "intens",-condition)

  desc_stats <- data_viz %>%
    dplyr::group_by(condition) %>%
    dplyr::summarise(
      mean = mean(intens),
      median = median(intens)
    )

  plot <- ggplot(data_viz, aes(x = intens, fill = condition)) +
    geom_density(alpha = 0.5) +
    geom_density(alpha = 0.4) +
    geom_vline(data = desc_stats, aes(xintercept = mean, color = condition, linetype = "Mean"), size = 1) +
    geom_vline(data = desc_stats, aes(xintercept = median, color = condition, linetype = "Median"), size = 1) +
    scale_linetype_manual(name = "Statistic", values = c("Mean" = "dashed", "Median" = "solid")) +
    labs(title = metabolite,
         x = "Intensity(Log10)",
         y = "Density",
         fill = "Condition",
         color = "Condition",
         linetype = "Statistic") +
    ggpubr::theme_pubr() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(n.breaks = 10)

  return(plot)

}

