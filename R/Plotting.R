
# MSEA -----------------------------------------------------------
#' Plot bootstrap enrichment analysis
#'
#' @param object A S2IsoMEr object after enrichment analysis.
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
barplot_MSEA_boot <- function (object, ...) {
  UseMethod("barplot_MSEA_boot", object)
}


#' @export
barplot_MSEA_boot.S2IsoMEr <- function(object, min.annotations = 2, q.value.cutoff = 0.1,
                                         bootstrap.fraction.cutoff = .5, by.statistic = 'ES'){
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
               aes(x = reorder(Term, NES), y = NES)) +
          coord_flip() +
          geom_bar(data = . %>% group_by(Term) %>%
                     summarise(NES = median(NES, na.rm = TRUE),
                               q.value_combined = first(q.value_combined)),
                   aes(fill = -log10(q.value_combined)),
                   color = NA, stat = "identity", alpha = 0.5) +
          geom_jitter(size = 1, width = 0.1, color = "gray30") +
          scale_fill_viridis_c(option = "viridis",
                               limits = c(0, max(-log10(.$q.value_combined), na.rm = TRUE)),
                               direction = 1) +
          geom_hline(yintercept = 0, linetype = "dotted") +
          labs(x = "", y = "Enrichment Score",
               fill = expression(-log[10]~italic(q)~value),
               subtitle = bquote(.(object$enrichment_analysis$comparison[2]) ~
                                   italic(vs.) ~ .(object$enrichment_analysis$comparison[1]))) +
          ggpubr::theme_pubr() +
          theme(plot.title = element_text(face = "bold", hjust = 0.5),
                axis.title.x = element_text(face = "bold"))
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

#' Plot Heatmap of MSEA Results Across Multiple Conditions
#'
#' @description This function generates a heatmap of combined MSEA results. The heatmap displays normalized enrichment scores (NES) across multiple conditions and highlights significant results based on a specified alpha cutoff.
#'
#' @param combined_MSEA_res A data frame containing combined MSEA results for each pairwise comparison from either enrichment results in \code{\link{Run_bootstrap_MSEA}} or \code{\link{Run_simple_MSEA}}
#' @param alpha_cutoff A numeric value specifying the cutoff for significance. Terms with p-values below this threshold will be marked as significant. Default is 0.05.
#'
#' @return A heatmap plot showing the NES values for each term across different condition pairs. Significant results are highlighted with an asterisk.
#'
#' @details The function checks for the presence of required columns and ensures there are enough conditions to generate a meaningful plot. The heatmap includes color coding for NES values, with significant results indicated by asterisks.
#'
#' @examples
#' \dontrun{
#'   plot_MSEA_Multi_cond(combined_MSEA_res = msea_results, alpha_cutoff = 0.05)
#' }
#'
#' @import pheatmap
#' @import RColorBrewer
#' @import grDevices
#' @export
plot_MSEA_Multi_cond = function(combined_MSEA_res,
                                alpha_cutoff = 0.05){

  missing_conds <- setdiff(c("condition.x", "condition.y"),
                           colnames(combined_MSEA_res))

  if (length(missing_conds) > 0) {
    stop(paste("The following columns are missing:", paste(missing_conds,
                                                           collapse = ", ")))
  }

  if (nrow(unique(combined_MSEA_res[c("condition.x", "condition.y")])) < 2) {
    stop("Not enough conditions to plot heatmap")
  }

  combined_MSEA_res <- combined_MSEA_res %>%
    rename_if_exists("n", "size") %>%
    rename_if_exists("pathway", "Term") %>%
    rename_if_exists("ES_median", "NES") %>%
    rename_if_exists("pval", "p_value") %>%
    rename_if_exists("p.value_combined", "p_value")

  combined_MSEA_res$signif = ifelse(combined_MSEA_res$p_value < alpha_cutoff,
                                    "*", "")

  heat_mat = combined_MSEA_res %>%
    dplyr::mutate(cond_label = paste0(condition.y, "_",condition.x)) %>%
    dplyr::select(cond_label, Term, NES) %>%
    tidyr::pivot_wider(names_from = cond_label, values_from = NES) %>%
    column_to_rownames("Term") %>%
    as.matrix()

  label_mat = combined_MSEA_res %>%
    dplyr::mutate(cond_label = paste0(condition.y, "_",condition.x)) %>%
    dplyr::select(cond_label, Term, signif) %>%
    tidyr::pivot_wider(names_from = cond_label, values_from = signif) %>%
    column_to_rownames("Term") %>%
    as.matrix()

  pheatmap::pheatmap(heat_mat,
                     color = grDevices::colorRampPalette(rev(
                       RColorBrewer::brewer.pal(n = 10, name ="RdBu")))(100),
                     cluster_rows = T, cluster_cols = F, na_col = "white",
                     treeheight_row = 5,
                     show_colnames = T,
                     display_numbers = label_mat,
                     fontsize = 14, angle_col = 45,
                     fontsize_number = 14, number_color = "lightblue",
                     legend = T,legend_breaks = c(-2,-1,0, 1, 2),
                     legend_labels = c("-2", "-1", "0", "1", "NES"))
}

# ORA ------------------------------------------------------------
#' Bar Plot for Simple ORA Results
#'
#' @description This function creates a bar plot for simple (no bootstraps) Over-Representation Analysis (ORA) results, indicating the significance of terms based on a q-value cutoff.
#'
#' @param ORA_simple_res A data frame containing the ORA results.
#' @param q_val_cutoff A numeric value specifying the q-value cutoff for significance. Default is 0.05.
#'
#' @return A ggplot object displaying the bar plot of ORA results.
#'
#' @details The function creates a bar plot where the terms are reordered by their scores in descending order. Bars are colored based on whether the term's score is above or below the negative log10 of the q-value cutoff. A dashed line indicates the cutoff for significance.
#'
#' @examples
#' \dontrun{
#'   plot <- barplot_ORA_simple(ORA_simple_res = ora_results, q_val_cutoff = 0.05)
#'   print(plot)
#' }
#'
#' @import ggplot2
#' @import ggpubr
#' @export
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

#' Bar Plot for ORA Bootstrapped Results
#'
#' @description This function creates a bar plot for Over-Representation Analysis (ORA) bootstrapped results, including error bars to show the range of enrichment scores (OR) and indicating the significance based on q-values.
#'
#' @param ORA_boot_res A list containing the ORA bootstrapped results. The list should include `unfiltered_enrich_res` and `clean_enrich_res` data frames.
#'
#' @return A ggplot object displaying the bar plot of ORA bootstrapped results.
#'
#' @details The function creates a bar plot where the terms are reordered by their median enrichment scores. Bars are filled with colors representing the combined q-values, and error bars show the minimum and maximum enrichment scores (ES_min and ES_max) across bootstraps. The number of true positives (n) is included in the term labels.
#'
#' @examples
#' \dontrun{
#'   plot <- barplot_ORA_boot(ORA_boot_res = ora_boot_results)
#'   print(plot)
#' }
#'
#' @import ggplot2
#' @import ggpubr
#' @import dplyr
#' @export
barplot_ORA_boot = function(ORA_boot_res, collapse_multi_cond = F){
  if (collapse_multi_cond){
    ORA_boot_res = collapse_ORA_boot_multi_cond(ORA_boot_res_list = ORA_boot_res)
  }

  ORA_boot_res = ORA_boot_res[c("unfiltered_enrich_res",
                                "clean_enrich_res")]

  ORA_boot_res = lapply(ORA_boot_res, function(x){
    if ("Condition" %in% colnames(x)){
      x
    }
    else{
      x = x %>% dplyr::mutate(Condition = "Query")
      x
    }
  })

  boot_summary = ORA_boot_res[["unfiltered_enrich_res"]] %>%
    dplyr::filter(Term %fin% ORA_boot_res[["clean_enrich_res"]]$Term) %>%
    dplyr::group_by(Term,Condition) %>%
    dplyr::summarise(ES_max = max(OR),
                     ES_min = min(OR))


  data_viz = ORA_boot_res[["clean_enrich_res"]] %>%
    dplyr::left_join(boot_summary) %>%
    dplyr::mutate(n = ceiling(n),
                  Term_label = paste0(Term, " (", n, ")")) %>%
    dplyr::select(Term_label, ES_median, ES_max, ES_min,
                  p.value_combined, q.value_combined,Condition)

  min_qval = min(data_viz$q.value_combine)

  plot <- ggplot(data_viz, aes(x = reorder(Term_label, ES_median,),
                               y = ES_median, fill = -log10(q.value_combined))) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
    scale_y_continuous(n.breaks = 10) +
    geom_errorbar(aes_string(ymin = "ES_min", ymax = "ES_max"),
                  width = 0.2, position = position_dodge(0.7)) +
    scale_fill_gradient(low = "#B3CDE3", high = "#FBB4AE", name = "-Log10(q.value)") +
    ggpubr::theme_pubr() +
    coord_flip() +
    xlab("") +
    ylab("Fold Enrichment Score") +
    facet_wrap(~Condition)
    # theme(plot.title = element_text(hjust = 0.5))

    return(plot)
}

#' Dot Plot for Over-Representation Analysis (ORA) Results
#'
#' @description This function creates a dot plot for ORA results, with options to color the dots by enrichment score or significance and to facet by a specified variable.
#'
#' @param ORA_res A data frame containing the ORA results from either filtered results in output of \code{\link{Run_bootstrap_ORA}} or direct output of \code{\link{Run_simple_ORA}}.
#' @param alpha_cutoff A numeric value specifying the alpha cutoff for significance. Default is 0.05.
#' @param color_by A character string indicating whether to color the dots by "Enrichment_score" or "Significance"
#' @param facet_by An optional character string specifying a column name by which to facet the plot.
#'
#' @return A ggplot object displaying the dot plot of ORA results.
#'
#' @details The function creates a dot plot where the size of the dots represents the size of the term (e.g., number of true positives), and the color represents either the enrichment score or the significance (-log10 of the p-value). If multiple conditions are present, the plot is adjusted accordingly. The plot can also be faceted by a specified variable.
#'
#' @examples
#' \dontrun{
#'   plot <- dotplot_ORA(plot_data = ora_results, alpha_cutoff = 0.05, color_by = "Significance", facet_by = "condition")
#'   print(plot)
#' }
#'
#' @import ggplot2
#' @import ggpubr
#' @import dplyr
#' @import DOSE
#' @export
dotplot_ORA = function(ORA_res, alpha_cutoff = 0.05,
                       color_by = "Significance",
                       facet_by = NULL){
  cols = c("#3B6FB6", "#D41645")

  plot_data <- ORA_res %>%
    rename_if_exists("TP", "size") %>%
    rename_if_exists("n", "size") %>%
    rename_if_exists("ES_median", "Enrichment_Score") %>%
    rename_if_exists("score", "Enrichment_Score") %>%
    rename_if_exists("p_value", "sig_score") %>%
    rename_if_exists("p.value_combined", "sig_score") %>%
    rename_if_exists("condition", "Condition")

  if (color_by == "Enrichment_score"){
    plot_data[["Score"]] = plot_data[["Enrichment_score"]]
  }
  else if (color_by == "Significance"){
    plot_data[["Score"]] = -log10(plot_data[["sig_score"]])
  }
  else{
    stop("Invalid argument for color_by. Select either Enrichment_score or Significance")
  }

  label_func = dotplot_label_func(n = 30)

  color_legend_title = switch(color_by, "Significance" = '-Log10(pvalue)',
                              "Enrichment_score" = 'Enrichment Score')

  if(length(unique(plot_data$Condition)) > 1){
    p = ggplot(plot_data, aes_string(x = "Condition",
                                     y = "Term", size = "size")) +
      geom_point(aes_string(color = "Score")) +
      scale_color_continuous(low = "#3B6FB6",
                             high = "#D41645", guide = guide_colorbar(reverse = TRUE)) +
      ylab(NULL) +
      # ggtitle("Trial") +
      DOSE::theme_dose(12) +
      scale_size_continuous(range = c(4, 9)) +
      scale_y_discrete(labels = label_func) +
      ggpubr::rotate_x_text(angle = 45) +
      guides(size = guide_legend(order = 1, title = 'Term/query overlap'),
             color = guide_colorbar(order = 2, title = color_legend_title)) +
      theme(axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 16),
            strip.text = element_text(size = 14))
  }
  else{
    p = ggplot(data = plot_data, aes(x = Enrichment_Score,
                                     y = reorder(Term, Enrichment_Score))) +
      geom_point(aes(size = size,color = Score)) +
      scale_size(range = c(3, 6)) +
      scale_color_gradient(low = cols[1], high = cols[2]) +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
      guides(size = guide_legend(title = 'Term/query overlap'),
             color = guide_colorbar(title = color_legend_title)) +
      labs(
        x = 'Enrichment Score',
        y = ''
      ) +
      ggpubr::theme_pubr() +
      theme(axis.text.x = element_text(angle = 0,hjust = 1))
  }

  if (!is.null(facet_by)){
    facet_col_sym = sym(facet_by)
    p = p +
      facet_wrap(vars(!!enquo(facet_col_sym))) +
      theme(
        panel.spacing = unit(x = 1, units = "lines"),
        strip.background = element_blank())
  }
  return(p)
}


# Generic/other --------------------------------------------------------

#' Ridge Plot of Bootstrap Enrichment Results for Terms of Interest
#'
#' @description This function creates a ridge plot visualizing the distribution of intersection sizes (TP) for specified terms of interest from bootstrap enrichment results.
#'
#' @param enrich_res A data frame containing the enrichment results per bootstrap, including columns for Term, n (TP), and fraction.
#' @param terms_of_interest A character vector specifying the terms of interest to be plotted.
#'
#' @return A ggplot object displaying the ridge plot of intersection sizes for the specified terms of interest.
#'
#' @details The function filters the enrichment results to include only the specified terms of interest and creates a ridge plot showing the distribution of intersection sizes (TP). The terms are labeled with their names and the fraction of bootstraps in which they appear, and the plot includes quantile lines to indicate distribution quartiles.
#'
#' @examples
#' \dontrun{
#'   plot <- ridge_bootstraps(enrich_res = enrichment_results,
#'                            terms_of_interest = c("Term1", "Term2", "Term3"))
#'   print(plot)
#' }
#'
#' @import ggplot2
#' @import ggridges
#' @import ggpubr
#' @export
ridge_bootstraps = function(enrich_res, terms_of_interest, condition = NULL){


  plot_data = enrich_res[enrich_res$Term %in% terms_of_interest,]
  colnames(plot_data)[which(colnames(plot_data) == "n")] = "TP"

  plot_data$Term = paste0(plot_data$Term, " (", plot_data$fraction * 100,
                          ")")

  if (!is.null(condition)){
    plot_data = plot_data[plot_data$Condition == condition,]
  }

  ggplot2::ggplot(plot_data, aes(x = TP, y = Term,
                                 fill = stat(density))) +
    ggridges::stat_density_ridges(geom = "density_ridges_gradient",calc_ecdf = T,
                                  scale = 0.95) +
    ggplot2::scale_fill_viridis_c(name = "Density", direction = 1) +
    ggplot2::scale_x_continuous(n.breaks = 15, limits = c(0,max(plot_data$TP) +2)) +
    ylab("") +
    xlab("Intersection size") +
    ggpubr::theme_pubr()
}

#' Plot intensity distribution of metabolite across conditions
#'
#' @param object A S2IsoMEr object after enrichment analysis.
#' @param metabolite A character indicating metabolite name
#'
#' @return A ggplot2 object.
#' @examples
#'
#' compare_metabo_distr(obj, metabolite="C6H8N2O3-H")
#'
#' @export
compare_metabo_distr = function(obj, metabolite, conds_of_interest = NULL){
  data = obj$scmatrix[metabolite,] %>%
    as.data.frame()
  colnames(data)[1] = "intens"
  data$condition = obj$conditions

  data$intens = convert_to_log10(data$intens)

  data_viz = data
  # data_viz = data %>%
  #   tidyr::gather(key = "metabo", value = "intens",-condition)

  if(!is.null(conds_of_interest)){
    data_viz = data_viz %>%
      dplyr::filter(condition %in% conds_of_interest)
  }

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
    scale_x_continuous(n.breaks = 10) +
    facet_wrap(~condition, scales = "free_y")

  return(plot)

}

