
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
#' \dontrun{
#' data("example_MSEA_obj")
#' p = barplot_MSEA_boot(object = example_MSEA_obj,
#'       q.value.cutoff = 0.2, by.statistic = "ES")
#' }
#'
#' @export
barplot_MSEA_boot <- function (object, min.annotations = 2, q.value.cutoff = 0.1,
                               bootstrap.fraction.cutoff = .5, by.statistic = 'ES') {
  UseMethod("barplot_MSEA_boot")
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
    dplyr::filter(!is.na(.data$p_value), !is.na(.data$NES))


  ## filter terms with lower presence than bootstrap.fraction.cutoff
  if (any(c(object$consider_isomers, object$consider_isobars))){
    if("fraction" %in% colnames(enrichment_analysis)){
      enrichment_analysis = enrichment_analysis %>%
        dplyr::filter(.data$fraction > bootstrap.fraction.cutoff)
    }
    else{
      enrichment_analysis <-
        enrichment_analysis %>%
        dplyr::group_by(.data$Term) %>%
        dplyr::mutate(total_bootstrap = length(.data$bootstrap),
                      fraction = .data$total_bootstrap / object$enrichment_analysis$n) %>%
        dplyr::ungroup() %>%
        dplyr::filter(.data$fraction > bootstrap.fraction.cutoff)
    }
    if ("n" %nin% colnames(enrichment_analysis)){
      enrichment_analysis <-
        enrichment_analysis %>%
        dplyr::group_by(.data$Term) %>%
        dplyr::mutate(n = stats::median(.data$n, na.rm = T))
    }
  }

  enrichment_analysis <-
    enrichment_analysis %>%
    dplyr::group_by(.data$bootstrap) %>%
    dplyr::mutate(q.value = stats::p.adjust(.data$p_value, method = "fdr")) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$Term) %>%
    dplyr::mutate(p.value_combined = metap::sumlog(.data$p_value)[["p"]],
                  q.value_combined = metap::sumlog(.data$q.value)[["p"]]) %>%
    dplyr::arrange(.data$q.value_combined) %>%
    dplyr::ungroup()%>%
    dplyr::filter(.data$n > min.annotations,
                  .data$q.value_combined < q.value.cutoff,
                  .data$Term != "all")


  switch(by.statistic, 'ES' = {

    if (dim(enrichment_analysis)[1] < 1) {
      stop("Not enough enriched terms to visualize")
    }

    enrichment_plot <-
      enrichment_analysis %>%
      ggplot2::ggplot(ggplot2::aes(x = stats::reorder(.data$Term, .data$NES), y = .data$NES)) +
      ggplot2::coord_flip() +
      ggplot2::geom_bar(data = . %>% dplyr::group_by(.data$Term) %>%
                          dplyr::summarise(NES = stats::median(.data$NES, na.rm = TRUE),
                                           q.value_combined = dplyr::first(.data$q.value_combined)),
                        ggplot2::aes(fill = -log10(.data$q.value_combined)),
                        color = NA, stat = "identity", alpha = 0.5) +
      ggplot2::geom_jitter(size = 1, width = 0.1, color = "gray30") +
      ggplot2::scale_fill_viridis_c(option = "viridis",
                                    limits = c(0, max(-log10(enrichment_analysis$q.value_combined), na.rm = TRUE)),
                                    direction = 1) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
      ggplot2::labs(x = "", y = "Enrichment Score",
                    fill = "-log10 q value",
                    subtitle = paste(object$enrichment_analysis$comparison[2],
                                     "vs.",
                                     object$enrichment_analysis$comparison[1])) +
      ggpubr::theme_pubr() +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                     axis.title.x = ggplot2::element_text(face = "bold"))
  },
  'q.value' = {

    enrichment_analysis <-
      enrichment_analysis %>%
      dplyr::group_by(.data$Term) %>%
      dplyr::mutate(NES = stats::median(.data$NES, na.rm = T),
             up_down = factor(ifelse(sign(.data$NES)>0, "UP","DOWN"),
                              levels = c("UP","DOWN")))

    if(dim(enrichment_analysis)[1] < 1){
      stop("Not enough enriched terms to visualize")
    }

    enrichment_plot <-
      enrichment_analysis  %>%
      ggplot2::ggplot(ggplot2::aes(x = stats::reorder(.data$Term, dplyr::desc(-log(.data$q.value, base = 10))),
                          y = -log(.data$q.value, base = 10))) +
      ggplot2::coord_flip() +
      ggplot2::geom_bar(data = . %>%
                 dplyr::group_by(.data$Term) %>%
                 dplyr::summarise(up_down = .data$up_down[1],
                           q.value = stats::median(-log(.data$q.value, base = 10)),
                           q.value_clipped = ifelse(.data$q.value > 10, 10, .data$q.value)),
                 ggplot2::aes(y = .data$q.value, fill = .data$q.value_clipped), stat = "identity") +
      ggplot2::geom_jitter(size = 1, width = .1, color = "gray30") +
      ggplot2::facet_grid(.data$up_down~.,  space = "free", scales = "free") +
      ggplot2::scale_fill_gradient2(low = "gray", mid = "gray",
                           high = "red",midpoint = -log(0.05, base = 10),
                           limits = c(0,10)) +
      ggplot2::geom_hline(yintercept = c(0,-log(0.05,base = 10)), linetype = 3) +
      ggplot2::labs(x = "", y = "-Log10 (q-value)",
           fill = "-Log10 (q-value)",
           subtitle = paste(object$condition.y, "vs.", object$condition.x)) +
      ggpubr::theme_pubr(legend = "right") +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5), axis.title.x = ggplot2::element_text(face = "bold"))
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
#' data("example_MSEA_multicond")
#' p = plot_MSEA_Multi_cond(combined_MSEA_res = example_MSEA_multicond,
#'                          alpha_cutoff = 0.05)
#' }
#'
#' @export
plot_MSEA_Multi_cond = function(combined_MSEA_res,
                                alpha_cutoff = 0.05){

  missing_conds <- base::setdiff(c("condition.x", "condition.y"),
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
    dplyr::mutate(cond_label = paste0(.data$condition.y, "_",.data$condition.x)) %>%
    dplyr::select(.data$cond_label, .data$Term, .data$NES) %>%
    tidyr::pivot_wider(names_from = .data$cond_label, values_from = .data$NES) %>%
    tibble::column_to_rownames("Term") %>%
    as.matrix()

  label_mat = combined_MSEA_res %>%
    dplyr::mutate(cond_label = paste0(.data$condition.y, "_",.data$condition.x)) %>%
    dplyr::select(.data$cond_label, .data$Term, .data$signif) %>%
    tidyr::pivot_wider(names_from = .data$cond_label, values_from = .data$signif) %>%
    tibble::column_to_rownames("Term") %>%
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
#' data("example_ORA_markers")
#' bg = Load_background(mol_type = "Metabo",bg_type = "main_class",feature_type = "sf")
#' enrich_res = Run_simple_ORA(marker_list = example_ORA_markers,background = bg)
#' p = barplot_ORA_simple(enrich_res, q_val_cutoff = 0.2)
#' }
#'
#'
#' @export
barplot_ORA_simple = function(ORA_simple_res, q_val_cutoff = 0.05){
  plot_data = ORA_simple_res
  plot_data$Significance = ifelse(plot_data$score < -log10(q_val_cutoff), "No", "Yes")

  plot_data %>%
    ggplot2::ggplot(ggplot2::aes(x = stats::reorder(.data$Term, .data$score, max), y = .data$score, fill = .data$Significance)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_text(ggplot2::aes(label = as.character(.data$TP),hjust = 1), position = ggplot2::position_dodge(width = 0.5),
              color = "black", size = 6) +
    ggplot2::coord_flip() +
    ggpubr::theme_pubr() +
    ggplot2::ylab("-Log10 (q-value)") +
    ggplot2::xlab("") +
    ggplot2::geom_hline(yintercept = -log10(.data$q_val_cutoff), linetype = "dashed", lwd = 1) +
    ggplot2::theme(axis.title = ggplot2::element_text(size = 16),
          axis.text = ggplot2::element_text(size = 14))
}

#' Bar Plot for ORA Bootstrapped Results
#'
#' @description This function creates a bar plot for Over-Representation Analysis (ORA) bootstrapped results, including error bars to show the range of enrichment scores (OR) and indicating the significance based on q-values.
#'
#' @param ORA_boot_res A list containing the ORA bootstrapped results. The list should include `unfiltered_enrich_res` and `clean_enrich_res` data frames.
#' @param collapse_multi_cond A logical value indicating whether to collapse results across multiple conditions.
#'   Default is \code{FALSE}.
#' @return A ggplot object displaying the bar plot of ORA bootstrapped results.
#'
#' @details The function creates a bar plot where the terms are reordered by their median enrichment scores. Bars are filled with colors representing the combined q-values, and error bars show the minimum and maximum enrichment scores (ES_min and ES_max) across bootstraps. The number of true positives (n) is included in the term labels.
#'
#' @examples
#' \dontrun{
#' data("example_ORA_obj")
#' data("example_ORA_custom_universe")
#' input_scm = example_ORA_obj$scmatrix
#' conds = example_ORA_obj$conditions
#' cond_x = "U"
#' cond_y = "F"
#' ORA_obj <- initEnrichment(
#'   scmatrix = input_scm,
#'   conditions = conds,
#'   enrichment_type = "ORA",
#'   annot_db = "HMDB",
#'   consider_isomers = TRUE,
#'   consider_isobars = TRUE,
#'   polarization_mode = "positive",
#'   background_type = "sub_class",
#'   molecule_type = "Metabo",
#'   condition.x = cond_x,
#'   condition.y = cond_y
#' )
#' ORA_res <- Run_enrichment(
#'   object = ORA_obj,
#'   custom_universe = example_ORA_custom_universe,
#'   report_ambiguity_scores = TRUE,
#'   DE_LFC_cutoff = 0,
#'   min.pct.diff = 0
#' )
#' p = barplot_ORA_boot(ORA_boot_res = ORA_res,collapse_multi_cond = T)
#' }
#'
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
    dplyr::filter(.data$Term %fin% ORA_boot_res[["clean_enrich_res"]]$Term) %>%
    dplyr::group_by(.data$Term,.data$Condition) %>%
    dplyr::summarise(ES_max = max(.data$OR),
                     ES_min = min(.data$OR))


  data_viz = ORA_boot_res[["clean_enrich_res"]] %>%
    dplyr::left_join(boot_summary) %>%
    dplyr::mutate(n = ceiling(.data$n),
                  Term_label = paste0(.data$Term, " (", .data$n, ")")) %>%
    dplyr::select(.data$Term_label, .data$ES_median, .data$ES_max, .data$ES_min,
                  .data$p.value_combined, .data$q.value_combined,.data$Condition)

  min_qval = min(data_viz$q.value_combine)

  plot <- ggplot2::ggplot(data_viz, ggplot2::aes(x = stats::reorder(.data$Term_label, .data$ES_median),
                               y = .data$ES_median, fill = -log10(.data$q.value_combined))) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(), width = 0.7) +
    ggplot2::scale_y_continuous(n.breaks = 10) +
    ggplot2::geom_errorbar(ggplot2::aes_string(ymin = "ES_min", ymax = "ES_max"),
                  width = 0.2, position = ggplot2::position_dodge(0.7)) +
    ggplot2::scale_fill_gradient(low = "#B3CDE3", high = "#FBB4AE", name = "-Log10(q.value)") +
    ggpubr::theme_pubr() +
    ggplot2::coord_flip() +
    ggplot2::xlab("") +
    ggplot2::ylab("Fold Enrichment Score") +
    ggplot2::facet_wrap(~.data$Condition)
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
#' data("example_ORA_obj")
#' data("example_ORA_custom_universe")
#' input_scm = example_ORA_obj$scmatrix
#' conds = example_ORA_obj$conditions
#' cond_x = "U"
#' cond_y = "F"
#' ORA_obj <- initEnrichment(
#'   scmatrix = input_scm,
#'   conditions = conds,
#'   enrichment_type = "ORA",
#'   annot_db = "HMDB",
#'   consider_isomers = TRUE,
#'   consider_isobars = TRUE,
#'   polarization_mode = "positive",
#'   background_type = "sub_class",
#'   molecule_type = "Metabo",
#'   condition.x = cond_x,
#'   condition.y = cond_y
#' )
#' ORA_res <- Run_enrichment(
#'   object = ORA_obj,
#'   custom_universe = example_ORA_custom_universe,
#'   report_ambiguity_scores = TRUE,
#'   DE_LFC_cutoff = 0,
#'   min.pct.diff = 0
#' )
#' multi_cond_res = collapse_ORA_boot_multi_cond(ORA_boot_res_list = ORA_res)
#' p = dotplot_ORA(ORA_res = multi_cond_res$clean_enrich_res)
#' }
#'
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
    p = ggplot2::ggplot(plot_data, ggplot2::aes_string(x = "Condition",
                                     y = "Term", size = "size")) +
      ggplot2::geom_point(ggplot2::aes_string(color = "Score")) +
      ggplot2::scale_color_continuous(low = "#3B6FB6",
                             high = "#D41645", guide = ggplot2::guide_colorbar(reverse = TRUE)) +
      ggplot2::ylab(NULL) +
      # ggtitle("Trial") +
      ggpubr::theme_pubr() +
      ggplot2::scale_size_continuous(range = c(4, 9)) +
      ggplot2::scale_y_discrete(labels = label_func) +
      ggpubr::rotate_x_text(angle = 45) +
      ggplot2::guides(size = ggplot2::guide_legend(order = 1, title = 'Term/query overlap'),
             color = ggplot2::guide_colorbar(order = 2, title = color_legend_title)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 14),
            axis.text.y = ggplot2::element_text(size = 14),
            axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_text(size = 16),
            strip.text = ggplot2::element_text(size = 14))
  }
  else{
    p = ggplot2::ggplot(data = plot_data, ggplot2::aes(x = .data$Enrichment_Score,
                                     y = stats::reorder(.data$Term, .data$Enrichment_Score))) +
      ggplot2::geom_point(ggplot2::aes(size = .data$size,color = .data$Score)) +
      ggplot2::scale_size(range = c(3, 6)) +
      ggplot2::scale_color_gradient(low = cols[1], high = cols[2]) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank()) +
      ggplot2::guides(size = ggplot2::guide_legend(title = 'Term/query overlap'),
             color = ggplot2::guide_colorbar(title = color_legend_title)) +
      ggplot2::labs(
        x = 'Enrichment Score',
        y = ''
      ) +
      ggpubr::theme_pubr() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0,hjust = 1))
  }

  if (!is.null(facet_by)){
    facet_col_sym = ggplot2::sym(facet_by)
    p = p +
      ggplot2::facet_wrap(ggplot2::vars(!!ggplot2::enquo(facet_col_sym))) +
      ggplot2::theme(
        panel.spacing = ggplot2::unit(x = 1, units = "lines"),
        strip.background = ggplot2::element_blank())
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
#' @param condition Optional. A character string specifying a particular condition to filter the enrichment results. If \code{NULL} (the default), results from all conditions will be included.
#' @return A ggplot object displaying the ridge plot of intersection sizes for the specified terms of interest.
#'
#' @details The function filters the enrichment results to include only the specified terms of interest and creates a ridge plot showing the distribution of intersection sizes (TP). The terms are labeled with their names and the fraction of bootstraps in which they appear, and the plot includes quantile lines to indicate distribution quartiles.
#'
#' @examples
#' \dontrun{
#' data("example_ORA_obj")
#' data("example_ORA_custom_universe")
#' input_scm = example_ORA_obj$scmatrix
#' conds = example_ORA_obj$conditions
#' cond_x = "U"
#' cond_y = "F"
#' ORA_obj <- initEnrichment(
#'   scmatrix = input_scm,
#'   conditions = conds,
#'   enrichment_type = "ORA",
#'   annot_db = "HMDB",
#'   consider_isomers = TRUE,
#'   consider_isobars = TRUE,
#'   polarization_mode = "positive",
#'   background_type = "sub_class",
#'   molecule_type = "Metabo",
#'   condition.x = cond_x,
#'   condition.y = cond_y
#' )
#' ORA_res <- Run_enrichment(
#'   object = ORA_obj,
#'   custom_universe = example_ORA_custom_universe,
#'   report_ambiguity_scores = TRUE,
#'   DE_LFC_cutoff = 0,
#'   min.pct.diff = 0
#' )
#' multi_cond_res = collapse_ORA_boot_multi_cond(ORA_boot_res_list = ORA_res)
#' p4 = ridge_bootstraps(enrich_res = multi_cond_res$unfiltered_enrich_res,
#'                      terms_of_interest = c(multi_cond_res$clean_enrich_res$Term),
#'                      condition = "upregulated")
#' }
#'
#' @export
ridge_bootstraps = function(enrich_res, terms_of_interest, condition = NULL){


  plot_data = enrich_res[enrich_res$Term %in% terms_of_interest,]
  colnames(plot_data)[which(colnames(plot_data) == "n")] = "TP"

  plot_data$Term = paste0(plot_data$Term, " (", plot_data$fraction * 100,
                          ")")

  if (!is.null(condition)){
    plot_data = plot_data[plot_data$Condition == condition,]
  }

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$TP, y = .data$Term,
                                 fill = ggplot2::after_stat(density))) +
    ggridges::stat_density_ridges(geom = "density_ridges_gradient",calc_ecdf = T,
                                  scale = 0.95) +
    ggplot2::scale_fill_viridis_c(name = "Density", direction = 1) +
    ggplot2::scale_x_continuous(n.breaks = 15, limits = c(0,max(plot_data$TP) +2)) +
    ggplot2::ylab("") +
    ggplot2::xlab("Intersection size") +
    ggpubr::theme_pubr()
}

#' Plot intensity distribution of metabolite across conditions
#'
#' @param obj A S2IsoMEr object after enrichment analysis.
#' @param metabolite A character indicating metabolite name
#' @param conds_of_interest Optional. A character vector specifying which conditions to include in the plot. If `NULL` (the default), all conditions are included.
#'
#' @return  A `ggplot2` object displaying the density plot of metabolite intensities across conditions.
#' @examples
#' \dontrun{
#' data("example_ORA_obj")
#' data("example_ORA_custom_universe")
#' input_scm = example_ORA_obj$scmatrix
#' conds = example_ORA_obj$conditions
#' cond_x = "U"
#' cond_y = "F"
#' ORA_obj <- initEnrichment(
#'   scmatrix = input_scm,
#'   conditions = conds,
#'   enrichment_type = "ORA",
#'   annot_db = "HMDB",
#'   consider_isomers = TRUE,
#'   consider_isobars = TRUE,
#'   polarization_mode = "positive",
#'   background_type = "sub_class",
#'   molecule_type = "Metabo",
#'   condition.x = cond_x,
#'   condition.y = cond_y
#' )
#' ORA_res <- Run_enrichment(
#'   object = ORA_obj,
#'   custom_universe = example_ORA_custom_universe,
#'   report_ambiguity_scores = TRUE,
#'   DE_LFC_cutoff = 0,
#'   min.pct.diff = 0
#' )
#' TP_ions = get_TP_markers_per_Term(ORA_boot_df = ORA_res$downregulated$unfiltered_enrich_res,
#'                                  term_of_interest = "Glycerophosphocholines")
#' TP_ions = map_TP_markers_to_ions(markers = TP_ions,scm_ions = rownames(input_scm))
#' p = compare_metabo_distr(ORA_obj, metabolite = TP_ions[5],
#'                          conds_of_interest = c(cond_x, cond_y))
#' }
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
      dplyr::filter(.data$condition %in% conds_of_interest)
  }

  desc_stats <- data_viz %>%
    dplyr::group_by(.data$condition) %>%
    dplyr::summarise(
      mean = mean(.data$intens),
      median = stats::median(.data$intens)
    )

  plot <- ggplot2::ggplot(data_viz, ggplot2::aes(x = .data$intens, fill = .data$condition)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::geom_density(alpha = 0.4) +
    ggplot2::geom_vline(data = desc_stats, ggplot2::aes(xintercept = .data$mean, color = .data$condition, linetype = "Mean"), size = 1) +
    ggplot2::geom_vline(data = desc_stats, ggplot2::aes(xintercept = .data$median, color = .data$condition, linetype = "Median"), size = 1) +
    ggplot2::scale_linetype_manual(name = "Statistic", values = c("Mean" = "dashed", "Median" = "solid")) +
    ggplot2::labs(title = metabolite,
         x = "Intensity(Log10)",
         y = "Density",
         fill = "Condition",
         color = "Condition",
         linetype = "Statistic") +
    ggpubr::theme_pubr() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::scale_x_continuous(n.breaks = 10) +
    ggplot2::facet_wrap(~condition, scales = "free_y")

  return(plot)

}

