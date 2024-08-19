#' Wrapper function for Simple Over-Representation Analysis (ORA)
#'
#' @description This function performs a simple over-representation analysis (ORA) on a given list of metabolites and background.
#'
#' @param marker_list A vector or list of marker sets to analyze.
#' @param background A list of named terms where each term is a character vector of metabolites.
#' @param custom_universe An optional vector specifying a custom universe of terms. If not provided all metabolites in background will be used as universe.
#' @param alpha_cutoff A numeric value indicating the alpha cutoff for significance.
#' @param min_intersection An integer specifying the minimum intersection required between input list and a given term in background.
#' @return A dataframe containing the ORA results for each group of markers.
#' @examples
#' \dontrun{
#' data("example_ORA_markers")
#' bg = Load_background(mol_type = "Metabo",bg_type = "main_class",feature_type = "sf")
#' enrich_res = Run_simple_ORA(marker_list = example_ORA_markers,background = bg)
#' }
#' @export
Run_simple_ORA = function(marker_list, background, custom_universe = NULL,
                          alpha_cutoff = 0.05, min_intersection = 3){

  if (!is.list(marker_list)){
    q = sub("[.+-].*","", marker_list) %>% unique()
    # ORA_conting = hyper_geom_enrich(query = q, term_list = background,
    #                            universe = custom_universe)
    marker_list  = list("Condition" = q)
  }
  else{
    marker_list = lapply(marker_list, function(x){
      sub("[-+].*","", x) %>% unique()
    })
  }

  if (!is.null(custom_universe)){
    pathway_list_slim <- sapply(background, function(i){
      i[i %in% custom_universe]
    }, simplify = F)
    background <- pathway_list_slim[sapply(pathway_list_slim, length) > 0]
    univ = custom_universe
  }
  else{
    pathway_list_slim <- sapply(background, function(i){
      length(intersect(i, unlist(marker_list)))
      }, simplify = T)
    background <- background[which(pathway_list_slim != 0)]
    univ = unlist(background) %>% unique()
  }

  background = background[which(names(background) != "all")]

  ORA_res = decouple_ORA_wrapper(marker_list = marker_list,
                                 term_list = background,
                                 universe = univ,ORA_boot = F)

  ORA_conting = ORA_res[[2]]
  ORA_res = ORA_res[[1]]
  ORA_final = ORA_conting %>% dplyr::left_join(ORA_res)
  ORA_final = ORA_final %>%
    mutate(q.value = p.adjust(p_value, "BH")) %>%
    dplyr::filter(TP >= min_intersection, p_value < alpha_cutoff)

  colnames(ORA_final)[which(colnames(ORA_final) == "source")] = "Term"

  return(ORA_final)
}

#' Wrapper function for Bootstrap Over-Representation Analysis (ORA)
#'
#' @description This function performs a bootstrap-based over-representation analysis (ORA) on a given list of metabolites and background.
#'
#' @param marker_list A vector or list of marker sets to analyze.
#' @param background A list of named terms where each term is a character vector of metabolites.
#' @param custom_universe An optional vector specifying a custom universe of terms. If not provided all metabolites in background will be used as universe.
#' @param alpha_cutoff A numeric value indicating the alpha cutoff for significance.
#' @param min_intersection An integer specifying the minimum intersection required between input list and a given term in background.
#' @param consider_isobars A logical indicating whether to consider isobars.
#' @param polarization_mode A parameter for polarization mode. If provided, it will be used to consider isomers and isobars.
#' @param mass_range_ppm A numeric value for mass range in parts per million (ppm).
#' @param annot_db A character string specifying the annotation database(s) used for annotation. Check  \code{\link{build_iso_bg}} for more information.
#' @param annot_custom_db An optional custom annotation database. If provided, it will be used alongside or instead of the specified annotation database.
#' @param use_LION A logical indicating whether to use LION ontology for select isomers/isobars for background.
#' @param endogenous_only A logical indicating whether to consider only endogenous compounds. Applies only for HMDB and CoreMetabolome annotation databases.
#' @param pathway_assoc_only A logical indicating whether to consider only pathway-associated compounds. Applies only for HMDB and CoreMetabolome annotation databases.
#' @param remove_expected_predicted A logical indicating whether to remove expected and predicted annotations. Applies only for HMDB and CoreMetabolome annotation databases.
#' @param annot_list An optional list of annotations. If provided, it will be used instead of generating isomers from built-in formula to molecule mapping.
#' @param annot_weights An optional list of annotation weights. If provided, it will be used in the ambiguity score calculation and bootstrap analysis.
#' @param n_bootstraps An integer specifying the number of bootstrap iterations.
#' @param boot_fract_cutoff A numeric value specifying the minimum bootstrap fraction cutoff required. Default to 0.5 .
#' @param q.val_cutoff A numeric value specifying the q-value cutoff for significance.
#' @param selected_terms An optional vector of selected terms to focus on. If provided, the analysis will be restricted to these terms.
#' @param adjust_contingency A logical indicating whether to adjust the contingency table to account for isomeric ambiguity.
#' @param report_ambiguity_scores A logical indicating whether to report ambiguity scores. If TRUE, ambiguity scores will be included in the results.
#' @return A list containing the ORA results for each group of markers. Both filtered and per-bootstrap results are provided.
#' @examples
#' \dontrun{
#' data("example_ORA_markers")
#' data("example_ORA_obj")
#' object = example_ORA_obj
#' enrich_res <- Run_bootstrap_ORA(
#'   marker_list = example_ORA_markers,
#'   background = object$pathway_list,
#'   polarization_mode = object$polarization_mode,
#'   mass_range_ppm = object$mass_range_ppm,
#'   annot_db = object$Annotation_database,
#'   annot_custom_db = object$Custom_database,
#'   use_LION = ifelse(stringr::str_detect(object$background_name, "LION"), TRUE, FALSE),
#'   endogenous_only = object$endogenous_only,
#'   pathway_assoc_only = object$pathway_assoc_only,
#'   remove_expected_predicted = object$remove_expected_predicted,
#'   annot_weights = object$annotation.weights,
#'   consider_isobars = object$consider_isobars,
#'   annot_list = object$annotations
#' )
#' }
#' @export
Run_bootstrap_ORA = function(marker_list, background, custom_universe = NULL,
                             alpha_cutoff = 0.05, min_intersection = 3,
                             consider_isobars = T,polarization_mode = NA, mass_range_ppm = 3,
                             annot_db = "HMDB", annot_custom_db = NULL,
                             use_LION = F, endogenous_only = T,
                             pathway_assoc_only = F,
                             remove_expected_predicted = T,
                             annot_list = NULL,
                             annot_weights = NULL,
                             n_bootstraps = 50,
                             boot_fract_cutoff = 0.5,q.val_cutoff = 0.2,
                             selected_terms = NULL,
                             adjust_contingency = T,
                             report_ambiguity_scores = F){

  if(!is.null(custom_universe)){
    univ_iso = get_metabo_iso(sf_vec = custom_universe, consider_isobars = consider_isobars, polarization_mode = polarization_mode,
                              mass_range_ppm = mass_range_ppm,annot_db = annot_db,
                              annot_custom_db = annot_custom_db, use_LION = use_LION, endogenous_only = endogenous_only,
                              pathway_assoc_only = pathway_assoc_only, remove_expected_predicted = remove_expected_predicted)
    custom_universe = univ_iso %>% unlist() %>% unique()
  }
  else{
    univ_iso = NULL
  }

  if (!is.list(marker_list)){
    # q = sub("[-+].*","", marker_list) %>% unique()
    q = gsub("\\+|\\-",".",marker_list)
    marker_list = list("query" = q)
  }
  ORA_boot_all_grps = list()
  for (grp in 1:length(marker_list)){
    q = marker_list[[grp]]
    q = gsub("\\+|\\-",".",q)
    # q = sub("[-+].*","", q) %>% unique()

    message(paste0("\n", "Getting Isomers and Isobars", "\n"))

    if(!is.null(annot_list)){
      q = q[which(q %in% names(annot_list))]
      iso_list = annot_list
      iso_list = iso_list[q]
    }
    else{
      iso_list = get_metabo_iso(sf_vec = q, consider_isobars = consider_isobars, polarization_mode = polarization_mode,
                                mass_range_ppm = mass_range_ppm,annot_db = annot_db,
                                annot_custom_db = annot_custom_db, use_LION = use_LION, endogenous_only = endogenous_only,
                                pathway_assoc_only = pathway_assoc_only, remove_expected_predicted = remove_expected_predicted)
    }

    if(report_ambiguity_scores){
      if(!is.null(univ_iso)){
        univ_ambig = calc_ambiguity(input_iso_list = univ_iso, weights = NULL)
      }
      else{
        univ_ambig = NULL
      }
      iso_ambig = calc_ambiguity(input_iso_list = iso_list, weights = annot_weights)
      ambig_scores = list("Query ambiguity" = iso_ambig,
                          "Universe ambiguity" = univ_ambig)
    }
    else{
      ambig_scores = NULL
    }

    boot_list = metabo_bootstrap(annot_list = iso_list, annot_weights = annot_weights,
                                 n_bootstraps = n_bootstraps)

    final_res = simplify_hypergeom_bootstrap(bootstrap_list = boot_list,
                                             term_list = background,
                                             universe = custom_universe,
                                             boot_fract_cutoff = boot_fract_cutoff, min_intersection = min_intersection,
                                             q.val_cutoff = q.val_cutoff,selected_terms = selected_terms,
                                             alpha_cutoff = alpha_cutoff, pass_adjust = !adjust_contingency)
    if (report_ambiguity_scores){
      ORA_boot_all_grps[[names(marker_list)[grp]]] = c(final_res, ambig_scores)
    }
    else{
      ORA_boot_all_grps[[names(marker_list)[grp]]] = final_res
    }

  }
  return(ORA_boot_all_grps)
}

simplify_hypergeom_bootstrap = function(bootstrap_list,term_list,universe = NULL,
                                        boot_fract_cutoff = 0.5,
                                        min_intersection = 3, q.val_cutoff = 0.2,
                                        selected_terms = NULL,
                                        alpha_cutoff = 0.05,
                                        pass_adjust = F){

  if (!is.null(universe)){
    pathway_list_slim <- sapply(term_list, function(i){
      i[i %in% universe]
    }, simplify = F)
    term_list <- pathway_list_slim[sapply(pathway_list_slim, length) > 0]
    univ = universe
  }
  else{
    univ = unlist(term_list) %>% unique()
  }
  if (!is.null(selected_terms)){
    term_list = term_list[which(names(term_list) %in% selected_terms)]
  }

  enrich_res = decouple_ORA_wrapper(marker_list = bootstrap_list, term_list = term_list,
                                    universe = univ, pass_adjust = pass_adjust,ORA_boot = T)

  boot_conting_res = enrich_res[[2]]
  enrich_res = enrich_res[[1]]

  colnames(enrich_res)[which(colnames(enrich_res) == "source")] = "Term"
  colnames(enrich_res)[which(colnames(enrich_res) == "condition")] = "bootstrap"

  colnames(boot_conting_res)[which(colnames(boot_conting_res) == "source")] = "Term"
  colnames(boot_conting_res)[which(colnames(boot_conting_res) == "condition")] = "bootstrap"

  observed = boot_conting_res$TP / (boot_conting_res$TP + boot_conting_res$FP)
  expected = (boot_conting_res$TP + boot_conting_res$FN) / (boot_conting_res$TP + boot_conting_res$FP +
                                                              boot_conting_res$FN + boot_conting_res$TN)
  boot_conting_res$OR = observed / expected


  boot_enrich_res = boot_conting_res %>%
    dplyr::left_join(enrich_res, by = c("Term","bootstrap")) %>%
    dplyr::mutate(padj = p.adjust(p_value, "BH")) %>%
    dplyr::group_by(Term) %>%
    dplyr::mutate(fraction = length(Term) / length(bootstrap_list)) %>%
    dplyr::ungroup()


  final_enrich_res = boot_enrich_res %>%
    dplyr::filter(fraction > boot_fract_cutoff) %>%
    dplyr::group_by(bootstrap) %>%
    dplyr::mutate(q.value = p.adjust(p_value, method = "fdr"))  %>%
    dplyr::group_by(Term) %>%
    dplyr::summarise(n = median(TP, na.rm = T),
                     ES_median = median(OR, na.rm = T),
                     ES_sd = sd(OR, na.rm = T),
                     p.value_combined = metap::sumlog(p_value)[["p"]],
                     q.value_combined = metap::sumlog(q.value)[["p"]],
                     fraction.bootstrap.presence = median(fraction, na.rm = T)) %>%
    dplyr::arrange(q.value_combined) %>%
    dplyr::filter(n >= min_intersection,
                  q.value_combined < q.val_cutoff,
                  p.value_combined < alpha_cutoff,
                  Term != "") %>%
    dplyr::ungroup() %>%
    as.data.frame()

  return(list("unfiltered_enrich_res" = boot_enrich_res,
              "clean_enrich_res" = final_enrich_res))

}
