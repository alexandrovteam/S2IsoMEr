

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
                                    universe = univ, pass_adjust = pass_adjust)

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
#' @export
Run_simple_ORA = function(marker_list, background, custom_universe = NULL,
                          alpha_cutoff = 0.05, min_intersection = 3){

  if (!is.list(marker_list)){
    q = sub("[-+].*","", marker_list) %>% unique()
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

  if (!is.list(marker_list)){
    q = sub("[-+].*","", marker_list) %>% unique()
    marker_list = list("query" = q)
  }
  ORA_boot_all_grps = list()
  for (grp in 1:length(marker_list)){
    q = marker_list[[grp]]
    q = sub("[-+].*","", q) %>% unique()

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
      ORA_boot_all_grps[[grp]] = c(final_res, ambig_score)
    }
    else{
      ORA_boot_all_grps[[grp]] = final_res
    }

  }
  return(ORA_boot_all_grps)
}
