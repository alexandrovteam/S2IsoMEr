decouple_ORA_wrapper = function(marker_list,term_list, n_background, seed = 42){
  # n_up = nrow(input_mat)
  # n_bottom = 0

  # withr::with_seed(seed, {
  #   targets <- decoupleR:::.ora_slice_targets(input_mat, n_up, n_bottom, with_ties)
  # })

  ORA_res = decoupleR:::.ora_analysis(regulons = term_list, targets = marker_list,
                                      n_background = n_background)
  return(ORA_res)
}

hyper_geom_enrich = function(query,term_list, universe = NULL,
                             fisher_alternative = "greater",
                             min_path_size = 3){
  if (!is.null(universe)){
    all_bg_mols = universe
  }
  else{
    all_bg_mols = unlist(term_list) %>% unique()
  }

  query = query[which(query %in% all_bg_mols)]

  final_bg = lapply(term_list, function(x){
    x[x %in% all_bg_mols]
  })
  final_bg = final_bg[lengths(final_bg) >= min_path_size]

  final_res = sapply(final_bg, function(x){
    TP = length(intersect(query, x))
    FP = length(which(query %nin% x))
    FN = length(which(x %nin% query))

    TN = all_bg_mols[which(all_bg_mols %nin% x)]
    TN = length(TN[which(TN %nin% query)])

    # fisher_mat = matrix(c(TP, FP, FN, TN), nrow = 2, ncol = 2, byrow = T)
    # fisher_res = hypergea::hypergeom.test(fisher_mat, alternative = fisher_alternative)

    term_res = data.frame(TP = TP, FP = FP, FN = FN , TN = TN)
    term_res
  },USE.NAMES = T,simplify = F)

  terms = names(final_res)
  final_res = final_res %>% dplyr::bind_rows()
  final_res$term = terms

  final_res = final_res[,c(5,1:4)]

  return(final_res)
}

get_metabo_iso = function(sf_vec, consider_isobars = T,
                          polarization_mode = NA, mass_range_ppm = 3,
                          only_HMDB = F){

  annotation_formulas_adduct <- gsub("\\+|\\-",".",sf_vec)
  annotation_formulas <- gsub("\\..+$","",annotation_formulas_adduct)
  annotation_adduct <- gsub("^.+\\.","",annotation_formulas_adduct)
  annotation_list <-
    lapply(annotation_formulas, function(annotation_formula_i){
      if (only_HMDB){
        res = metaspace_databases$name[which(metaspace_databases$formula == annotation_formula_i &
                                               metaspace_databases$db == "HMDB")]
      }
      else{
        res = metaspace_databases$name[metaspace_databases$formula == annotation_formula_i]
      }
      res
    })
  is_adduct = any(c(paste0("pos", annotation_adduct),
                    paste0("neg", annotation_adduct)) %in% colnames(metabo_exact_masses))
  if (consider_isobars &
      is_adduct &
      polarization_mode %in% c("negative", "positive")){
    switch(polarization_mode,
           'positive' = {
             col_name <- paste0("pos",annotation_adduct)

             exact_masses_slim <- metabo_exact_masses[,c(1,which(grepl("pos",colnames(metabo_exact_masses))))]
             colnames(exact_masses_slim) <- gsub("^pos","",colnames(exact_masses_slim))
             exact_masses_slim <- exact_masses_slim %>% tidyr::pivot_longer(cols = -1, values_to = "mass", names_to = "adduct")
             exact_masses_slim$formula_adduct <- paste0(exact_masses_slim$formula,".",exact_masses_slim$adduct)
           },
           'negative' = {
             col_name <- paste0("neg",annotation_adduct)

             exact_masses_slim <- metabo_exact_masses[,c(1,which(grepl("neg",colnames(metabo_exact_masses))))]
             colnames(exact_masses_slim) <- gsub("^neg","",colnames(exact_masses_slim))
             exact_masses_slim <- exact_masses_slim %>% tidyr::pivot_longer(cols = -1, values_to = "mass", names_to = "adduct")
             exact_masses_slim$formula_adduct <- paste0(exact_masses_slim$formula,".",exact_masses_slim$adduct)
           })
    exact_masses_annotations <-
      mapply(annotation_formulas_i = annotation_formulas,
             col_name_i = col_name,
             function(annotation_formulas_i, col_name_i) {
               mass <-
                 metabo_exact_masses[metabo_exact_masses$formula == annotation_formulas_i, col_name_i]
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


  }
  else{
    isobars_list = NULL
  }

  if (!is.null(isobars_list)){
    annotation_list <-
      lapply(seq_along(annotation_formulas), function(i){
        c(isomer = annotation_list[[i]],
          isobar = metaspace_databases$name[metaspace_databases$formula %in% gsub("\\..+$","",isobars_list[[i]])])

      })
  }
  names(annotation_list) = sf_vec

  return(annotation_list)
}

metabo_bootstrap = function(annot_list, annot_weights = NULL,
                            n_bootstraps = 50,sample_core_metab_only = T){
  if (!is.null(annot_weights)) {
    if (length(annot_weights) != length(annot_list)) {
      stop("annotations and annotation.weights do not have the same length")
    } else {
      ## extra check
      if (!all(sapply(annot_weights, length) == sapply(annot_list, length))) {
        stop("annotations and annotation.weights do not have the same length ")
      }
    }
  }
  if (sample_core_metab_only){
    annot_list <- sapply(annot_list, function(i){
      i[i %in% core_metab$name]
    }, simplify = F)
  }

  bootstrapped_sublist <- pbapply::pblapply(seq(n_bootstraps),function(n_i){
    sapply(seq(length(annot_list)),function(row_number_i){
      molecules_to_sample <- annot_list[[row_number_i]]
      if(length(molecules_to_sample)==0){
        molecules_to_sample <- names(annot_list)[[row_number_i]]
      }

      if(!is.null(annot_weights)){
        weights_to_sample <- annot_weights[[row_number_i]]
      } else {
        weights_to_sample <- rep(x = 1, times = length(molecules_to_sample))
      }
      if(length(molecules_to_sample)!=1){
        sample(x = molecules_to_sample,size = 1,prob = weights_to_sample)
      }
      else {
        molecules_to_sample
      }

    })
  })
  return(bootstrapped_sublist)
}

simplify_hypergeom_bootstrap = function(bootstrap_list,term_list,universe = NULL,
                                        core_metabo_only = F,
                                        boot_fract_cutoff = 0.5,
                                        min_annot = 3, q.val_cutoff = 0.2,
                                        selected_terms = NULL,
                                        alpha_cutoff = 0.05){

  if (!is.null(universe)){
    pathway_list_slim <- sapply(term_list, function(i){
      i[i %in% universe]
    }, simplify = F)
    term_list <- pathway_list_slim[sapply(pathway_list_slim, length) > 0]
  }
  if (!is.null(selected_terms)){
    term_list = term_list[which(names(term_list) %in% selected_terms)]
  }

  enrich_res = decouple_ORA_wrapper(marker_list = bootstrap_list, term_list = term_list,
                                    n_background = ifelse(is.null(universe),
                                                          length(unique(unlist(term_list))),
                                                          length(unique(unlist(universe)))))
  colnames(enrich_res)[which(colnames(enrich_res) == "source")] = "term"
  colnames(enrich_res)[which(colnames(enrich_res) == "condition")] = "bootstrap"


  boot_conting_res = pbapply::pblapply(seq(length(bootstrap_list)),function(n_i){
    enrich_res = hyper_geom_enrich(query = unique(bootstrap_list[[n_i]]),
                                   term_list = term_list, universe = universe,
                                   core_metabo_only = core_metabo_only)
    enrich_res$bootstrap = n_i
    #enrich_res = enrich_res %>% dplyr::filter(!is.na(OR), !is.na(pval))
    if(nrow(enrich_res) == 0){
      return(NULL)
    }
    else{
      return(enrich_res)
    }
  })
  boot_conting_res = boot_conting_res %>% dplyr::bind_rows()

  observed = boot_conting_res$TP / (boot_conting_res$TP + boot_conting_res$FP)
  expected = (boot_conting_res$TP + boot_conting_res$FN) / (boot_conting_res$TP + boot_conting_res$FP +
                                                              boot_conting_res$FN + boot_conting_res$TN)
  boot_conting_res$OR = observed / expected


  boot_enrich_res = boot_conting_res %>% dplyr::left_join(enrich_res, by = c("term","bootstrap"))

  boot_enrich_res = boot_enrich_res %>%
    mutate(padj = p.adjust(p_value, "BH")) %>%
    dplyr::filter(TP >= min_annot, p_value < alpha_cutoff)

  boot_enrich_res <- boot_enrich_res %>%
    dplyr::group_by(term) %>%
    dplyr::mutate(fraction = length(term) / length(bootstrap_list)) %>%
    dplyr::filter(fraction > boot_fract_cutoff)

  final_enrich_res <- boot_enrich_res %>% group_by(bootstrap) %>%
    dplyr::mutate(q.value = p.adjust(p_value, method = "fdr"))  %>%
    dplyr::group_by(term) %>%
    dplyr::summarise(n = median(TP, na.rm = T),
                     ES_median = median(OR, na.rm = T),
                     ES_sd = sd(OR, na.rm = T),
                     p.value_median = median(p_value, na.rm = T),
                     p.value_sd = sd(p_value, na.rm = T),
                     q.value_median = median(q.value, na.rm = T),
                     q.value_sd = sd(q.value, na.rm = T),
                     fraction.bootstrap.presence = median(fraction, na.rm = T)) %>%
    dplyr::arrange(q.value_median)

  final_enrich_res = final_enrich_res %>%
    dplyr::filter(n > min_annot, q.value_median < q.val_cutoff, term != "") %>%
    ungroup() %>% as.data.frame()

  return(list("unfiltered_enrich_res" = boot_enrich_res,
              "clean_enrich_res" = final_enrich_res))

}

Run_simple_ORA = function(marker_list, background, custom_universe = NULL,
                          alpha_cutoff = 0.05, min_intersection = 3){
  if (!is.list(marker_list)){
    q = sub("[-+].*","", marker_list) %>% unique()
    ORA_conting = hyper_geom_enrich(query = q, term_list = background,
                               universe = custom_universe)
    ORA_res = decouple_ORA_wrapper(marker_list = list("Condition" = q),
                                   term_list = background,
                                   n_background = ifelse(is.null(custom_universe),
                                                         length(unique(unlist(term_list))),
                                                         length(unique(unlist(custom_universe)))))
    ORA_final = ORA_conting %>% dplyr::left_join(ORA_res, by = c("term" = "source"))
    ORA_final = ORA_final %>%
      mutate(padj = p.adjust(p_value, "BH")) %>%
      dplyr::filter(TP >= min_intersection, p_value < alpha_cutoff)
    return(ORA_final)
  }
  else{
    ORA_res = decouple_ORA_wrapper(marker_list = marker_list,
                                   term_list = background,
                                   n_background = ifelse(is.null(custom_universe),
                                                         length(unique(unlist(term_list))),
                                                         length(unique(unlist(custom_universe)))))
    ORA_conting_all = list()
    for (grp in 1:length(marker_list)){
      q = marker_list[[grp]]
      q = sub("[-+].*","", q) %>% unique()
      ORA_conting = hyper_geom_enrich(query = q, term_list = background,
                                 universe = custom_universe)
      ORA_conting_all[[names(marker_list)[grp]]] = ORA_conting
    }
    ORA_conting_all = ORA_conting_all %>% dplyr::bind_rows(.id = "Group")
    ORA_final = ORA_conting_all %>% dplyr::left_join(ORA_res, by = c("term" = "source"))
    ORA_final = ORA_final %>%
      mutate(padj = p.adjust(p_value, "BH")) %>%
      dplyr::filter(TP >= min_intersection, p_value < alpha_cutoff)

    return(ORA_final)
  }
}


Run_bootstrap_ORA = function(marker_list, background, custom_universe = NULL,
                             alpha_cutoff = 0.05, min_intersection = 3,
                             consider_isobars = T,polarization_mode = NA, mass_range_ppm = 3,
                             only_HMDB = F,annot_weights = NULL,
                             n_bootstraps = 50,sample_core_metab_only = T,
                             boot_fract_cutoff = 0.5,q.val_cutoff = 0.2,
                             selected_terms = NULL){

  if (!is.list(marker_list)){
    q = sub("[-+].*","", marker_list) %>% unique()
    marker_list = list("query" = q)
  }
  ORA_boot_all_grps = list()
  for (grp in 1:length(marker_list)){
    q = marker_list[[grp]]
    q = sub("[-+].*","", q) %>% unique()

    message(paste0("\n", "Getting Isomers and Isobars", "\n"))

    iso_list = get_metabo_iso(sf_vec = q, consider_isobars = consider_isobars, polarization_mode = polarization_mode,
                              mass_range_ppm = mass_range_ppm, only_HMDB = only_HMDB)

    boot_list = metabo_bootstrap(annot_list = iso_list, annot_weights = annot_weights,
                                 n_bootstraps = n_bootstraps, sample_core_metab_only = sample_core_metab_only)

    final_res = simplify_hypergeom_bootstrap(bootstrap_list = boot_list,
                                             term_list = background,
                                             core_metabo_only = sample_core_metab_only,universe = custom_universe,
                                             boot_fract_cutoff = boot_fract_cutoff, min_annot = min_intersection,
                                             q.val_cutoff = q.val_cutoff,selected_terms = selected_terms,
                                             alpha_cutoff = alpha_cutoff)
    ORA_boot_all_grps[[grp]] = final_res
  }
  return(ORA_boot_all_grps)
}


