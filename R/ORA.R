#TODO make sure all variable names are the same between simple and bootstraps


.ora_analysis <- function(regulons, targets, universe, ...) {

  # NSE vs. R CMD check workaround
  p.value <- NULL

  message("\nRunning ORA analysis ... \n")

  n = length(names(regulons)) * length(names(targets))

  pb = progress::progress_bar$new(total = n, incomplete = " ")

  tidyr::expand_grid(source = names(regulons), condition = names(targets)) %>%
    dplyr::rowwise(source, condition) %>%
    dplyr::summarise(.ora_fisher_exact_test(
      dat = adjust_conting_iso(
        expected = regulons[[source]],
        observed = targets[[condition]],
        universe = universe
      ),
      pbar = pb,
      ...
    ),
    .groups = "drop"
    ) %>%
    dplyr::select(source, condition,
           p_value = p.value, everything()
    ) %>%
    dplyr::mutate(score = -log10(p_value)) %>%
    tibble::add_column(statistic = "ora", .before = 1) %>%
    dplyr::select(statistic, source, condition, score, p_value,
                  TP, FP, FN, TN)
}
.ora_fisher_exact_test <- function(dat,pbar, ...) {
  pbar$tick()
  conting = ora_conting_decoupleR(dat, as_matrix = F) %>% as_tibble()
  rlang::exec(
    .fn = stats::fisher.test,
    x = matrix(as.integer(conting), nrow = 2, ncol = 2),
    y = NULL,
    alternative='greater',
    !!!list(...)
  ) %>%
    broom::glance() %>%
    cbind(conting)
}
ora_conting_decoupleR = function(dat, as_matrix = T) {

  observed = dat[["obs"]]
  expected = dat[["exp"]]
  n_background = dat[["n_bg"]]

  true_positive <- intersect(observed, expected) %>% length()
  false_negative <- setdiff(expected, observed) %>% length()
  false_positive <- setdiff(observed, expected) %>% length()
  true_negative <- (n_background -
                      true_positive - false_positive - false_negative)
  if (as_matrix){
    conting = c(true_positive, false_positive, false_negative, true_negative) %>%
      matrix(nrow = 2, ncol = 2, byrow = FALSE)
  }
  else{
    conting = data.frame(TP = true_positive, FP = false_positive,
                         FN = false_negative , TN = true_negative)
  }
  return(conting)
}
decouple_ORA_wrapper = function(marker_list,term_list, universe,
                                min_intersect = 3, seed = 42){
  # n_up = nrow(input_mat)
  # n_bottom = 0

  # withr::with_seed(seed, {
  #   targets <- decoupleR:::.ora_slice_targets(input_mat, n_up, n_bottom, with_ties)
  # })

  ORA_res = .ora_analysis(regulons = term_list, targets = marker_list,
                          universe = universe)

  ORA_conting = ORA_res %>% dplyr::select(source, condition, TP, FP, FN, TN)
  ORA_res = ORA_res %>% dplyr::select(statistic, source, condition, score, p_value)

  return(list(ORA_res, ORA_conting))
}

adjust_conting_iso = function(observed, expected,universe){
  # if (check_feat_type(observed) == "sf"){
  #   return(list("obs" = observed,
  #               "exp" = expected,
  #               "n_bg" = unique(universe)))
  # }
  FN = setdiff(expected, observed)
  TP = intersect(expected, observed)

  # obs_iso = metaspace_databases[which(metaspace_databases$name %in% observed),]
  obs_iso = metaspace_databases[metaspace_databases$name %fin% observed,]


  # TP_iso = obs_iso[which(obs_iso$name %in% TP),]
  TP_iso = obs_iso[obs_iso$name %fin% TP,]

  # univ_iso = metaspace_databases[which(metaspace_databases$name %in% universe),]
  univ_iso = metaspace_databases[metaspace_databases$name %fin% universe,]


  # FN_iso = metaspace_databases[which(metaspace_databases$name %in% FN),]
  # FN_iso = FN_iso[which(FN_iso$formula %nin% TP_iso$formula),]

  FN_iso = metaspace_databases[metaspace_databases$name %fin% FN,]
  FN_iso = FN_iso[FN_iso$formula %nin% TP_iso$formula,]

  # TN_iso = univ_iso[which(univ_iso$formula %nin% unique(c(FN_iso$formula,
  #                                                     obs_iso$formula))),]

  FN_iso = FN_iso[!duplicated(FN_iso$formula),]
  # TN_iso = TN_iso[!duplicated(TN_iso$formula),]

  TN_n = length(unique(univ_iso$formula)) - length(observed) - nrow(FN_iso)

  n_background = (length(observed) + nrow(FN_iso) + TN_n)

  return(list("obs" = observed,
              "exp" = c(unique(FN_iso$name), unique(TP_iso$name)),
              "n_bg" = n_background))
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
                    paste0("neg", annotation_adduct)) %in% colnames(exact_masses))
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
  names(bootstrapped_sublist) = c(1:length(bootstrapped_sublist))
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
    univ = universe
  }
  else{
    univ = unlist(term_list) %>% unique()
  }
  if (!is.null(selected_terms)){
    term_list = term_list[which(names(term_list) %in% selected_terms)]
  }



  enrich_res = decouple_ORA_wrapper(marker_list = bootstrap_list, term_list = term_list,
                                    universe = univ,
                                    min_intersect = min_annot)

  boot_conting_res = enrich_res[[2]]
  enrich_res = enrich_res[[1]]

  colnames(enrich_res)[which(colnames(enrich_res) == "source")] = "term"
  colnames(enrich_res)[which(colnames(enrich_res) == "condition")] = "bootstrap"

  colnames(boot_conting_res)[which(colnames(boot_conting_res) == "source")] = "term"
  colnames(boot_conting_res)[which(colnames(boot_conting_res) == "condition")] = "bootstrap"


  # boot_conting_res = pbapply::pblapply(seq(length(bootstrap_list)),function(n_i){
  #   enrich_res = hyper_geom_enrich(query = unique(bootstrap_list[[n_i]]),
  #                                  term_list = term_list, universe = universe)
  #   enrich_res$bootstrap = n_i
  #   #enrich_res = enrich_res %>% dplyr::filter(!is.na(OR), !is.na(pval))
  #   if(nrow(enrich_res) == 0){
  #     return(NULL)
  #   }
  #   else{
  #     return(enrich_res)
  #   }
  # })
  # boot_conting_res = boot_conting_res %>% dplyr::bind_rows()

  observed = boot_conting_res$TP / (boot_conting_res$TP + boot_conting_res$FP)
  expected = (boot_conting_res$TP + boot_conting_res$FN) / (boot_conting_res$TP + boot_conting_res$FP +
                                                              boot_conting_res$FN + boot_conting_res$TN)
  boot_conting_res$OR = observed / expected


  boot_enrich_res = boot_conting_res %>%
    dplyr::left_join(enrich_res, by = c("term","bootstrap")) %>%
    dplyr::mutate(padj = p.adjust(p_value, "BH")) %>%
    dplyr::group_by(term) %>%
    dplyr::mutate(fraction = length(term) / length(bootstrap_list)) %>%
    dplyr::ungroup()

  final_enrich_res = boot_enrich_res %>%
    dplyr::filter(TP >= min_annot, p_value < alpha_cutoff,
                  fraction > boot_fract_cutoff)

  final_enrich_res <- final_enrich_res %>%
    group_by(bootstrap) %>%
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
#' @export
Run_simple_ORA = function(marker_list, background, custom_universe = NULL,
                          alpha_cutoff = 0.05, min_intersection = 3){

  if (!is.null(custom_universe)){
    pathway_list_slim <- sapply(background, function(i){
      i[i %in% custom_universe]
    }, simplify = F)
    background <- pathway_list_slim[sapply(pathway_list_slim, length) > 0]
  }

  if (!is.list(marker_list)){
    q = sub("[-+].*","", marker_list) %>% unique()
    # ORA_conting = hyper_geom_enrich(query = q, term_list = background,
    #                            universe = custom_universe)
    ORA_res = decouple_ORA_wrapper(marker_list = list("Condition" = q),
                                   term_list = background,
                                   n_background = length(unique(unlist(background))),
                                   min_intersect = min_intersection)
  }
  else{
    ORA_res = decouple_ORA_wrapper(marker_list = marker_list,
                                   term_list = background,
                                   n_background = length(unique(unlist(background))),
                                   min_intersect = min_intersection)
  }
  ORA_conting = ORA_res[[2]]
  ORA_res = ORA_res[[1]]
  ORA_final = ORA_conting %>% dplyr::left_join(ORA_res)
  ORA_final = ORA_final %>%
    mutate(padj = p.adjust(p_value, "BH")) %>%
    dplyr::filter(TP >= min_intersection, p_value < alpha_cutoff)

  colnames(ORA_final)[which(colnames(ORA_final) == "source")] = "term"

  return(ORA_final)
}

#' @export
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

    if(!is.null(custom_universe)){
      univ_iso = get_metabo_iso(sf_vec = custom_universe, consider_isobars = consider_isobars, polarization_mode = polarization_mode,
                                mass_range_ppm = mass_range_ppm, only_HMDB = only_HMDB)
      custom_universe = univ_iso %>% unlist() %>% unique()
    }

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


