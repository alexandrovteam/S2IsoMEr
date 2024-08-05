#' Build Isomer/Isobar Background
#'
#' @description This function constructs an isomeric/isobaric background database based on specified annotation databases and additional filtering criteria.
#'
#' @param annot_db A character vector specifying the annotation databases to use. Options include "HMDB", "LipidMaps", "SwissLipids", and "CoreMetabolome".
#' @param annot_custom_db A custom annotation database. If not NULL, this custom database will be used instead of the specified annotation databases.
#' @param use_LION A logical indicating whether to use the LION ontology or internal LipidMaps backgrounds if "LipidMaps" is included in `annot_db`.
#' @param endogenous_only A logical indicating whether to include only endogenous metabolites. Applies only to "HMDB" and "CoreMetabolome".
#' @param pathway_assoc_only A logical indicating whether to include only metabolites associated with pathways. Applies only to "HMDB" and "CoreMetabolome".
#' @param remove_expected_predicted A logical indicating whether to remove metabolites with a status of "expected" or "predicted" in HMDB and CoreMetabolome databases.
#' @return A data frame containing the isomers and isobers background based on the specified criteria.
#' @details This function builds an isomer/isobar background by either using the specified annotation databases or a provided custom annotation database. Various filtering options are available to refine the background.
#' @examples
#' \dontrun{
#' iso_bg <- build_iso_bg(annot_db = "HMDB", use_LION = FALSE, endogenous_only = TRUE,
#'                        pathway_assoc_only = FALSE, remove_expected_predicted = TRUE)
#' custom_db <- data.frame(ID = 1:3, Name = c("Met1", "Met2", "Met3"), Endogenous = "Yes", Pathway_assoc = "Yes", HMDB_status = "known")
#' iso_bg_custom <- build_iso_bg(annot_custom_db = custom_db)
#' }
#' @export
build_iso_bg = function(annot_db = "HMDB",annot_custom_db = NULL,
                        use_LION = F, endogenous_only = T,
                        pathway_assoc_only = F,remove_expected_predicted = T){
  if (!is.null(annot_custom_db)){
    iso_bg = .check_annot_custom_db(annot_custom_db)
    iso_bg$db = "CustomDB"
  }
  else{
    metasp_r_idx = c()
    for (i in annot_db){
      if (i == "LipidMaps"){
        if (use_LION){
          indices = stringr::str_which(metaspace_databases$db,
                                       "LION_LipidMaps")
        }
        else{
          indices = stringr::str_which(metaspace_databases$db,
                                       "Ramp_LipidMaps")
        }
      }
      else{
        indices = switch(i, "HMDB" = stringr::str_which(metaspace_databases$db,
                                                        "HMDB"),
                         "SwissLipids" = stringr::str_which(metaspace_databases$db,
                                                            "swisslipids"),
                         "CoreMetabolome" = stringr::str_which(metaspace_databases$db,
                                                               "CoreMetabolome"))
      }
      metasp_r_idx = c(metasp_r_idx, indices)
    }
    metasp_r_idx = unique(metasp_r_idx)
    iso_bg = metaspace_databases[metasp_r_idx,]
    rownames(iso_bg) = NULL

    if(endogenous_only){
      iso_bg = iso_bg[iso_bg$Endogenous == "Yes",]
    }
    if (pathway_assoc_only){
      iso_bg = iso_bg[iso_bg$Pathway_assoc == "Yes",]
    }
    if (remove_expected_predicted){
      iso_bg = iso_bg[iso_bg$HMDB_status %nin% c("expected", "predicted"),]
    }
  }

  return(iso_bg)
}

#' Get Metabolite Isomers and Isobars
#'
#' @description This function retrieves metabolite isomers and isobars for a given list of sum formula, considering various parameters such as adducts, polarization mode, and mass range.
#'
#' @param sf_vec A character vector of metabolite formulas with optional adducts. The formulas should be formatted as `formula[+adduct]` or `formula[-adduct]`.
#' @param consider_isobars A logical indicating whether to consider isobars. If TRUE, the function will identify and include isobars based on the provided mass range.
#' @param polarization_mode A character string specifying the ionization mode. Options are "positive" or "negative". Only relevant if `consider_isobars` is TRUE.
#' @param mass_range_ppm A numeric value specifying the mass range in parts per million (ppm) for identifying isobars.
#' @param annot_db A character vector specifying the annotation databases to use. Options include "HMDB", "LipidMaps", "SwissLipids", and "CoreMetabolome".
#' @param annot_custom_db A custom annotation database. If not NULL, this custom database will be used instead of the specified annotation databases.
#' @param use_LION A logical indicating whether to use the LION version of the LipidMaps database if "LipidMaps" is included in `annot_db`.
#' @param endogenous_only A logical indicating whether to include only endogenous metabolites. Applies only to "HMDB" and "CoreMetabolome".
#' @param pathway_assoc_only A logical indicating whether to include only metabolites associated with pathways. Applies only to "HMDB" and "CoreMetabolome".
#' @param remove_expected_predicted A logical indicating whether to remove metabolites with a status of "expected" or "predicted" in HMDB and CoreMetabolome databases.
#' @return A named list where each element corresponds to a formula from `sf_vec`. Each element is a character vector containing metabolite names matching the formula and its isobars, if `consider_isobars` is TRUE.
#' @details This function processes a vector of metabolite formulas and retrieves the corresponding metabolite names from the specified annotation databases. If `consider_isobars` is TRUE, it will also find isobars based on the specified mass range and polarization mode.
#' @examples
#' \dontrun{
#' sf_vec <- c("C6H12O6+H", "C6H12O6-H")
#' results <- get_metabo_iso(sf_vec, consider_isobars = TRUE, polarization_mode = "positive", mass_range_ppm = 5)
#' }
#' @export
get_metabo_iso = function(sf_vec, consider_isobars = T,
                          polarization_mode = NA, mass_range_ppm = 3,
                          annot_db = "HMDB", annot_custom_db = NULL,
                          use_LION = F, endogenous_only = T,
                          pathway_assoc_only = F,
                          remove_expected_predicted = T){

  annotation_formulas_adduct <- gsub("\\+|\\-",".",sf_vec)
  annotation_formulas <- gsub("\\..+$","",annotation_formulas_adduct)
  annotation_adduct <- gsub("^.+\\.","",annotation_formulas_adduct)

  iso_bg = build_iso_bg(annot_db = annot_db ,annot_custom_db = annot_custom_db,
                        use_LION = use_LION, endogenous_only = endogenous_only,
                        pathway_assoc_only = pathway_assoc_only,remove_expected_predicted = remove_expected_predicted)

  annotation_list <-
    lapply(annotation_formulas, function(annotation_formula_i){
      iso_bg$name[iso_bg$formula == annotation_formula_i]
    })
  is_adduct = any(c(paste0("pos", annotation_adduct),
                    paste0("neg", annotation_adduct)) %in% colnames(exact_masses))
  if (consider_isobars &
      is_adduct &
      polarization_mode %in% c("negative", "positive")){
    switch(polarization_mode,
           'positive' = {
             col_name <- paste0("pos",annotation_adduct)

             exact_masses_slim <- exact_masses[,c(1,which(grepl("pos",colnames(exact_masses))))]
             colnames(exact_masses_slim) <- gsub("^pos","",colnames(exact_masses_slim))
             exact_masses_slim <- exact_masses_slim %>% tidyr::pivot_longer(cols = -1, values_to = "mass", names_to = "adduct")
             exact_masses_slim$formula_adduct <- paste0(exact_masses_slim$formula,".",exact_masses_slim$adduct)
           },
           'negative' = {
             col_name <- paste0("neg",annotation_adduct)

             exact_masses_slim <- exact_masses[,c(1,which(grepl("neg",colnames(exact_masses))))]
             colnames(exact_masses_slim) <- gsub("^neg","",colnames(exact_masses_slim))
             exact_masses_slim <- exact_masses_slim %>% tidyr::pivot_longer(cols = -1, values_to = "mass", names_to = "adduct")
             exact_masses_slim$formula_adduct <- paste0(exact_masses_slim$formula,".",exact_masses_slim$adduct)
           })
    exact_masses_annotations <-
      mapply(annotation_formulas_i = annotation_formulas,
             col_name_i = col_name,
             function(annotation_formulas_i, col_name_i) {
               mass <-
                 exact_masses[exact_masses$formula == annotation_formulas_i, col_name_i]
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
          isobar = iso_bg$name[iso_bg$formula %in% gsub("\\..+$","",isobars_list[[i]])])

      })
  }
  names(annotation_list) = sf_vec

  return(annotation_list)
}

#' Perform Bootstrap Sampling on Metabolite Annotations
#'
#' @description This function performs bootstrap sampling on metabolite annotations, allowing for resampling of molecules with or without associated weights. It generates multiple bootstrap samples for further analysis.
#'
#' @param annot_list A list of character vectors where each vector contains metabolites or molecules for a particular annotation. Each element in the list corresponds to a specific sum formula.
#' @param annot_weights A list of numeric vectors where each vector contains weights associated with the metabolites in `annot_list`. Must have the same length as `annot_list`. If NULL, all metabolites are equally weighted.
#' @param n_bootstraps An integer specifying the number of bootstrap samples to generate.
#' @return A list of length `n_bootstraps`, where each element is a list of bootstrap samples for each sum formula. Each bootstrap sample is a character vector of sampled metabolites.
#' @details This function is used to create multiple bootstrap samples of metabolite annotations by randomly sampling metabolites from each sum formula. If weights are provided, the sampling is performed according to these weights; otherwise, the sampling is uniform.
#'
#' @examples
#' \dontrun{
#' annot_list <- list(
#'   Sf1 = c("Metabolite1", "Metabolite2", "Metabolite3"),
#'   Sf2 = c("MetaboliteA", "MetaboliteB")
#' )
#' annot_weights <- list(
#'   Sf1 = c(0.1, 0.2, 0.7),
#'   Sf2 = c(0.5, 0.5)
#' )
#' bootstrapped_samples <- metabo_bootstrap(annot_list, annot_weights, n_bootstraps = 100)
#' }
#' @export
metabo_bootstrap = function(annot_list, annot_weights = NULL,
                            n_bootstraps = 50){
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
  # if (sample_core_metab_only){
  #   annot_list <- sapply(annot_list, function(i){
  #     i[i %in% core_metab$name]
  #   }, simplify = F)
  # }

  bootstrapped_sublist <- pbapply::pblapply(seq(n_bootstraps), function(n_i) {
    sapply(seq(length(annot_list)), function(row_number_i) {
      molecules_to_sample <- annot_list[[row_number_i]]
      mols <- NULL

      if(length(molecules_to_sample) == 0) {
        mols <- names(annot_list)[row_number_i]
      } else {
        if(!is.null(annot_weights)) {
          weights_to_sample <- annot_weights[[row_number_i]]
        } else {
          weights_to_sample <- rep(1, length(molecules_to_sample))
        }

        if(length(molecules_to_sample) != 1) {
          mols <- sample(x = molecules_to_sample, size = 1, prob = weights_to_sample)
        } else {
          mols <- molecules_to_sample
        }
      }

      names(mols) <- names(annot_list)[row_number_i]
      mols

    }, simplify = TRUE)
  })
  names(bootstrapped_sublist) = c(1:length(bootstrapped_sublist))
  return(bootstrapped_sublist)
}

adjust_conting_iso = function(observed, obs_iso,
                              expected, exp_iso,
                              universe_iso, pass = F){


  if (pass){
    return(list("obs" = observed,
                "exp" = expected,
                "n_bg" = length(unique(universe_iso$name))))
  }

  #https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
  #expected = m
  #observed = k
  #universe = m+n

  #TP = x
  #FP = k-x
  #FN = m-x
  #TN = n-(k-x)

  #Sanity checks:
  #TP + FP = Observed
  #TP + FN = Expected
  #FP + TN = Universe - expected
  #FN + TN = Universe - observed

  FN = setdiff(expected, observed)
  TP = intersect(expected, observed)

  TP_iso = obs_iso[obs_iso$name %fin% TP,]
  FN_iso = exp_iso[exp_iso$name %fin% FN,]

  #Since in bootstrapping only 1 isomer is selected for observed while in background,
  #there are multiple isomers of the same formula for a given term, It's removed
  #from false negatives if any of the isomers are observed as true positve.

  FN_iso = FN_iso[FN_iso$formula %nin% TP_iso$formula,]
  FN_iso = FN_iso[!duplicated(FN_iso$formula),]

  FP_iso = length(unique(obs_iso$formula)) - length(unique(TP_iso$formula))

  TN_n = (length(unique(universe_iso$formula)) - length(unique(obs_iso$formula)) ) - FP_iso

  n_background = (length(observed) + nrow(FN_iso) + TN_n)

  return(list("obs" = observed,
              "exp" = c(unique(FN_iso$name), unique(TP_iso$name)),
              "n_bg" = n_background))
}

.check_annot_custom_db = function(df){
  colnames_common = intersect(colnames(df), c("formula", "name"))
  if(length(colnames_common) != 2){
    stop("Invalid format for custom annotation database. Please make sure that it is a dataframe with 2 columns : formula and molecule name")
  }
  return(df)
}
