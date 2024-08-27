`%fin%` <- function(x, table) {
  fastmatch::fmatch(x, table, nomatch = 0L) > 0L
}

"%nin%" = Negate("%fin%")

mismatch_df = function (x, y, on = NULL){
  if (is.null(on)) {
    on <- dplyr::intersect(names(x), names(y))
    message("Matching on: ", paste(on, collapse = ", "))
  }
  keys <- plyr::join.keys(x, y, on)
  x[keys$x %nin% keys$y, , drop = FALSE]
}

#' @importFrom rlang :=
rename_if_exists <- function(data, old_name, new_name) {
  if (old_name %in% colnames(data)) {
    data <- dplyr::rename(data, !!new_name := !!old_name)
  }
  data
}

list_backgrounds = function(feat_type){
  all_bgs = c("Lipid_LION","Lipid_super_class", "Lipid_main_class", "Lipid_sub_class","Lipid_pathways",
    "Metabo_super_class", "Metabo_main_class", "Metabo_sub_class", "Metabo_pathways")
  all_bgs = paste0(all_bgs, "_" ,feat_type)
  all_bgs = all_bgs[all_bgs %nin% c("Lipid_LION_sf")]
  return(all_bgs)
}

check_feat_type = function(feats){
  if (length(intersect(feats,metaspace_databases$name)) > 1 ){
    feat_type = "name"
  }
  else{
    feat_type = "sf"
  }
  return(feat_type)
}

#' Load enrichment background
#'
#' Load_background() retrieves a list of terms and associated metabolites for enrichment based on requested
#' molecule type, background type and input feature type.
#'
#' @param mol_type Character indicating molecule type. Can be either `Lipid` for Lipids and `Metabo` for small molecules
#' @param bg_type Character indicating background type :
#' \describe{
#'  \item{`"LION"`}{LION ontology for lipids only.}
#'  \item{`"main_class"`}{Level 1 classification - finer classification compared to super_class.}
#'  \item{`"super_class"`}{Level 0 classification.}
#'  \item{`"sub_class"`}{Level 2 classification - finer classification compared to main_class.}
#'  \item{`"pathways"`}{Biological pathways based on KEGG, Reactome, and SMPDB.}
#' }
#' @param feature_type Character indicating input feature type. Can be either `sf` for sum formula or ion
#' and `name` for molecule name.
#'
#' @return List where each element name correspond to a given term and values are
#' either metabolite names or sum formulas.
#' @examples
#' \dontrun{
#' bg <-
#' Load_background(mol_type = "Metabo",
#'                    bg_type = "main_class",
#'                    feature_type = "name")
#' }
#'
#'
#' @export
Load_background = function(mol_type = c("Lipid", "Metabo"),
                           bg_type = c("LION","main_class","super_class", "sub_class",
                                       "pathways"),
                           feature_type = c("sf","name")){

  mol_type = match.arg(mol_type)
  bg_type = match.arg(bg_type)
  feature_type = match.arg(feature_type)

  if (mol_type == "Metabo" & bg_type == "LION"){
    stop("LION ontology is for Lipids only, please check molecule type")
  }
  bg_name = paste0(mol_type, "_", bg_type, "_", feature_type)
  if (bg_name %nin% list_backgrounds(feature_type)){
    stop("Background not found, please check list_backgrounds()")
  }
  # if(bg_name == "Lipid_LION_name"){
  #
  # }
  if (feature_type == "sf"){
   bg = switch(bg_name,
               "Lipid_super_class_sf" = LipidMaps_category_sf,
               "Lipid_main_class_sf" = LipidMaps_main_class_sf,
               "Lipid_sub_class_sf" = LipidMaps_sub_class_sf,
               "Lipid_pathways_sf" = LipidMaps_pathway_sf,
               "Metabo_super_class_sf" = Metabo_super_class_sf,
               "Metabo_main_class_sf" = Metabo_class_sf,
               "Metabo_sub_class_sf" = Metabo_sub_class_sf,
               "Metabo_pathways_sf" = Metabo_pathway_sf)
  }
  if (feature_type == "name"){
    bg = switch(bg_name,
                "Lipid_super_class_name" = LipidMaps_category_name,
                "Lipid_main_class_name" = LipidMaps_main_class_name,
                "Lipid_sub_class_name" = LipidMaps_sub_class_name,
                "Lipid_pathways_name" = LipidMaps_pathway_name,
                "Lipid_LION_name" = pathway_list_LION,
                "Metabo_super_class_name" = Metabo_super_class_name,
                "Metabo_main_class_name" = Metabo_class_name,
                "Metabo_sub_class_name" = Metabo_sub_class_name,
                "Metabo_pathways_name" = Metabo_pathway_name)
  }
  return(bg)
}


#' Set conditions for enrichment analysis
#'
#' @param object A S2IsoMEr object.
#' @param condition.x A optional character describing the reference condition.
#' @param condition.y A optional character describing condition to interest.
#'
#' @return An object of class `S2IsoMEr`.
#' @examples
#' \dontrun{
#' data("example_MSEA_obj")
#' example_MSEA_obj = setConditions(example_MSEA_obj, condition.x = 'U', condition.y = "FI")
#' }
#'
#' @export
setConditions <- function (object, condition.x = NULL, condition.y = NULL) {
  UseMethod("setConditions", object)
}

#' @export
setConditions.S2IsoMEr <- function(object, condition.x = NULL, condition.y = NULL){
  if (is.null(condition.x) & is.null(condition.y)){
    stop("no condition identifiers submitted, with no defaults")
  }

  if (!is.null(condition.x) ){
    if(!any(condition.x == unique(object$conditions))){
      stop("submitted condition not found in dataset")
    } else {
      object$condition.x <- condition.x
    }

  }

  if (!is.null(condition.y) ){
    if(!any(condition.y == unique(object$conditions))){
      stop("submitted condition not found in dataset")
    } else {
      object$condition.y <- condition.y
    }

  }
  return(object)
}

calc_ambiguity = function(input_iso_list, weights = NULL){
  if (!is.null(weights)) {
    if (length(weights) != length(input_iso_list)) {
      stop("Input isomer list and weights do not have the same length")
    }
    else {
      ## extra check
      if (!all(sapply(weights, length) == sapply(input_iso_list, length))) {
        stop("Input isomer list and weights do not have the same length for each element")
      }
    }
  }

  n_sf = length(input_iso_list)

  message("Calculating ambiguity ...")
  ambig_score = pbapply::pblapply(seq(n_sf),function(n_i){
    n_mol_names = length(input_iso_list[[n_i]])
    if(n_mol_names < 2){
      entropy = 0
    }
    else{
      if(!is.null(weights)){
        mol_weights = weights[[n_i]]
      }
      else {
        mol_weights = rep(x = 1, times = n_mol_names)
      }

      probs = mol_weights / sum(mol_weights)
      entropy = -sum(probs * log(probs, base = 2))
    }

    entropy
  })
  names(ambig_score) = names(input_iso_list)
  return(ambig_score)
}

#' Report filters applied to bootstrap-based enrichment results
#'
#' @description This function filters enrichment results based on various criteria, such as minimum intersection, significance thresholds, and bootstrapping fractions, reporting which terms passed / didn't pass which filter.
#'
#' @param unfiltered_df A data frame containing unfiltered enrichment results from \code{\link{Run_bootstrap_ORA}} or \code{\link{Run_bootstrap_MSEA}}.
#' @param enrich_type A character string specifying the type of enrichment analysis. Must be either "ORA" (Over-Representation Analysis) or "MSEA" (Metabolite Set Enrichment Analysis).
#' @param min_intersection An integer specifying the minimum number of true positives (TP) required for a term to pass the filter.
#' @param alpha_cutoff A numeric value specifying the p-value cutoff for significance. Default is 0.05.
#' @param q.val_cutoff A numeric value specifying the q-value cutoff for significance. Default is 0.2.
#' @param boot_fract_cutoff A numeric value specifying the minimum fraction of bootstraps in which a term must be present to be considered.
#'
#' @return A data frame with binary columns indicating which terms pass the specified criteria. The data frame columns include:
#'   \itemize{
#'     \item Term
#'     \item min_TP (minimum true positives)
#'     \item significant_adj_boot (significant adjusted bootstrap p-value)
#'     \item significant_adj_terms (significant adjusted term p-value)
#'     \item pass_boot_fraction (pass bootstrapping fraction)
#'     \item pass_all_filts (Whether term passes all filters)
#'   }
#'
#' @details The function adjusts p-values using the False Discovery Rate (FDR) method and calculates combined p-values using the `metap::sumlog` function. It then applies several filters to determine which terms pass all criteria.
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
#' enrich_ORA_summary <- passed_filters_per_term(
#'   unfiltered_df = ORA_res$upregulated$unfiltered_enrich_res,
#'   enrich_type = "ORA"
#' )
#' }
#'
#' @export
passed_filters_per_term = function(unfiltered_df,
                                 enrich_type = c("ORA", "MSEA"),
                                 min_intersection = 3,
                                 alpha_cutoff = 0.05,
                                 q.val_cutoff = 0.2,
                                 boot_fract_cutoff = 0.5){
  if(!enrich_type %in% c("ORA","MSEA")){
    stop("Invalid enrichment type. Please use either ORA or MSEA")
  }
  else{
    if (enrich_type == "MSEA"){
      colnames(unfiltered_df)[which(colnames(unfiltered_df) == "n")] = "TP"
    }
  }


  pass_filts = unfiltered_df %>%
    dplyr::group_by(.data$bootstrap) %>%
    dplyr::mutate(q.value = stats::p.adjust(.data$p_value, method = "fdr"))  %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$Term) %>%
    dplyr::mutate(n = stats::median(.data$TP, na.rm = T),
                  p.value_combined = metap::sumlog(.data$p_value)[["p"]],
                  q.value_combined = metap::sumlog(.data$q.value)[["p"]],
                  min_TP = ifelse(.data$n >= min_intersection, 1, 0),
                  significant_adj_boot = ifelse(.data$p.value_combined < alpha_cutoff, 1, 0),
                  significant_adj_terms = ifelse(.data$q.value_combined < alpha_cutoff, 1, 0),
                  pass_boot_fraction = ifelse(.data$fraction > boot_fract_cutoff,1,0)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data$Term != "") %>%
    dplyr::select(.data$Term, .data$min_TP, .data$significant_adj_boot,
                  .data$significant_adj_terms, .data$pass_boot_fraction) %>%
    dplyr::distinct()

  pass_all = rowSums(pass_filts[,-1]) / (ncol(pass_filts) - 1)
  pass_all[pass_all < 1] = 0
  pass_filts$pass_all_filts = pass_all

  return(pass_filts)
}


match_LUT_to_Term = function(df, LUT){
  if (all(df$Term %in% LUT$name)){
    return(df)
  }
  else{
    updated_df = df %>%
      dplyr::left_join(LUT, by = c("Term" = "ID")) %>%
      dplyr::select(-.data$Term) %>%
      dplyr::rename("Term" = "name") %>%
      dplyr::relocate(.data$Term)

    return(updated_df)
  }
}


collapse_ORA_boot_multi_cond = function(ORA_boot_res_list){
  unfiltered_res = lapply(ORA_boot_res_list, function(x){x[["unfiltered_enrich_res"]]}) %>%
    dplyr::bind_rows(.id = "Condition")
  clean_res = lapply(ORA_boot_res_list, function(x){x[["clean_enrich_res"]]}) %>%
    dplyr::bind_rows(.id = "Condition")

  return(list("unfiltered_enrich_res" = unfiltered_res,
              "clean_enrich_res" = clean_res))
}

dotplot_label_func = function(n){
  function(str) {
    str <- gsub("_", " ", str)
    yulab.utils::str_wrap(str, n)
  }
}

convert_to_log10 <- function(vec) {
  # Check if the values are in log10 scale
  is_log10 <- function(x) all(log10(10^x) == x)

  # Apply the check and convert if necessary
  if (!is_log10(vec)) {
    vec <- log10(vec)
  }

  # Replace -Inf values with 0
  vec[is.infinite(vec) & vec < 0] <- 0

  return(vec)
}


