###


#' Generate bmetenrichr enrichment object
#'
#' initEnrichment() creates object to perform bootstrapping metabolite set enrichment analysis
#'
#' @param scmatrix A numeric matrix of n metabolites (rows) and m cells or measurments (columns).
#' @param conditions A vector of length m with condition identifiers.
#' @param enrichment_type A character specifying whether Overrepresentation analysis (ORA) or Metabolite Set Enrichment Analysis (MSEA) is performed
#' @param annot_db A character or character vector specifying which annotation database(s) were used. Current databases include ("CoreMetabolome", "HMDB","SwissLipids","LipidMaps"). Multiple databases can be selected
#' @param annot_custom_db An optional dataframe from which the isomers will be sampled. It should have 2 columns : formula and molecule name. If provided, `annot_db` will be ignored.
#' @param endogenous_only A logical indicating whether to only consider endogenous metabolites (default = TRUE).
#' @param pathway_assoc_only A logical indicating whether to only consider metabolites associated with a biological pathway (default = FALSE)
#' @param remove_expected_predicted A logical indicating whether to remove expected and predicted isomers based on HMDB status (default = TRUE)
#' @param annotations An optional custom list of length n, with each element contains a vector of isomer names. If not specified,
#' bmetenrichr uses the CoreMetabolome, LIPIDMAPS, SwissLipids, and HMDB databases from METASPACE (https://metaspace2020.eu/) to generate an annotation list automatically.
#' @param annotation.weights An optional list of length n, each element contains a vector of isomer weights. Only when annotations is provided as list.
#' @param consider_isobars A logical indicating whether to include isobars (default = FALSE)
#' @param consider_isomers A logical indicating whether to include isomers (default = TRUE)
#' @param mass_range_ppm A numeric indicating the mass range in ppm (default: mass_range_ppm = 3). Molecular formulas + adducts within this range will be treated as isobars. Only required when isobars = TRUE.
#' @param polarization_mode A character with either 'positive' (default) or 'negative'. Only required when isobars = TRUE. When set to 'positive', included adducts are '+H', '+Na', and '+K'.
#' When set to 'negative', included adducts are '-H', '+Cl'.
#' @param include An optional logical vector of length n indicating whether to include the annotations in the analysis.
#' @param molecule_type A character specifying the feature type for the enrichment background, either Metabolites or Lipids. Valid choices are c("Metabo", "Lipid").
#' @param background_type A character specifying the background type for enrichment, choose one from c("LION","main_class","super_class","sub_class","pathways").
#' @param custom_bg A named list with character vectors of metabolite names. Default is NULL
#' @param condition.x first condition identifier for pairwise comparison.
#' @param condition.y second condition identifier for pairwise comparison.
#' @param termsOfInterest A character containing 'selection' (for default LION-term selection), 'all', or a vector of term names (see 'pathway').
#' @param ranking.by A character of either 't.test' or 'wilcox.test', to rank metabolites for the respective statistic.Ignored if [enrichment_type] is 'ORA'.
#' @param gsea.method A character of either 'ks_signed' or 'fgsea'. Ignored if [enrichment_type] is 'ORA'.
#'
#' @return An object of class bmetenrich.
#' @examples
#' myTestRun <-
#' initEnrichment(scmatrix = scMatrix,
#'                    annotations = my_annotations,
#'                    annotation.weights = my_weights,
#'                    conditions = my_conditions,
#'                    condition.x = "A",
#'                    condition.y = "B" )
#'
#'
#' @export
initEnrichment <- function(scmatrix,
                           conditions,
                           enrichment_type = c("ORA", "MSEA"),
                           annot_db = "HMDB",
                           annot_custom_db = NULL,
                           endogenous_only = T,
                           pathway_assoc_only = F,
                           remove_expected_predicted = T,
                           annotations = NULL,
                           annotation.weights = NULL,
                           consider_isomers = TRUE,
                           consider_isobars = FALSE,
                           mass_range_ppm = 3,
                           polarization_mode = "positive",
                           include = NULL,
                           molecule_type = c("Lipid", "Metabo"),
                           background_type = c("LION","main_class","super_class",
                                               "sub_class","pathways"),
                           custom_bg = NULL,
                           termsOfInterest = "all",
                           condition.x = NULL,
                           condition.y = NULL,
                           ranking.by = c("wilcox.test","t.test",
                                          "BWS", "logFC"),
                           gsea.method = c("fgsea","ks_signed")){


  if (!is.null(annotations) & dim(scmatrix)[1] !=  length(annotations)){
    stop("single-cell matrix and annotations do not have the same length")
  }
  if (!is.null(annotations) & !is.null(include) & length(annotations) != length(include)){
    stop("annotations and include do not have the same length")
  }
  if (dim(scmatrix)[2] !=  length(conditions)){
    stop("single-cell matrix and conditions do not have the same length")
  }
  if (consider_isobars & is.list(annotations)){
    stop("isobars = TRUE in combination with a custom list of annotations is not supported")
  }
  if(consider_isobars){
    if(!polarization_mode %in% c('positive','negative')){
      stop("polarization_mode should be either 'positive' or 'negative'")
    } else {
      message(paste0("polarization_mode is: ", polarization_mode))
    }
  }


  if (!is.null(annotation.weights)) {
    ## are the weights in the right format?

    if (!is.null(annotations) & length(annotation.weights) != length(annotations)) {
      stop("annotations and annotation.weights do not have the same length")
      }
    else {
      ## extra check
      if (!is.null(annotations) & !all(sapply(annotation.weights, length) == sapply(annotations, length))) {
        stop("annotations and annotation.weights do not have the same length ")
      }
    }
  }

  if (!is.null(condition.x) ){
    if(!any(condition.x == unique(conditions))){
      stop("submitted condition.x not found in dataset")
    }
  }

  if (!is.null(condition.y) ){
    if(!any(condition.y == unique(conditions))){
      stop("submitted condition.y not found in dataset")
    }
  }

  if (is.null(annot_custom_db) & annot_db %nin% c("CoreMetabolome",
                                                  "HMDB","SwissLipids","LipidMaps")){
    stop("Annotation database is invalid and no custom database was provided\n Please check the documentation for more info")
  }

  if ((molecule_type == "Lipid" & annot_db %in% c("CoreMetabolome", "HMDB")) ||
      (molecule_type == "Metabo" & annot_db %in% c("LipidMaps", "SwissLipids"))){
    warning("Background molecule type doesn't match choice of annotation database\n
            Results might be unreliable, interpret with caution")
  }
  if (molecule_type == "Lipid" & endogenous_only){
    warning("Endogenous only has less coverage with lipids \n
            We recommend switching it to FALSE for lipid annotations")
  }


  feat_type = check_feat_type(rownames(scmatrix))

  if (any(c(consider_isobars, consider_isomers))){
    feat_type = "name"
  }

  #NOTE pathway will be built based on molecule type (Lipids / Metabo) and background type, otherwise custom
  if (is.null(custom_bg)) {
    pathway_list = Load_background(mol_type = molecule_type,
                                   bg_type = background_type,
                                   feature_type = feat_type)
    pathway_list$all <- unique(unlist(pathway_list))
    if (background_type != "LION"){
      LUT <- data.frame(ID = names(pathway_list),
                        name = names(pathway_list))
    }
    else{
      LUT <- LION_LUT
      }
    }
  else {
    if (is.list(custom_bg)) {
      pathway_list <- custom_bg
      pathway_list$all <- unique(unlist(pathway_list))

      ## make an *ad-hoc* LUT
      LUT <- data.frame(ID = names(pathway_list),
                        name = names(pathway_list))
    }
    else {
      stop("Custom background is not in the right format")
    }
  }

  #NOTE annotations can be constructed based on scmatrix, so the argument was changed to only accept optional custom list from user
  if (!is.null(annotations) & is.list(annotations)){
    ## use provided list when used as input
    annotation_list <- annotations
  }
  else{
    annotations = rownames(scmatrix)

    ## when vector of molecular formulas is provided, generate list with molecular names

    ## change + or - adduct into .
    annotation_formulas_adduct <- gsub("\\+|\\-",".",annotations)

    annotation_formulas <- gsub("\\..+$","",annotation_formulas_adduct)   ## remove adduct

    annotation_adduct <- gsub("^.+\\.","",annotation_formulas_adduct)

    cat("\nParsing isomers...\n")


    if (!is.null(annot_custom_db)){
      iso_bg = annot_custom_db
      iso_bg$db = "CustomDB"
    }
    else{
      metasp_r_idx = c()
      for (i in annot_db){
        if (i == "LipidMaps"){
          if (background_type == "LION"){
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


    #NOTE Adding Isomers for Metabo and Lipids
    if (consider_isomers){
      annotation_list <-
        lapply(annotation_formulas, function(annotation_formula_i){
          iso_bg$name[which(iso_bg$formula == annotation_formula_i)]
        })
    }
    else{
      annotation_list = vector("list", length(annotation_formulas))
      names(annotation_list) = annotation_formulas
      annotation_list = lapply(annotation_list, function(x){character(0)})
    }
  }

  if(consider_isobars){
    cat("Parsing potential isobars...\n\n")

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


    annotation_list <-
      lapply(seq_along(annotation_formulas), function(i){
        c(isomer = annotation_list[[i]],
          isobar = iso_bg$name[which(iso_bg$formula %in% gsub("\\..+$","",isobars_list[[i]]))])
      })
  }
  else {   ## isobars == FALSE
    isobars_list <- NULL
  }


  if(termsOfInterest == "selection" & background_type == "LION"){

    ## filter pathway_list by terms of interest
    # termsSelection <- read.csv(file = 'data-raw/LION_selection.csv')
    pathway_list <- pathway_list[names(pathway_list) %in% c("all", termsSelection$LION_ID)]

  }
  else if(termsOfInterest == "all"){
    termsSelection <- names(pathway_list)
  }
  else if(termsOfInterest == "selection" & background_type != "LION"){
    stop("Setting termsOfInterest = 'selection' is only valid for Background_type = 'LION'\nUse termsOfInterest = 'all' instead")
  }
  else {
    pathway_list <- pathway_list[termsOfInterest]
    termsSelection <- names(pathway_list)
  }


  #NOTE added enrichment_type, background_name and Nullified ranking and GSEA method if ORA enrichment is chosen
  object <-
    structure(
      list(scmatrix = scmatrix,
           enrichment_type = match.arg(enrichment_type),
           Annotation_database = annot_db,
           Custom_database = annot_custom_db,
           annotations = annotation_list,
           annotation.weights = annotation.weights,
           isobars_list = isobars_list,
           conditions = conditions,
           include = include,
           consider_isomers = consider_isomers,
           consider_isobars = consider_isobars,
           background_name = ifelse(!is.null(custom_bg), "custom",
                                    paste0(molecule_type, "_", background_type, "_",
                                           feat_type)),
           pathway_list = pathway_list,
           LUT = LUT,
           termsSelection = termsSelection,
           condition.x = condition.x,
           condition.y = condition.y,
           ranking.by = ifelse(enrichment_type == "MSEA", match.arg(ranking.by), NULL),
           gsea.method = ifelse(enrichment_type == "MSEA", match.arg(gsea.method), NULL)),
      class = "bmetenrich")

  print(object)
  return(object)
}


#' @export
print.bmetenrich <- function(object){
  cat("single-cell metabolomics matrix of", dim(object$scmatrix)[1], "metabolites and",
      dim(object$scmatrix)[2], "cells\n")
  cat("active pathway:", object$background_type ,"\n\n")

  cat("conditions:", paste(unique(object$conditions), collapse = ", "),"\n\n")

  cat("condition.x:", object$condition.x,"\n")
  cat("condition.y:", object$condition.y,"\n")
}


#' Run Metabolite/Lipids enrichment analysis for single cell metabolomics
#'
#' Run_enrichment() A wrapper to run either ORA or MSEA based on the initialized enrichment object
#'
#' @param object A numeric matrix of n metabolites (rows) and m cells or measurments (columns).
#' @param Run_DE A logical indicating whether to run differential analysis using limma's rank Sum Test With Correlation. Ignored if enrichment type is 'MSEA'.
#' @param DE_LFC_cutoff A numeric indicating the minimum log2 fold change for differential analysis
#' @param min.pct.diff A numeric indicating the minimum percentage difference between samples/cells in both conditions for a marker to be considered differentially abundant.
#' @return Data.frame with enrichment results
#' @examples
#' myTestRun <-
#' Run_enrichment(object = object,
#'                    Run_DE = TRUE)
#'
#'
#' @export
Run_enrichment <- function(object, Run_DE = FALSE,
                           DE_pval_cutoff = 0.05, DE_LFC_cutoff = 1,
                           min.pct.diff = 0.1,...){
  if(object$enrichment_type == "ORA"){
    cond_x_cells = colnames(object$scmatrix)[which(object$conditions == condition.x)]
    cond_y_cells = colnames(object$scmatrix)[which(object$conditions == condition.y)]
    LFC = apply(object$scmatrix, 1, function(x){
      log2(mean(x[cond_y_cells])) - log2(mean(x[cond_x_cells]))
    })
    pct.exp_cond_x = apply(object$scmatrix, 1, function(x){
      length(which(x[cond_x_cells] != 0)) / length(cond_x_cells)
    })
    pct.exp_cond_y = apply(object$scmatrix, 1, function(x){
      length(which(x[cond_y_cells] != 0)) / length(cond_y_cells)
    })
    if (Run_DE){

      pvals = seurat_wilcoxDETest(data.use = object$scmatrix,
                                  cells.1 = cond_x_cells,
                                  cells.2 = cond_y_cells)
      DE_res = pvals
      colnames(DE_res)[1] = "p.val"
      DE_res$LFC = LFC
      DE_res$pct.exp.1 = pct.exp_cond_x
      DE_res$pct.exp.2 = pct.exp_cond_y
      DE_res$pct.diff = DE_res$pct.exp.2 - DE_res$pct.exp.1
      sel_markers = DE_res %>% dplyr::filter(p.val < DE_pval_cutoff,
                                             abs(LFC) >= DE_LFC_cutoff,
                                             abs(pct.diff) >= min.pct.diff)
    }
    else{
      DE_res = cbind(LFC, pct.exp_cond_x, pct.exp_cond_y) %>% as.data.frame()
      colnames(DE_res) = c("LFC", "pct.exp.1", "pct.exp.2")
      DE_res$pct.diff = DE_res$pct.exp.2 - DE_res$pct.exp.1
      sel_markers = DE_res %>% dplyr::filter(abs(LFC) >= DE_LFC_cutoff,
                                             abs(pct.diff) >= min.pct.diff)
    }
    if (nrow(sel_markers) == 0){
      message("No differentially markers are found, maybe consider less strict cutoffs")
      return(NULL)
    }
    final_markers = list("upregulated" = rownames(sel_markers)[which(sel_markers$LFC > 0)],
                         "downregulated" = rownames(sel_markers)[which(sel_markers$LFC < 0)])
    if (any(c(object$consider_isomers, object$consider_isobars))){
      enrich_res = Run_bootstrap_ORA(marker_list = final_markers,
                                     background = object$pathway_list, ...)
    }
    else{
      enrich_res = Run_simple_ORA(marker_list = final_markers,
                                     background = object$pathway_list, ...)
    }
  }
  else{
    if (any(c(object$consider_isomers, object$consider_isobars))){
      enrich_res = Run_bootstrap_MSEA(object = object, ...)
    }
    else{
      enrich_res = Run_simple_MSEA(object = object, ranking_by = object$ranking.by,...)
    }
  }
  return(enrich_res)
}

