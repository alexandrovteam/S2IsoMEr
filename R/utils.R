`%fin%` <- function(x, table) {
  fastmatch::fmatch(x, table, nomatch = 0L) > 0L
}

"%nin%" = Negate("%fin%")

to_ORA_query = function(dat, metadata = NULL, DE_data = T,
                        min.pct.exp = 0.5){
  if (DE_data){
    #TODO Check if we need to filer by LFC and pval first
    final_markers_per_cond = dat %>%
      dplyr::group_by(condition) %>%
      dplyr::summarise(
        markers = rlang::set_names(list(markers), condition[1]),
        .groups = "drop"
      ) %>%
      dplyr::pull(markers)
  }
  else{
    if (is.matrix(dat) & !is.null(metadata)){
      message("We recommend performing differential expression before ORA enrichment")
      final_markers_per_cond = list()
      for (cond in unique(metadata$condition)){
        cond_cells = metadata$cell[which(metadata$condition == cond)]
        subset = dat[,which(colnames(dat) %in% cond_cells)]
        pct_detected = data.frame(marker = rownames(subset),
                                  pct.exp = apply(subset, 1, function(x){
                                    length(which(x != 0)) / length(x)
                                  }),
                                  condition = cond)
        markers_per_cond = pct_detected$marker[which(pct_detected$pct.exp >= min.pct.exp)]
        if (length(markers_per_cond) == 0){
          final_markers_per_cond[[cond]] = character(0)
        }
        else{
          final_markers_per_cond[[cond]] = markers_per_cond
        }
      }
    }
    else if (is.data.frame(dat)){
      final_markers_per_cond = dat %>%
        dplyr::group_by(condition) %>%
        dplyr::summarise(
          markers = rlang::set_names(list(markers), condition[1]),
          .groups = "drop"
        ) %>%
        dplyr::pull(markers)
    }
    else if (is.matrix(dat) & is.null(metadata)){
      stop("Cannot map cells to condition, please provide cell metadata")
    }
  }

  return(final_markers_per_cond)
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
                "Metabo_super_class_name" = Metabo_super_class_name,
                "Metabo_main_class_name" = Metabo_class_name,
                "Metabo_sub_class_name" = Metabo_sub_class_name,
                "Metabo_pathways_name" = Metabo_pathway_name)
  }
  return(bg)
}


#' Set conditions for enrichment analysis
#'
#' @param object A bmetenrichr object.
#' @param condition.x A optional character describing the reference condition.
#' @param condition.y A optional character describing condition to interest.
#'
#' @return An object of class bmetenrich.
#'
#' @return An object of class bmetenrich.
#' @examples
#'
#' setConditions(myTestRun, condition.x = 'CON', condition.y = "TREATMENT")
#'
#' @export
setConditions <- function (object, ...) {
  UseMethod("setConditions", object)
}

#' @export
setConditions.bmetenrich <- function(object, condition.x = NULL, condition.y = NULL){
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

