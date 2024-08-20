.ora_analysis <- function(regulons, targets, universe,pass_adjust = F, ...) {

  # NSE vs. R CMD check workaround
  p.value <- NULL

  message("\nRunning ORA analysis ... \n")

  n = length(names(regulons)) * length(names(targets))

  all_obs = unique(unlist(targets))
  all_exp = unique(unlist(regulons))

  all_obs_iso = metaspace_databases[metaspace_databases$name %fin% all_obs,]
  all_exp_iso = metaspace_databases[metaspace_databases$name %fin% all_exp,]

  univ_iso = metaspace_databases[metaspace_databases$name %fin% universe,]
  obs_iso_list = lapply(targets, function(i){all_obs_iso[all_obs_iso$name %fin% i,]})
  exp_iso_list = lapply(regulons, function(i){all_exp_iso[all_exp_iso$name %fin% i,]})

  pb = progress::progress_bar$new(total = n, incomplete = " ")

  tidyr::expand_grid(source = names(regulons), condition = names(targets)) %>%
    dplyr::rowwise(.data$source, .data$condition) %>%
    dplyr::summarise(.ora_fisher_exact_test(
      dat = adjust_conting_iso(
        expected = regulons[[.data$source]],
        exp_iso = exp_iso_list[[.data$source]],
        observed = targets[[.data$condition]],
        obs_iso = obs_iso_list[[.data$condition]],
        universe_iso = univ_iso,
        pass = pass_adjust
      ),
      pbar = pb,
      ...
    ),
    .groups = "drop"
    ) %>%
    dplyr::select(.data$source, .data$condition,
                  p_value = p.value, dplyr::everything()
    ) %>%
    dplyr::mutate(score = -log10(.data$p_value)) %>%
    tibble::add_column(statistic = "ora", .before = 1) %>%
    dplyr::select(.data$statistic, .data$source, .data$condition, .data$score,
                  .data$p_value,
                  .data$TP, .data$FP, .data$FN, .data$TN, .data$TP_markers)
}
.ora_analysis_simple <- function(regulons, targets, universe,pass_adjust = F, ...) {

  # NSE vs. R CMD check workaround
  p.value <- NULL

  message("\nRunning ORA analysis ... \n")

  n = length(names(regulons)) * length(names(targets))

  pb = progress::progress_bar$new(total = n, incomplete = " ")

  n_univ = length(unique(universe))

  tidyr::expand_grid(source = names(regulons), condition = names(targets)) %>%
    dplyr::rowwise(.data$source, .data$condition) %>%
    dplyr::summarise(.ora_fisher_exact_test(
      dat = list("obs" = targets[[.data$condition]],
                 "exp" = regulons[[.data$source]],
                 "n_bg" = n_univ),
      pbar = pb,
      ...
      ),
    .groups = "drop"
    ) %>%
    dplyr::select(.data$source, .data$condition,
                  p_value = p.value, dplyr::everything()
    ) %>%
    dplyr::mutate(score = -log10(.data$p_value)) %>%
    tibble::add_column(statistic = "ora", .before = 1) %>%
    dplyr::select(.data$statistic, .data$source, .data$condition, .data$score,
                  .data$p_value,
                  .data$TP, .data$FP, .data$FN, .data$TN,.data$TP_markers)
}

.ora_fisher_exact_test <- function(dat,pbar, ...) {
  pbar$tick()
  conting = ora_conting_decoupleR(dat, as_matrix = F) %>% tibble::as_tibble()
  rlang::exec(
    .fn = stats::fisher.test,
    x = matrix(as.integer(conting[,1:4]), nrow = 2, ncol = 2),
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

  TP_mols = observed[observed %fin% expected]

  true_positive <- length(TP_mols)
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
                         FN = false_negative , TN = true_negative,
                         TP_markers = ifelse(true_positive != 0,
                                             paste(names(TP_mols), collapse = ";"),
                                             NA))
  }
  return(conting)
}

#' Run ORA Analysis with Optional Bootstrapping
#'
#' @description This function performs over-representation analysis (ORA) on a given list of markers using either a bootstrapped or a simple method.
#' It is a wrapper around the `run_ora` function from the `decoupleR` package. This function is also used internally in the `Run_bootstrap_ORA` and `Run_simple_ORA` functions.
#'
#' @param marker_list A list of metabolites per condition. If there are no conditions, please provide it as list("condition" = c("metabo_1", "metabo_2", ...))
#' @param term_list A list of named terms where each term is a character vector of metabolites.
#' @param universe A character vector specifying the universe of metabolites.
#' @param pass_adjust A logical indicating whether to pass the adjustment parameter to the ORA function. If TRUE, the contingency table won't be adjusted. Only required for bootstrap-ORA.
#' @param seed An integer specifying the seed for reproducibility.
#' @param ORA_boot A logical indicating whether to use bootstrapped ORA. If TRUE, the function will use the bootstrapped ORA method; otherwise, it will use a simple ORA method.
#' @return A list containing two elements:
#' \item{ORA_res}{A data frame with the ORA results, including columns for statistic, source, condition, score, and p-value.}
#' \item{ORA_conting}{A data frame with the ORA contingency table, including columns for source, condition, TP, FP, FN, and TN.}
#' @details This function is a wrapper around the `run_ora` function from the `decoupleR` package. For more details, please refer to the [decoupleR documentation](https://saezlab.github.io/decoupleR/reference/run_ora.html).
#' @references
#' Badia-i-Mompel, P., Nagai, J. S., & Saez-Rodriguez, J. (2022). decoupleR: A flexible tool to handle various modes of biological network analysis. *Bioinformatics Advances*, 2(1), vbac016. [https://doi.org/10.1093/bioadv/vbac016](https://academic.oup.com/bioinformaticsadvances/article/2/1/vbac016/6544613)
#' @examples
#' \dontrun{
#' marker_list <- c("B2", "B5", "B8")
#' term_list <- list(
#'   "Term1" = c("B1", "B2", "B3", "B4"),
#'   "Term2" = c("B4", "B5", "B6", "B9")
#' )
#' universe <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9")
#' results <- decouple_ORA_wrapper(
#'   marker_list, term_list, universe,
#'   pass_adjust = TRUE, seed = 42, ORA_boot = FALSE
#' )
#' }
#' @export
decouple_ORA_wrapper = function(marker_list,term_list, universe,
                                pass_adjust = F, seed = 42, ORA_boot = T){

  if(ORA_boot){
    ORA_res = .ora_analysis(regulons = term_list, targets = marker_list,
                            universe = universe, pass_adjust = pass_adjust)
  }
  else{
    ORA_res = .ora_analysis_simple(regulons = term_list, targets = marker_list,
                            universe = universe)
  }


  ORA_conting = ORA_res %>% dplyr::select(.data$source, .data$condition, .data$TP,
                                          .data$FP, .data$FN, .data$TN, .data$TP_markers)
  ORA_res = ORA_res %>% dplyr::select(.data$statistic,
                                      .data$source, .data$condition, .data$score,
                                      .data$p_value)

  return(list(ORA_res, ORA_conting))
}
