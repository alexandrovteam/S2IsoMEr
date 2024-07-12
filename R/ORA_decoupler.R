.ora_analysis <- function(regulons, targets, universe,pass_adjust = F, ...) {

  # NSE vs. R CMD check workaround
  p.value <- NULL

  message("\nRunning ORA analysis ... \n")

  n = length(names(regulons)) * length(names(targets))

  univ_iso = metaspace_databases[metaspace_databases$name %fin% universe,]

  pb = progress::progress_bar$new(total = n, incomplete = " ")

  tidyr::expand_grid(source = names(regulons), condition = names(targets)) %>%
    dplyr::rowwise(source, condition) %>%
    dplyr::summarise(.ora_fisher_exact_test(
      dat = adjust_conting_iso(
        expected = regulons[[source]],
        observed = targets[[condition]],
        universe_iso = univ_iso,
        pass = pass_adjust
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
                                pass_adjust = F, seed = 42){

  ORA_res = .ora_analysis(regulons = term_list, targets = marker_list,
                          universe = universe, pass_adjust = pass_adjust)

  ORA_conting = ORA_res %>% dplyr::select(source, condition, TP, FP, FN, TN)
  ORA_res = ORA_res %>% dplyr::select(statistic, source, condition, score, p_value)

  return(list(ORA_res, ORA_conting))
}
