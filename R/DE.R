PackageCheck <- function(..., error = TRUE) {
  pkgs <- unlist(x = c(...), use.names = FALSE)
  package.installed <- vapply(
    X = pkgs,
    FUN = requireNamespace,
    FUN.VALUE = logical(length = 1L),
    quietly = TRUE
  )
  if (error && any(!package.installed)) {
    stop(
      "Cannot find the following packages: ",
      paste(pkgs[!package.installed], collapse = ', '),
      ". Please install"
    )
  }
  invisible(x = package.installed)
}


#' @importFrom pbapply pbsapply
#' @importFrom future.apply future_sapply
#' @importFrom future nbrOfWorkers
#' @importFrom limma rankSumTestWithCorrelation
seurat_wilcoxDETest <- function(data.use,cells.1,cells.2,verbose = TRUE,...){
  data.use <- data.use[, c(cells.1, cells.2), drop = FALSE]
  j <- seq_len(length.out = length(x = cells.1))

  my.sapply <- ifelse(
    test = verbose && nbrOfWorkers() == 1,
    yes = pbsapply,
    no = future_sapply
  )
  overflow.check <- ifelse(
    test = is.na(x = suppressWarnings(length(x = data.use[1, ]) * length(x = data.use[1, ]))),
    yes = FALSE,
    no = TRUE
  )
  limma.check <- PackageCheck("limma", error = FALSE)
  if (limma.check[1] && overflow.check) {
    p_val <- my.sapply(
      X = 1:nrow(x = data.use),
      FUN = function(x) {
        return(min(2 * min(limma::rankSumTestWithCorrelation(index = j, statistics = data.use[x, ])), 1))
      }
    )
  } else {
    if (getOption('Seurat.limma.wilcox.msg', TRUE) && overflow.check) {
      message(
        "For a more efficient implementation of the Wilcoxon Rank Sum Test,",
        "\n(default method for FindMarkers) please install the limma package",
        "\n--------------------------------------------",
        "\ninstall.packages('BiocManager')",
        "\nBiocManager::install('limma')",
        "\n--------------------------------------------",
        "\nAfter installation of limma, Seurat will automatically use the more ",
        "\nefficient implementation (no further action necessary).",
        "\nThis message will be shown once per session"
      )
      options(Seurat.limma.wilcox.msg = FALSE)
    }
    group.info <- data.frame(row.names = c(cells.1, cells.2))
    group.info[cells.1, "group"] <- "Group1"
    group.info[cells.2, "group"] <- "Group2"
    group.info[, "group"] <- factor(x = group.info[, "group"])
    data.use <- data.use[, rownames(x = group.info), drop = FALSE]
    p_val <- my.sapply(
      X = 1:nrow(x = data.use),
      FUN = function(x) {
        return(wilcox.test(data.use[x, ] ~ group.info[, "group"], ...)$p.value)
      }
    )
  }
  return(data.frame(p_val, row.names = rownames(x = data.use)))
}
