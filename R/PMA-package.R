#' Penalized Multivariate Analysis
#'
#' This package is called __PMA__, for __P__enalized __M__ultivariate
#' __A__nalysis.  It implements three methods: A penalized matrix
#' decomposition, sparse principal components analysis, and sparse
#' canonical correlations analysis. All are described in the reference below.
#' The main functions are: `PMD`, `CCA` and `SPC`.
#'
#' The first, `PMD`, performs a penalized matrix decomposition.  `CCA`
#' performs sparse canonical correlation analysis. `SPC` performs sparse
#' principal components analysis.
#'
#' There also are cross-validation functions for tuning parameter selection for
#' each of the above methods: `SPC.cv`, `PMD.cv`, `CCA.permute`. And `PlotCGH` produces
#' nice plots for DNA copy number data.
#'
#' @name PMA-package
#' @aliases PMA-package PMA
#' @docType package
#' @useDynLib PMA
#' @references
#' Witten D. M., Tibshirani R.,  and Hastie, T. (2009) \doi{10.1093/biostatistics/kxp008}.
#' @keywords package
NULL

