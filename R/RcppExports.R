# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Functions that compute Enrichmnet Score
#'
#' @param stat NumericVector
#' @param geneOrd CharacterVector
#' @param Pmiss double
#' @param Nr double
#' @param poz NumericVector
#' @export
#' @rdname ES
ES <- function(stat, geneOrd, Pmiss, Nr, poz) {
    .Call('gseak_ES', PACKAGE = 'gseak', stat, geneOrd, Pmiss, Nr, poz)
}

#' ginGS
#'
#' @param g CharacterVector
#' @param gs CharacterVector
#' @param stat NumericVector
#' @export
#' @rdname ginGS
ginGS <- function(g, gs, stat) {
    .Call('gseak_ginGS', PACKAGE = 'gseak', g, gs, stat)
}

#' Function to compute NES
#'
#' @param x ES_p vector
#' @param mp mean value of ES_p>=0
#' @param mn mean value of ES_p<0
#' @export
#' @rdname NES
NES <- function(x, mp, mn) {
    .Call('gseak_NES', PACKAGE = 'gseak', x, mp, mn)
}

#' permutation
#'
#' @param x NumericVector
#' @param n int
#' @export
#' @rdname permutation
permutation <- function(x, n) {
    .Call('gseak_permutation', PACKAGE = 'gseak', x, n)
}

