#' ncGTW: A package for detecting and aligning the misaligned features in
#' LC-MS data
#'
#' The purpose of ncGTW is to help XCMS for LC-MS data alignment. Currently,
#' ncGTW can detect the misaligned feature groups by XCMS, and the user can
#' choose to realign these feature groups by ncGTW or not.
#'
#' @references Wu, Chiung-Ting, et al. Alignment of LC-MS Profiles by
#' Neighbor-wise Compound-specific Graphical Time Warping with Misalignment
#' Detection. bioRxiv, 2019, 715334.
#'
#' @docType package
#' @name ncGTW-package
#' @aliases ncGTW
#'
#' @import xcms
#' @import BiocParallel
#' @import methods
#' @importFrom grDevices dev.off png rgb
#' @importFrom graphics matplot title
#' @importFrom stats approx convolve cor kmeans median p.adjust pnorm sd t.test
#'   var weighted.mean
NULL
