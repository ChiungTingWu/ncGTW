#' Class "ncGTWparam"
#'
#' An S4 class for storing the needed paramters of ncGTW alignment.
#'
#' @slot downSample A factor of downsampling. The larger, the faster speed
#'   of alignment, but the accuracy may decrease. The default is 2.
#' @slot stpRat A factor to control the maximum RT shift of a point in
#'   alignment, and the maximum shift is determined by \code{stpRat} * "The RT
#'   range of the feature". If \code{maxStp} is set, then \code{stpRat} would be
#'   neglected. The default is 0.6.
#' @slot maxStp A value determines the maximum RT shift of a point in ncGTW
#'   alignment. If the user wants to decide the maximum shif by the RT range of
#'   the feature, this argument should be NaN. The default is NaN.
#' @slot strNum A value controls how many neighboring warping functions are
#'   connected to each warping function in ncGTW graph. There are two samples
#'   corresponding to a warping function, and at least one sample should be the
#'   same in another warping function to be considered as a neighbor controlled
#'   by \code{strNum}.
#' @slot diaNum A value controls how many neighboring warping functions are
#'   connected to each warping function in ncGTW graph. There are two samples
#'   corresponding to a warping function, and the two samples could also be
#'   different to another warping function to be considered as a neighbor
#'   controlled by \code{diaNum}.
#' @slot nor A value controls p-norm to compute the distance between the
#'   points on the profiles, and the default is 1 (Manhattan norm).
#' @details This function initializes the needed paramters of ncGTW alignment
#' with defaults, so this function could be called without any input. The
#' alignment should be fine with all default parameters. If the computing time
#' is an issue, the user could consider increase \code{downSample} and/or
#' decrease \code{stpRat} for a faster speed. If the alignment result is not
#' good enough, one can consider increase \code{strNum} and/or \code{diaNum}
#' to integrate more neighboring information to increase the quality of
#' alignment, but the speed may drop.

ncGTWparam <- setClass("ncGTWparam", slots=c(downSample="numeric",
                        stpRat="numeric", maxStp="numeric",
                        strNum="numeric", diaNum="numeric", nor="numeric")
)


#' Class "ncGTWinput"
#'
#' An S4 class for storing the inputs of ncGTW alignment.
#'
#' @slot groupInfo  A vector of the information of the feature.
#' @slot profiles   A raw data matrix in which each row is a sample profile.
#' @slot rtRaw      A raw RT matrix in which each row is the corresponding
#'   sample RT in \code{profiles}.

ncGTWinput <- setClass("ncGTWinput", slots=c(groupInfo="vector",
    profiles="matrix", rtRaw="matrix")
)

valid_ncGTWinput <- function(object) {
    errors <- character()

    cName <- c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax")
    if (!all(cName %in% names(object@groupInfo))){
        msg <- paste0("'groupInfo' must contain names: ",
                    paste(cName,collapse=', '), '.')
        errors <- c(errors, msg)
    }

    if (!all(dim(object@profiles) == dim(object@rtRaw))) {
        msg <- "The dimension of 'profiles' and 'rtRaw' must be the same."
        errors <- c(errors, msg)
    }

    if (length(errors) == 0) TRUE else errors
}

setValidity("ncGTWinput", valid_ncGTWinput)


#' Class "ncGTWoutput"
#'
#' An S4 class for storing the outputs of ncGTW alignment.
#'
#' @slot alignData  A matrix in which each row is a sample profile after
#'   downsampling.
#' @slot scanRange  A downsampled RT matrix in which each row is the
#'   corresponding sample RT in \code{alignData}.
#' @slot path       A list of the same length as the sample number, in which
#'   each element is a matrix of the alignment result of the corresponding
#'   sample.
#' @slot downSample      The factor of downsampling when perform ncGTW
#'   alignment.

ncGTWoutput <- setClass("ncGTWoutput", slots=c(alignData="matrix",
    scanRange="matrix", ncGTWpath="list", downSample="numeric")
)

valid_ncGTWoutput <- function(object) {
    errors <- character()

    if (!all(dim(object@alignData) == dim(object@scanRange))) {
        msg <- "The dimension of 'alignData' and 'scanRange' must be the same."
        errors <- c(errors, msg)
    }

    if (nrow(object@alignData) != length(object@ncGTWpath)) {
        msg <- paste("The row number of 'alignData' and the length of",
                "'ncGTWpath' must be the same (sample number).")
        errors <- c(errors, msg)
    }

    if (object@downSample < 1) {
        msg <- paste("The value of 'downSample' must be larger or equal to 1.")
        errors <- c(errors, msg)
    }

    if (length(errors) == 0) TRUE else errors
}

setValidity("ncGTWoutput", valid_ncGTWoutput)


#' Class "ncGTWwarp"
#'
#' An S4 class for storing the realigned RT and the peaks with adjusted RT of
#' ncGTW alignment.
#'
#' @slot rtncGTW    A list of the same length as the sample number, in which
#'   each element is a vector of the realigned RT of the corresponding sample.
#' @slot ncGTWpeaks A matrix containing peak data with adjusted RT.

ncGTWwarp <- setClass("ncGTWwarp", slots=c(rtncGTW="list", ncGTWpeaks="matrix"))

valid_ncGTWwarp <- function(object) {
    errors <- character()

    cName <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into", "intb",
                "maxo", "sn", "sample")
    if (!all(cName %in% colnames(object@ncGTWpeaks))){
        msg <- paste0("'ncGTWpeaks' must contain column names: ",
            paste(cName,collapse=', '), '.')
        errors <- c(errors, msg)
    }

    if (length(errors) == 0) TRUE else errors
}

setValidity("ncGTWwarp", valid_ncGTWwarp)
