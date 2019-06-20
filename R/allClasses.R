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


#' Class "ncGTWwarp"
#'
#' An S4 class for storing the realigned RT and the peaks with adjusted RT of
#' ncGTW alignment.
#'
#' @slot rtncGTW    A list of the same length as the sample number, in which
#'   each element is a vector of the realigned RT of the corresponding sample.
#' @slot ncGTWpeaks A matrix containing peak data with adjusted RT.

ncGTWwarp <- setClass("ncGTWwarp", slots=c(rtncGTW="list", ncGTWpeaks="matrix"))
