#' ncGTWinput-accessors
#'
#' Accessors to the feature and profiles loaded by \code{\link{loadProfile}}.
#' @param object a \code{\link{ncGTWinput}} object.
#' @return \code{groupInfo} returns a vector of the information of the loaded
#' feature.
#' @export
#' @rdname ncGTWinput-accessors
#' @aliases groupInfo
#' @examples
#' # obtain data
#' data('xcmsExamples')
#' xcmsLargeWin <- xcmsExamples$xcmsLargeWin
#' xcmsSmallWin <- xcmsExamples$xcmsSmallWin
#' ppm <- xcmsExamples$ppm
#'
#' # detect misaligned features
#' excluGroups <- misalignDetect(xcmsLargeWin, xcmsSmallWin, ppm)
#'
#' # obtain the paths of the sample files
#' filepath <- system.file("extdata", package = "ncGTW")
#' file <- list.files(filepath, pattern="mzxml", full.names=TRUE)
#'
#' tempInd <- matrix(0, length(file), 1)
#' for (n in seq_along(file)){
#'     tempCha <- file[n]
#'     tempLen <- nchar(tempCha)
#'     tempInd[n] <- as.numeric(substr(tempCha, regexpr("example", tempCha) + 7,
#'         tempLen - 6))
#' }
#' # sort the paths by data acquisition order
#' file <- file[sort.int(tempInd, index.return = TRUE)$ix]
#' \dontrun{
#' # load the sample profiles
#' ncGTWinputs <- loadProfile(file, excluGroups)
#'
#' gInfo <- groupInfo(ncGTWinputs[[1]])
#' prof <- profiles(ncGTWinputs[[1]])
#' rtR <- rtRaw(ncGTWinputs[[1]])
#' }

setMethod("groupInfo", "ncGTWinput", function(object) object@groupInfo)


#' @return \code{profiles} returns a raw data matrix in which each row is a
#' sample profile.
#' @export
#' @rdname ncGTWinput-accessors
#' @aliases profiles

setMethod("profiles", "ncGTWinput", function(object) object@profiles)


#' @return \code{rtRaw} returns a raw RT matrix in which each row is the
#' corresponding sample RT.
#' @export
#' @rdname ncGTWinput-accessors
#' @aliases rtRaw

setMethod("rtRaw", "ncGTWinput", function(object) object@rtRaw)



#' ncGTWoutput-accessors
#'
#' Accessors to the alignment information and result by ncGTW.
#' @param object a \code{\link{ncGTWoutput}} object.
#' @return \code{alignData} returns a matrix in which each row is a sample
#' profile after downsampling.
#' @export
#' @rdname ncGTWoutput-accessors
#' @aliases alignData
#' @examples
#' # obtain data
#' data('xcmsExamples')
#' xcmsLargeWin <- xcmsExamples$xcmsLargeWin
#' xcmsSmallWin <- xcmsExamples$xcmsSmallWin
#' ppm <- xcmsExamples$ppm
#'
#' # detect misaligned features
#' excluGroups <- misalignDetect(xcmsLargeWin, xcmsSmallWin, ppm)
#'
#' # obtain the paths of the sample files
#' filepath <- system.file("extdata", package = "ncGTW")
#' file <- list.files(filepath, pattern="mzxml", full.names=TRUE)
#'
#' tempInd <- matrix(0, length(file), 1)
#' for (n in seq_along(file)){
#'     tempCha <- file[n]
#'     tempLen <- nchar(tempCha)
#'     tempInd[n] <- as.numeric(substr(tempCha, regexpr("example", tempCha) + 7,
#'         tempLen - 6))
#' }
#' # sort the paths by data acquisition order
#' file <- file[sort.int(tempInd, index.return = TRUE)$ix]
#' \dontrun{
#' # load the sample profiles
#' ncGTWinputs <- loadProfile(file, excluGroups)
#'
#' # initialize the parameters of ncGTW alignment with default
#' ncGTWparam <- initncGTWparam()
#'
#' # run ncGTW alignment
#' ncGTWoutputs <- vector('list', length(ncGTWinputs))
#' for (n in seq_along(ncGTWinputs))
#'     ncGTWoutputs[[n]] <- ncGTWalign(ncGTWinputs[[n]], xcmsLargeWin, 5,
#'         ncGTWparam = ncGTWparam)
#'
#' data <- alignData(ncGTWoutputs[[1]])
#' rt <- scanRange(ncGTWoutputs[[1]])
#' paths <- ncGTWpath(ncGTWoutputs[[1]])
#' downSam <- downSample(ncGTWoutputs[[1]])
#' }

setMethod("alignData", "ncGTWoutput", function(object) object@alignData)


#' @return \code{scanRange} returns a downsampled RT matrix in which each row is
#' the corresponding sample RT in \code{data}.
#' @export
#' @rdname ncGTWinput-accessors
#' @aliases scanRange

setMethod("scanRange", "ncGTWoutput", function(object) object@scanRange)


#' @return \code{ncGTWpath} returns a list of the same length as the sample
#' number, in which each element is a matrix of the alignment result of the
#' corresponding sample.
#' @export
#' @rdname ncGTWinput-accessors
#' @aliases ncGTWpath

setMethod("ncGTWpath", "ncGTWoutput", function(object) object@ncGTWpath)


#' @return \code{downSample} returns the factor of downsampling when perform
#' ncGTW alignment.
#' @export
#' @rdname ncGTWinput-accessors
#' @aliases downSample

setMethod("downSample", "ncGTWoutput", function(object) object@downSample)



#' ncGTWwarp-accessors
#'
#' Accessors to the realigned RT and the peaks with adjusted RT of ncGTW
#' alignment.
#' @param object a \code{\link{ncGTWwarp}} object.
#' @return \code{rtncGTW} returns a list of the same length as the sample
#' number, in which each element is a vector of the realigned RT of the
#' corresponding sample.
#' @export
#' @rdname ncGTWwarp-accessors
#' @aliases rtncGTW
#' @examples
#' # obtain data
#' data('xcmsExamples')
#' xcmsLargeWin <- xcmsExamples$xcmsLargeWin
#' xcmsSmallWin <- xcmsExamples$xcmsSmallWin
#' ppm <- xcmsExamples$ppm
#'
#' # detect misaligned features
#' excluGroups <- misalignDetect(xcmsLargeWin, xcmsSmallWin, ppm)
#'
#' # obtain the paths of the sample files
#' filepath <- system.file("extdata", package = "ncGTW")
#' file <- list.files(filepath, pattern="mzxml", full.names=TRUE)
#'
#' tempInd <- matrix(0, length(file), 1)
#' for (n in seq_along(file)){
#'     tempCha <- file[n]
#'     tempLen <- nchar(tempCha)
#'     tempInd[n] <- as.numeric(substr(tempCha, regexpr("example", tempCha) + 7,
#'         tempLen - 6))
#' }
#' # sort the paths by data acquisition order
#' file <- file[sort.int(tempInd, index.return = TRUE)$ix]
#' \dontrun{
#' # load the sample profiles
#' ncGTWinputs <- loadProfile(file, excluGroups)
#'
#' # initialize the parameters of ncGTW alignment with default
#' ncGTWparam <- initncGTWparam()
#'
#' # run ncGTW alignment
#' ncGTWoutputs <- vector('list', length(ncGTWinputs))
#' ncGTWoutputs[[1]] <- ncGTWalign(ncGTWinputs[[1]], xcmsLargeWin, 5,
#'         ncGTWparam = ncGTWparam)
#'
#' # adjust RT with the realignment results from ncGTW
#' ncGTWres <- xcmsLargeWin
#' adjustRes <- adjustRT(ncGTWres, ncGTWinputs[[1]], ncGTWoutputs[[1]], ppm)
#' rt <- rtncGTW(adjustRes)
#' peaks <- ncGTWpeaks(adjustRes)
#' }

setMethod("rtncGTW", "ncGTWwarp", function(object) object@rtncGTW)


#' @return \code{rtncGTW} returns a matrix containing peak data with adjusted
#' RT.
#' @export
#' @rdname ncGTWwarp-accessors
#' @aliases ncGTWpeaks

setMethod("ncGTWpeaks", "ncGTWwarp", function(object) object@ncGTWpeaks)
