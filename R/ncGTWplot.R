#' Plot profiles for each peak group
#'
#' This function plots sample profiles with loaded information.
#' @param ncGTWinput An object return by \code{\link{loadProfile}} of sample
#'   profiles for plotting.
#' @param sampleRt A list of the same length as the sample number in which each
#'   element is a vector corresponding to the sample raw/adjusted RT for
#'   plotting.
#' @param sampleInd Indicate which samples should be plotted, and the default is
#'   \code{1:dim(ncGTWinput$rtRaw)[1]}.
#' @param ind A user defined index, and the default is \code{NULL}.
#' @param savePath The path to save the plots, and the default is \code{NULL}
#'   (do not save anything).
#' @param show Show the plot in R or not, and the default is \code{TRUE}.
#' @param sub Show more information on the plot or not, and the default is
#'   \code{TRUE}.
#' @param filter Apply a Gaussian filter for demonstration or not, and the
#'   default is \code{FALSE}.
#'
#' @details This function plots the extracted ion chromatogram obtained by
#'   \code{\link{loadProfile}}. The user can decide to save the figure, show the
#'   figure, or apply a Gaussian filter on the data by parameter setting.
#' @return A plot to the current device.
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
#' for (n in 1:length(file)){
#'   tempCha <- file[n]
#'   tempLen <- nchar(tempCha)
#'   tempInd[n] <- as.numeric(substr(tempCha, regexpr("example", tempCha) + 7, tempLen - 6))
#' }
#' # sort the paths by data acquisition order
#' file <- file[sort.int(tempInd, index.return = TRUE)$ix]
#'
#' # load the sample profiles
#' ncGTWinputs <- loadProfile(file, excluGroups)
#'
#' # plot all loaded features
#' for (n in 1:length(ncGTWinputs))
#'   plotGroup(ncGTWinputs[[n]], xcmsLargeWin@rt$raw)
#' @export

plotGroup <- function(ncGTWinput, sampleRt, sampleInd = 1:dim(ncGTWinput$rtRaw)[1],
                      ind = NULL, savePath = NULL, show = TRUE, sub = TRUE, filter = FALSE){
  samNum <- dim(ncGTWinput$rtRaw)[1]
  profiles <- ncGTWinput$profiles
  if (filter)
    for (n in 1:samNum)
      profiles[n, ] <- gaussFilter(profiles[n, ])

  rtRange <- matrix(0, samNum, dim(ncGTWinput$rtRaw)[2])
  for (n in 1:samNum)
    rtRange[n, ] <- sampleRt[[n]][ncGTWinput$rtRaw[n, ]]

  colVec <- matrix(0, samNum, 1)
  for (n in 1:dim(colVec)[1])
    colVec[n] <- rgb(n / samNum, 1 - n / samNum,
                       abs(samNum / 2 - abs(n - samNum / 2)) / (samNum / 2))

  mzmed <- round(ncGTWinput$groupInfo['mzmed'], 2)
  groupInd <- ncGTWinput$groupInfo['index']
  tit <- paste("Extracted Ion Chromatogram:", mzmed, "m/z")
  if (!is.null(ind)){
    subt <- paste0("Group ", groupInd, " (", ind, ")", "   Color: Green -> Purple -> Red")
  } else{
    subt <- paste0("Group ", groupInd,"   Color: Green -> Purple -> Red")
  }

  if (show){
    matplot(t(rtRange[sampleInd, , drop = FALSE]), t(profiles[sampleInd, , drop = FALSE]),
            type = 'l', col = colVec[sampleInd], lty = 1,
            xlab="rt (seconds)", ylab="Intensity")
    if (sub){
      title(main = paste("Extracted Ion Chromatogram:", mzmed, "m/z"), sub = subt)
    } else{
      title(main = paste("Extracted Ion Chromatogram:", mzmed, "m/z"))
    }
  }
  if (length(savePath)!=0){
    if (substr(savePath, nchar(savePath), nchar(savePath)) != '/')
      savePath <- paste0(savePath, '/')
    if (ind){
      filePath <- paste0(savePath, "group", groupInd, "_", ind, ".png")
    } else{
      filePath <- paste0(savePath, "group", groupInd, ".png")
    }
    png(filename = filePath, width = 1080, height = 720)
    matplot(t(rtRange[sampleInd, , drop = FALSE]), t(profiles[sampleInd, , drop = FALSE]),
            type = 'l', col = colVec[sampleInd], lty = 1,
            xlab="rt (seconds)", ylab="Intensity")
    if (sub){
      title(main = paste("Extracted Ion Chromatogram:", mzmed, "m/z"), sub = subt)
    } else{
      title(main = paste("Extracted Ion Chromatogram:", mzmed, "m/z"))
    }
    dev.off()
  }
}
