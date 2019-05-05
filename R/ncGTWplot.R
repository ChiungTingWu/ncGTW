#' Plot profiles for each peak group
#'
#' This function calculates the p-value of each peak group in the
#' \code{\link[xcms::xcmsSet-class]{xcmsSet}} with the smaller "bw" parameter,
#' and finds the corresponding peak group in the
#' \code{\link[xcms::xcmsSet-class]{xcmsSet}} with the larger "bw" parameter.
#' @param xcmsLargeWin A \code{\link[xcms::xcmsSet-class]{xcmsSet}} object with
#'     a larger bw, usually the maximum expected retension time drift.
#' @param xcmsSmallWin A \code{\link[xcms::xcmsSet-class]{xcmsSet}}
#'     object with a smaller bw, usually the resolution of the retension time.
#' @details This function includes two major steps to determine a peak group is
#' misaligned or not.
#' @return A \code{\link[xcms::xcmsSet-class]{xcmsSet}} object with all
#' detected misaligned peak groups.
#' @examples
#' add(1, 1)
#' add(10, 1)
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
  if (ind){
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
