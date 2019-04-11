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
