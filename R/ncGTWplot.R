plotGroup <- function(ncGTWinput, sampleRt, sampleInd = 1:dim(ncGTWinput$rtRaw)[1]){
  samNum <- dim(ncGTWinput$rtRaw)[1]
  profiles <- ncGTWinput$profiles
  rtRange <- matrix(0, samNum, dim(ncGTWinput$rtRaw)[2])
  for (n in 1:samNum)
    rtRange[n, ] <- sampleRt[[n]][ncGTWinput$rtRaw[n, ]]
  
  colVec <- matrix(0, samNum, 1)
  for (n in 1:dim(colVec)[1])
    colVec[n] <- rgb(n / samNum, 1 - n / samNum, 
                       abs(samNum / 2 - abs(n - samNum / 2)) / (samNum / 2))
  matplot(t(rtRange[sampleInd, , drop = FALSE]), t(profiles[sampleInd, , drop = FALSE]), 
          type = 'l', col = colVec[sampleInd], lty = 1,
          xlab="rt (seconds)", ylab="Intensity")
  mzmed <- round(ncGTWinput$groupInfo['mzmed'], 2)
  title(main = paste("Extracted Ion Chromatogram:", mzmed, "m/z"),
        sub = "Color: Red -> Purple -> Green")
}
