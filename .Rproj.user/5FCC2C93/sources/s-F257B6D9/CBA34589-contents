


parInfos <- vector('list', parSect)
for (n in 1:parSect){
  parInfos[[n]]$num <- n
  parInfos[[n]]$parSpec <- parSpec[1:parNum[n], , n]
  parInfos[[n]]$parNum <- parNum[n]
  parInfos[[n]]$parInd <- parInd[ , n]
}
parInfo = parInfos[[1]]



ncGTW1stOutput = bplapply(parInfos, ncGTW1stLayer, parSect, xcmsLargeWin, groupInd, scanRange,
         mir, strNum, diaNum, noiseVar, maxStp, dia, logt, nor,
         mu, sigma, biP, weiP, rangeThre, BPPARAM = SnowParam(7))



matplot(t(ncGTW1stOutput[[1]]$parWarped), type = 'l')
matplot(t(parWarped[[1]]), type = 'l')

ncGTW1stLayer(parInfo, parSect, xcmsLargeWin, groupInd, scanRange,
              mir, strNum, diaNum, noiseVar, maxStp, dia, logt, nor,
              mu, sigma, biP, weiP, rangeThre)

