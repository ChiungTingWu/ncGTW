gaussFilter <- function(prof, sig = 1){
  sz <- ceiling(sig * 6)    # length of gaussian filter vector
  if (sz < 2){
    sz <- 2
  }
  if (sz %% 2 != 0){
    sz <- sz + 1
  }
  x <- seq(-sz / 2, sz / 2, length = sz + 1)
  filterVec <- exp(-x ^ 2 / (2 * sig ^ 2))
  filterVec <- filterVec / sum (filterVec) # normalize
  filtered <- convolve(prof, filterVec, type = 'open')
  return(filtered[(sz / 2 + 1):(length(filtered) - sz / 2)])
}

rt2scan <- function(rt, rtAll)
  return(which.min(abs(rtAll - rt)))

smoTest <- function(xcmsLargeWin, groupInd, dataSub, scanRange, sampleInd, path2){
  # groupInd = 176;
  # sampleInd = parInd[1:parNum[n], n]
  peaks <- xcmsLargeWin@peaks
  groupidx <- xcmsLargeWin@groupidx
  rtXCMS <- xcmsLargeWin@rt$corrected
  rtRaw <- xcmsLargeWin@rt$raw

  prePeaks <- peaks[groupidx[[groupInd]], ]
  prePeaks <- round(prePeaks[is.element(prePeaks[ , 'sample'], sampleInd), , drop = FALSE], digits = 4)

  if (length(prePeaks) == 0){
    return(matrix(-1, 3, 3))
  } else{
    if (dim(prePeaks)[1] > 1)
      prePeaks <- prePeaks[!duplicated(prePeaks[,c('rt', 'rtmax', 'rtmin')]), , drop = FALSE]
    prePeakInd <- prePeaks[ , 'sample']
    prePeakMed <- prePeaks[, 'rt']

    sampleCount <- table(prePeakInd)
    groupNum <- max(sampleCount)
    groupSam <- as.numeric(names(sampleCount)[which.max(sampleCount)])
    if (groupNum != 1)
      groupSam <- as.numeric(names(sampleCount)[which(sampleCount == groupNum)])

    if (length(groupSam)>1){
      maxRange <- 0
      maxInd <- 0
      for (ind in 1:length(groupSam)){
        samPeaks <- prePeaks[prePeakInd ==  groupSam[ind], 'rt']
        if (max(samPeaks) - min(samPeaks) > maxRange){
          maxRange <- max(samPeaks) - min(samPeaks)
          maxInd <- groupSam[ind]
        }
      }
      groupSam <- if (maxRange == 0) groupSam[1] else maxInd
    }

    kmeansPreInd <- kmeans(prePeaks[, c('rt', 'rtmax', 'rtmin'), drop = FALSE],
                           prePeaks[prePeakInd == groupSam, c('rt', 'rtmax', 'rtmin'), drop = FALSE])

    oriPeakGroup <- vector('list', groupNum)
    XCMSPeakGroup <- vector('list', groupNum)
    ncGTWPeakGroup <- vector('list', groupNum)

    for (n in 1:groupNum)
      XCMSPeakGroup[[n]] <- prePeakMed[kmeansPreInd$cluster == n]

    ncGTWPeakMed <- prePeakMed * 0
    oriPeakMed <- prePeakMed * 0

    for (n in 1:length(ncGTWPeakMed)){
      samInd <- prePeaks[n, 'sample']
      samSubInd <- which(sampleInd == prePeaks[n, 'sample'])

      indDif <- abs(scanRange[samInd, ] - rt2scan(prePeakMed[n], rtXCMS[[samInd]]))
      minIndDif <- min(indDif)
      medInd <- which(indDif == minIndDif)
      medInd <- medInd[which.max(dataSub[samSubInd, medInd])]

      if (medInd - 5 < 1){
        staInd <- 1
      } else {
        staInd <- medInd - 5
      }
      if (medInd + 5 > dim(dataSub)[2]){
        endInd <- dim(dataSub)[2]
      } else{
        endInd <- medInd + 5
      }
      apexRange <- staInd:endInd
      apexInd <- apexRange[which.max(dataSub[samSubInd, apexRange])]
      oriPeakMed[n] <- rtRaw[[samInd]][scanRange[samInd, apexInd]]

      samPath <- path2[[samSubInd]]
      ncGTWPeakMed[n] <- rtRaw[[samInd]][scanRange[samInd, round(mean(
        samPath[which(samPath[ , 2] == apexInd), 1]))]]

    }
    oriPeakRt <- cbind(oriPeakMed, prePeaks[, 'rtmin'] - prePeakMed + oriPeakMed,
                       prePeaks[, 'rtmax'] - prePeakMed + oriPeakMed)
    ncGTWPeakRt <- cbind(ncGTWPeakMed, prePeaks[, 'rtmin'] - prePeakMed + ncGTWPeakMed,
                       prePeaks[, 'rtmax'] - prePeakMed + ncGTWPeakMed)

    kmeansOriInd <- kmeans(oriPeakRt, oriPeakRt[prePeakInd == groupSam, , drop = FALSE])
    kmeansncGTWInd <- kmeans(ncGTWPeakRt, ncGTWPeakRt[prePeakInd == groupSam, , drop = FALSE])

    for (n in 1:groupNum){
      oriPeakGroup[[n]] <- oriPeakMed[kmeansOriInd$cluster == n]
      ncGTWPeakGroup[[n]] <- ncGTWPeakMed[kmeansncGTWInd$cluster == n]
    }

    statResult <- matrix(0, 3, 2)
    statResult[1, 1] <- sum(sapply(oriPeakGroup, var), na.rm = TRUE)
    statResult[2, 1] <- sum(sapply(XCMSPeakGroup, var), na.rm = TRUE)
    statResult[3, 1] <- sum(sapply(ncGTWPeakGroup, var), na.rm = TRUE)
    statResult[1, 2] <- max(sapply(oriPeakGroup, function(x) range(x)[2] - range(x)[1]))
    statResult[2, 2] <- max(sapply(XCMSPeakGroup, function(x) range(x)[2] - range(x)[1]))
    statResult[3, 2] <- max(sapply(ncGTWPeakGroup, function(x) range(x)[2] - range(x)[1]))
  }
  return(statResult)

}

meanCorOl <- function(ncGTWinput, sampleRt){
  samNum <- dim(ncGTWinput$rtRaw)[1]
  pointNum <- dim(ncGTWinput$rtRaw)[2]
  profiles <- ncGTWinput$profiles
  rtRange <- matrix(0, samNum, pointNum)
  for (n in 1:samNum){
    profiles[n, ] <- gaussFilter(profiles[n, ])
    rtRange[n, ] <- sampleRt[[n]][ncGTWinput$rtRaw[n, ]]
  }
  proInter <- matrix(0, samNum, pointNum * 10)
  interX <- seq(max(rtRange[ , 1]), min(rtRange[ , pointNum]), length.out = pointNum * 10)
  for (n in 1:samNum)
    proInter[n, ] <- approx(rtRange[n, ], profiles[n, ], interX, yleft = NA, yright = NA)$y
  corM <- cor(t(proInter))
  olM <- matrix(0, samNum, samNum)
  for (i in 1:samNum)
    for (j in i:samNum){
      olM[i, j] <- sum(pmin(proInter[i, ], proInter[j, ])) / min(sum(proInter[i, ]), sum(proInter[j, ]))
      olM[j, i] <- olM[i, j]
    }

  return(list(cor = mean(corM), ol = mean(olM)))

}

compCV <- function(XCMSresFilled, ncGTWresFilled){
  groupNum <- dim(XCMSresFilled@groups)[1]
  sampleNum <- max(XCMSresFilled@peaks[, 'sample'])

  ncGTWpeaks <- matrix(0, groupNum, sampleNum)
  XCMSpeaks <- matrix(0, groupNum, sampleNum)

  ncGTWcv <- matrix(0, groupNum, 1)
  XCMScv <- matrix(0, groupNum, 1)

  for (n in 1:groupNum)
  {
    ncGTWgroupPeaks <- ncGTWresFilled@peaks[ncGTWresFilled@groupidx[[n]], ]
    XCMSgroupPeaks <- XCMSresFilled@peaks[XCMSresFilled@groupidx[[n]], ]

    ncGTWonePeak <- matrix(0, sampleNum, dim(ncGTWresFilled@peaks)[2])
    colnames(ncGTWonePeak) <- colnames(ncGTWgroupPeaks)
    XCMSonePeak <- matrix(0, sampleNum, dim(XCMSresFilled@peaks)[2])
    colnames(XCMSonePeak) <- colnames(XCMSgroupPeaks)

    for (m in 1:sampleNum)
    {
      ncGTWpeakInd <- which(ncGTWgroupPeaks[, 'sample'] == m)
      ncGTWpeakInd <- ncGTWpeakInd[which.max(ncGTWgroupPeaks[ncGTWpeakInd, 'into'])]
      ncGTWonePeak[m, ] <- ncGTWgroupPeaks[ncGTWpeakInd, ]

      XCMSpeakInd <- which(XCMSgroupPeaks[, 'sample'] == m)
      XCMSpeakInd <- XCMSpeakInd[which.max(XCMSgroupPeaks[XCMSpeakInd, 'into'])]
      XCMSonePeak[m, ] <- XCMSgroupPeaks[XCMSpeakInd, ]
    }
    ncGTWpeaks[n, ] <- ncGTWonePeak[ , 'into']
    XCMSpeaks[n, ] <- XCMSonePeak[ , 'into']

    ncGTWcv[n] <- sd(ncGTWpeaks[n, ])/mean(ncGTWpeaks[n, ])
    XCMScv[n] <- sd(XCMSpeaks[n, ])/mean(XCMSpeaks[n, ])
  }
  return(list(ncGTWcv = ncGTWcv, XCMScv = XCMScv))
}
