label2path <- function(cut, gtwInfo){
  #label2path4Aosokin Convert label of src and sink to path in primal graph
  # For the graph representation of Aosokin's codes
  # Src and sink edges do not use explicit node names, but other edges do

  # cs = find(labels==0);
  ct <- which(cut == 1)

  nEdgeGrid <- gtwInfo$nEdgeGrid
  nNodeGrid <- gtwInfo$nNodeGrid
  nPix <- gtwInfo$nPix
  pEdgeSS <- gtwInfo$pEdgeSS
  pEdgeEE <- gtwInfo$pEdgeEE
  dEdgeIntSS <- gtwInfo$dEdgeIntSS
  ssTmpPos <- gtwInfo$ssTmpPos
  ss <- gtwInfo$ss
  ee <- gtwInfo$ee

  # cut for within grid
  # not suitable do this pixel by pixel due to the spatial edge
  dtmp <- matrix(0, dim(ee)[1], 2)
  ia <- matrix(is.element(ee[ ,1:2], ct), dim(ee)[1], 2)
  dtmp[ia] <- 1
  isCutEE <- dtmp[ ,1] != dtmp[ , 2]

  # cut to mapping pattern in primal graph
  resPath <- vector('list', nPix)
  for (nn in 1:nPix){
    # cuts within grid
    idx1 <- nEdgeGrid * (nn - 1) + 1
    idx2 <- nEdgeGrid * nn
    cutNow <- isCutEE[idx1:idx2]
    resEE <- pEdgeEE[cutNow, ]

    # src and sink cuts
    idx1 <- nNodeGrid * (nn - 1) + 1
    idx2 <- nNodeGrid * nn
    s0 <- ss[idx1:idx2, ]
    b0 <- cut[idx1:idx2]

    # idxSrcCut = find(s0(:,1)>0 & b0==1);  % nodes that cuts
    idxSrcCut <- which(ssTmpPos[ , 1] != (nPix * nPix + 1) & b0 == TRUE)
    ia <- is.element(dEdgeIntSS[ , 2], idxSrcCut)  # get corresponding edges, src -> node
    resSrc <- pEdgeSS[ia, ]  # the primal edges

    # idxSinkCut = find(s0(:,2)>0 & b0==0);
    idxSinkCut <- which(ssTmpPos[ , 2] != nPix * nPix + 1 & b0 == 0)
    ia <- is.element(dEdgeIntSS[ , 1], idxSinkCut)  # node -> sink
    resSink <- pEdgeSS[ia, ]

    tempPath <- rbind(resEE, resSrc, resSink)
    sorted <- sort.int(tempPath[ , 1], index.return = TRUE)

    resPath[[nn]] <- tempPath[sorted$ix, ]
  }
  return(resPath)
}

warpCurve <- function(curve, path){

  nTps <- length(curve)
  warped <- matrix(0, 1, nTps)
  p0 <- path[ , 1:2]
  idxValid <- (p0[ , 1] >= 1) & (p0[ , 1] <= nTps) & (p0[ , 2] >= 1) & (p0[ , 2] <= nTps)
  p0 <- p0[idxValid, ]

  for (tt in 1:dim(p0)[1]){
    pRef <- p0[tt, 2]
    pTst <- p0[tt, 1]
    warped[pTst] <- max(warped[pTst], curve[pRef])
  }
  return(warped)
}


pathCombine <- function(parPath, path2, parInd){
  dataNum <- max(parInd)
  temp2path <- vector('list', dataNum)

  for (tempInd in 1:dataNum){
    ind <- which(parInd == tempInd, TRUE)
    tempPath1 <- parPath[[ind[2]]][[ind[1]]]
    tempPath2 <- path2[[ind[2]]]
    newPath <- c(0,0)

    for (pixInd in 2:dim(tempPath1)[1]){
      tempPix <- tempPath1[pixInd, 1]
      tempPix2 <- tempPath2[which(tempPath2[ , 2] == tempPix)]

      tempM <- matrix(0, length(tempPix2), 2)
      tempM[ , 1] <- tempPix2
      tempM[ , 2] <- tempPath1[pixInd, 2]
      newPath <- rbind(newPath, tempM)
    }

    for (pixInd in seq(dim(newPath)[1], 2, -1)){
      if ((newPath[pixInd - 1, 1] > newPath[pixInd, 1]) ||
          (newPath[pixInd - 1, 2] > newPath[pixInd, 2])){
        newPath[pixInd, ] <- 0
      }
    }

    newPath <- cbind(newPath, rbind(newPath[2:dim(newPath)[1], ], tempPath2[dim(tempPath2)[1], 3:4]))
    temp2path[[tempInd]] <- newPath
  }
  return(temp2path)
}


adjustRT <- function(xcmsLargeWin, ncGTWinput, ncGTWoutput, ppm){
  peaks <- xcmsLargeWin@peaks
  groups <- xcmsLargeWin@groups
  groupidx <- xcmsLargeWin@groupidx
  groupInd <- ncGTWinput$groupInfo['index']
  rtXCMS <- xcmsLargeWin@rt$corrected
  rtRaw <- xcmsLargeWin@rt$raw
  scanRangeOld <- ncGTWinput$rtRaw
  data <- ncGTWoutput$data
  parPath <- ncGTWoutput$parPath
  path2 <- ncGTWoutput$path2
  parInd <- ncGTWoutput$parInd
  scanRange <- ncGTWoutput$scanRange
  downSample <- ncGTWoutput$downSample

  XCMSPeaks <- cbind(peaks[groupidx[[groupInd]], ], peakInd = groupidx[[groupInd]])
  XCMSPeaks <- XCMSPeaks[!duplicated(XCMSPeaks[,c('rt', 'rtmax', 'rtmin')]), ]
  rmInd <- matrix(TRUE, dim(XCMSPeaks)[1], 1)
  for (n in 1:(dim(XCMSPeaks)[1] - 1)){
    XCMSSamPeaks <- XCMSPeaks[XCMSPeaks[ , 'sample'] == XCMSPeaks[n, 'sample'], , drop = FALSE]
    tempRmInd <- rmInd[XCMSPeaks[ , 'sample'] == XCMSPeaks[n, 'sample'], , drop = FALSE]
    if (abs(XCMSPeaks[n, 'rtmax'] - XCMSPeaks[n, 'rtmin']) > 15){
      rmInd[n] <- FALSE
    } else{
    if (!isEmpty(XCMSSamPeaks))
      for (m in 1:dim(XCMSSamPeaks)[1]){
        if (abs(XCMSPeaks[n, 'mz'] - XCMSSamPeaks[m, 'mz']) / XCMSPeaks[n, 1] < ppm/1000000 &&
          abs(XCMSPeaks[n, 'rt'] - XCMSSamPeaks[m, 'rt']) < 5){
          if (XCMSPeaks[n, 'maxo'] < XCMSSamPeaks[m, 'maxo'] && tempRmInd[m] == TRUE){
            rmInd[n] <- FALSE
            break
          }
        }

      }
    }
  }
  XCMSPeaks <- XCMSPeaks[rmInd, ]

  sampleCount <- table(XCMSPeaks[ , 'sample'])
  groupNum <- max(sampleCount)

  ncGTWpeaks <- XCMSPeaks
  ncGTWpeaks[ , c('rt', 'rtmin', 'rtmax')] <- 0
  peakScanRt <- matrix(0, dim(ncGTWpeaks)[1], 10)
  colnames(peakScanRt) <- c('scan_XCMS', 'scan_ncGTW', 'rt_ncGTW', 'scanmin_XCMS',
        'scanmin_ncGTW', 'rtmin_ncGTW', 'scanmax_XCMS', 'scanmax_ncGTW', 'rtmax_ncGTW', 'sample')
  path <- pathCombine(parPath, path2, parInd)

  warpOrder <- sort(XCMSPeaks[ , 'maxo'], index.return = TRUE, decreasing = TRUE)$ix
  for (nn in 1:dim(ncGTWpeaks)[1]){
    n <- warpOrder[nn]
    samInd <- XCMSPeaks[n, 'sample']
    samPath <- path[[samInd]]
    samRtXCMS <- rtXCMS[[samInd]][scanRange[samInd, ]]
    samRtRaw <- rtRaw[[samInd]][scanRange[samInd, ]]
    scanSubXCMS <- rt2scan(XCMSPeaks[n, 'rt'], samRtXCMS)
    scanSubncGTW <- round(mean(samPath[which(samPath[ ,2] == scanSubXCMS), 1]))
    ncGTWrt <- rtRaw[[samInd]][scanRange[samInd, scanSubncGTW]]

    peakScanRt[n, c('scan_XCMS', 'scan_ncGTW', 'rt_ncGTW')] <- c(scanSubXCMS, scanSubncGTW, ncGTWrt)
    inSamPeaks <- cbind(peakScanRt[peakScanRt[ , 'sample'] == samInd, , drop = FALSE],
                        peakInd = which(peakScanRt[ , 'sample'] == samInd))
    if (length(inSamPeaks) > 0)
      for (m in 1:dim(inSamPeaks)[1])
        if (XCMSPeaks[n, 'rt'] > XCMSPeaks[inSamPeaks[m, 'peakInd'], 'rt']){
          peakRtDif <- min(XCMSPeaks[n, 'rt'] - XCMSPeaks[inSamPeaks[m, 'peakInd'], 'rt'],
                           XCMSPeaks[n, 'rt'] - XCMSPeaks[n, 'rtmin'] +
                           XCMSPeaks[inSamPeaks[m, 'peakInd'], 'rtmax'] -
                           XCMSPeaks[inSamPeaks[m, 'peakInd'], 'rt'])
          if (peakScanRt[n, 'rt_ncGTW'] - inSamPeaks[m, 'rt_ncGTW'] < peakRtDif){
            peakScanRt[n, 'scan_ncGTW'] <- rt2scan(inSamPeaks[m, 'rt_ncGTW'] + peakRtDif, samRtXCMS)
            peakScanRt[n, 'rt_ncGTW'] <- samRtRaw[peakScanRt[n, 'scan_ncGTW']]
          }
        }
    peakScanRt[n, 'rtmin_ncGTW'] <- XCMSPeaks[n, 'rtmin'] - XCMSPeaks[n, 'rt'] + peakScanRt[n, 'rt_ncGTW']
    peakScanRt[n, 'rtmax_ncGTW'] <- XCMSPeaks[n, 'rtmax'] - XCMSPeaks[n, 'rt'] + peakScanRt[n, 'rt_ncGTW']
    peakScanRt[n, 'scanmin_ncGTW'] <- rt2scan(peakScanRt[n, 'rtmin_ncGTW'], samRtRaw)
    peakScanRt[n, 'scanmax_ncGTW'] <- rt2scan(peakScanRt[n, 'rtmax_ncGTW'], samRtRaw)
    peakScanRt[n, 'scanmin_XCMS'] <- peakScanRt[n, 'scanmin_ncGTW'] -
      peakScanRt[n, 'scan_ncGTW'] + peakScanRt[n, 'scan_XCMS']
    peakScanRt[n, 'scanmax_XCMS'] <- peakScanRt[n, 'scanmax_ncGTW'] -
      peakScanRt[n, 'scan_ncGTW'] + peakScanRt[n, 'scan_XCMS']
    peakScanRt[n, 'sample'] <- XCMSPeaks[n, 'sample']
    ncGTWpeaks[n, c('rt', 'rtmin', 'rtmax')] <- peakScanRt[n, c('rt_ncGTW', 'rtmin_ncGTW', 'rtmax_ncGTW')]
    peaks[ncGTWpeaks[n, 'peakInd'], c('rt', 'rtmin', 'rtmax')] <- ncGTWpeaks[n, c('rt', 'rtmin', 'rtmax')]

  }

  groupSam <- as.numeric(names(which(sampleCount == groupNum)))
  if (length(groupSam)>1){
    maxRange <- 0
    maxInd <- 0
    for (ind in 1:length(groupSam)){
      samPeaks <- XCMSPeaks[XCMSPeaks[ , 'sample'] ==  groupSam[ind], 'rt']
      if (max(samPeaks) - min(samPeaks) > maxRange){
        maxRange <- max(samPeaks) - min(samPeaks)
        maxInd <- groupSam[ind]
      }
    }
    groupSam <- if (maxRange == 0) groupSam[1] else maxInd
  } else groupSam <- groupSam[1]

#  if (groupNum > 4)
#    groupNum <- 4

  kcenter <- peakScanRt[peakScanRt[ , 'sample'] == groupSam,
                        c('rt_ncGTW', 'rtmin_ncGTW', 'rtmax_ncGTW'), drop = FALSE]
  intOrder <- sort(XCMSPeaks[peakScanRt[ , 'sample'] == groupSam, 'maxo'],
                    index.return = TRUE, decreasing = TRUE)$ix[1:groupNum]
  kcenter <- kcenter[intOrder, , drop = FALSE]


  kmeansncGTW <- kmeans(peakScanRt[ , c('rt_ncGTW', 'rtmin_ncGTW', 'rtmax_ncGTW')], kcenter)

  peakGroup <- vector('list', groupNum)
  tempPeaks <- matrix(NA, max(peaks[ , 'sample']), dim(peakScanRt)[2] + 2)
  colnames(tempPeaks) <- c(colnames(peakScanRt), 'detect', 'maxo')
  colnames(tempPeaks)[c(1, 4, 7)] <- c('scan_Raw', 'scanmin_Raw', 'scanmax_Raw')
  for (n in 1:groupNum)
    peakGroup[[n]] <- tempPeaks

  for (n in 1:dim(peakScanRt)[1])
    if (is.na(peakGroup[[kmeansncGTW$cluster[n]]][peakScanRt[n, 'sample'], 1])){
      peakGroup[[kmeansncGTW$cluster[n]]][peakScanRt[n, 'sample'], ] <-
        c(peakScanRt[n, ], 1, XCMSPeaks[n, 'maxo'])
    } else{
      if (XCMSPeaks[n, 'maxo'] > peakGroup[[kmeansncGTW$cluster[n]]][peakScanRt[n, 'sample'], 'maxo'])
        peakGroup[[kmeansncGTW$cluster[n]]][peakScanRt[n, 'sample'], ] <-
          c(peakScanRt[n, ], 1, XCMSPeaks[n, 'maxo'])
    }

  rmInd <- matrix(TRUE, length(peakGroup), 1)
  for (n in 1:groupNum){
    tempPeaks <- peakGroup[[n]]
    if (sum(is.na(tempPeaks[,1])) > dim(tempPeaks)[1] * (3 / 4) )
      rmInd[n] <- FALSE
  }
  peakGroup <- peakGroup[rmInd]
  groupNum <- length(peakGroup)
  peakGroupMed <- tempPeaks[1:groupNum, , drop = FALSE]
  for (n in 1:groupNum)
    for (m in 1:dim(peakGroupMed)[2])
      peakGroupMed[n, m] <- median(peakGroup[[n]][ , m], TRUE)

  for (n in 1:groupNum){
    tempPeaks <- peakGroup[[n]]
    for (m in 1:dim(tempPeaks)[1])
      if (is.na(tempPeaks[m, 1])){
        tempPeaks[m , c('scan_ncGTW', 'rt_ncGTW', 'scanmin_ncGTW',
                       'rtmin_ncGTW', 'scanmax_ncGTW', 'rtmax_ncGTW')] <-
          peakGroupMed[n, c('scan_ncGTW', 'rt_ncGTW', 'scanmin_ncGTW',
                            'rtmin_ncGTW', 'scanmax_ncGTW', 'rtmax_ncGTW')]
        samPath <- path[[m]]
        apexInd <- samPath[which(samPath[ , 1] == round(peakGroupMed[n, 'scan_ncGTW'])), 2]
        apexInd <- apexInd[which.max(data[n, apexInd])]

        tempPeaks[m, 'scan_Raw'] <- apexInd
        tempPeaks[m, 'scanmin_Raw'] <- peakGroupMed[n, 'scanmin_ncGTW'] +
          apexInd - peakGroupMed[n, 'scan_ncGTW'];
        tempPeaks[m, 'scanmax_Raw'] <- peakGroupMed[n, 'scanmax_ncGTW'] +
          apexInd - peakGroupMed[n, 'scan_ncGTW'];
        tempPeaks[m, 'sample'] <- m

      }
    peakGroup[[n]] <- tempPeaks
  }
  rtncGTW <- rtXCMS
  rtncGTWsub <- scanRangeOld * NA
  rtRawSub <- scanRangeOld * NA
  for (n in 1:length(rtRaw))
    rtRawSub[n, ] <- rtRaw[[n]][scanRangeOld[n, ]]

  tempShift <- 5
  while (any(rtRawSub[ ,2:dim(rtRawSub)[2]] - rtRawSub[ ,1:(dim(rtRawSub)[2] - 1)] > 5)){
    for (n in 1:length(rtRaw)){
      rtRawSub[n, ] <- rtRaw[[n]][scanRangeOld[n, ] - tempShift]
      tempShift <- tempShift + 5
    }
  }


  for (n in 1:groupNum){
    tempPeaks <- peakGroup[[n]]
    for (m in 1:dim(tempPeaks)[1]){
      rawSta <- tempPeaks[m, 'scanmin_Raw'] * downSample
      rawEnd <- tempPeaks[m, 'scanmax_Raw'] * downSample
      ncGTWSta <- tempPeaks[m, 'scanmin_ncGTW'] * downSample
      ncGTWEnd <- tempPeaks[m, 'scanmax_ncGTW'] * downSample

      if (rawEnd > dim(scanRangeOld)[2]){
        ncGTWEnd <- ncGTWEnd - (rawEnd - dim(scanRangeOld)[2])
        rawEnd <- dim(scanRangeOld)[2]
      }
      if (ncGTWEnd > dim(rtRawSub)[2]){
        rawEnd <- rawEnd - (ncGTWEnd - dim(rtRawSub)[2])
        ncGTWEnd <- dim(rtRawSub)[2]
      }
      if (rawSta < 1){
        ncGTWSta <- ncGTWSta + (1 - rawSta)
        rawSta <- 1
      }
      rawInd <- rawSta:rawEnd
      overlapFlag <- !is.na(rtncGTWsub[m, rawInd])
      rtncGTWsub[m, rawInd] <- rtRawSub[m, ncGTWSta:ncGTWEnd]
      if (sum(overlapFlag) != 0){
        olInd <- which(overlapFlag)
        olSta <- rawInd[olInd[1]] - 1
        if (olSta == 0)
          olSta <- 1
        olEnd <- rawInd[olInd[length(olInd)]] + 1
        if (olEnd > dim(scanRangeOld)[2])
          olEnd <- dim(scanRangeOld)[2]
        while (is.na(rtncGTWsub[m, olSta]))
          olSta <- olSta + 1
        while (is.na(rtncGTWsub[m, olEnd]))
          olEnd <- olEnd - 1

        while (rtncGTWsub[m, olSta] >= rtncGTWsub[m, olEnd]){
          olStaOld <- olSta
          olEndOld <- olEnd
          if ((olSta - 1 > 0) && !is.na(rtncGTWsub[m, olSta - 1]))
            olSta <- olSta - 1
          if((olEnd + 1 < dim(rtncGTWsub)[2]) && !is.na(rtncGTWsub[m, olEnd + 1]))
            olEnd <- olEnd + 1
          if (olSta == 0)
            olSta <- 1
          if (olEnd > dim(scanRangeOld)[2])
            olEnd <- dim(scanRangeOld)[2]
          if (olStaOld == olSta && olEndOld == olEnd){
            temp <- rtncGTWsub[m, olSta]
            rtncGTWsub[m, olSta] <- rtncGTWsub[m, olEnd]
            rtncGTWsub[m, olEnd] <- temp
          }


        }
        rtncGTWsub[m, olSta:olEnd] <-
          approx(c(olSta, olEnd), rtncGTWsub[m, c(olSta, olEnd)], n = length(olSta:olEnd))$y
      }
    }
  }

  for (n in 1:dim(rtncGTWsub)[1]){
    while (!is.na(which(is.na(rtncGTWsub[n, ]))[1])){
      ipSta <- which(is.na(rtncGTWsub[n, ]))[1]
      ipEnd <- which(!is.na(rtncGTWsub[n, ipSta:dim(rtncGTWsub)[2]]))[1]
      if (is.na(ipEnd)){
        ipEnd <- dim(rtncGTWsub)[2]
      } else {
        ipEnd <- ipEnd + ipSta - 2
      }
      if (ipSta == 1){
        ipStaRt <- rtncGTW[[n]][scanRangeOld[n, ipSta] - 1]
      } else{
        ipStaRt <- rtncGTWsub[n, ipSta - 1]
      }
      if (ipEnd == dim(rtncGTWsub)[2]){
        ipEndRt <- rtncGTW[[n]][scanRangeOld[n, ipEnd] + 1]
      } else{
        ipEndRt <- rtncGTWsub[n, ipEnd + 1]
      }

      if (ipStaRt > ipEndRt){
        if (ipSta == 1){
          nearInd <- rt2scan(ipEnd, rtncGTW[[n]])
          if (rtncGTW[[n]][nearInd] - ipEndRt > 0)
            nearInd <- nearInd - 1
          ipStaRt <- rtncGTW[[n]][ipSta - ipEnd + nearInd]
        }
      } else{
        while (ipStaRt > ipEndRt){
          ipSta <- ipSta - 1
          if (ipSta < 2)
            ipSta <- 2
          if (is.na(rtncGTWsub[n, ipSta - 1]))
            ipSta <- ipSta + 1
          ipStaRt <- rtncGTWsub[n, ipSta - 1]

          ipEnd <- ipEnd + 1
          if (ipEnd > dim(rtncGTWsub)[2] - 1)
            ipEnd <- dim(rtncGTWsub)[2] - 1
          if (is.na(rtncGTWsub[n, ipEnd + 1]))
            ipEnd <- ipEnd - 1
          ipEndRt <- rtncGTWsub[n, ipEnd + 1]
        }
      }

      ipRt <- approx(c(ipSta - 1, ipEnd + 1), c(ipStaRt, ipEndRt),
                     n = length((ipSta - 1):(ipEnd + 1)))$y
      rtncGTWsub[n, ipSta:ipEnd] <- ipRt[2:(length(ipRt) - 1)]

      if (ipSta == 1){
        ipEndTemp <- rt2scan(ipEndRt, rtncGTW[[n]])
        if (ipEndRt < rtncGTW[[n]][ipEndTemp])
          ipEndTemp <- ipEndTemp - 1
        ipStaTemp <- ipEndTemp + ipSta - ipEnd
        if (ipStaTemp < 0){
          ipEndTemp <- ipEndTemp + 1 - ipStaTemp
          ipStaTemp <- 1
        }
        rtncGTWsub[n, ipSta:ipEnd] <- rtncGTW[[n]][ipStaTemp:ipEndTemp]
      }

      if (ipEnd == dim(rtncGTWsub)[2]){
        ipStaTemp <- rt2scan(ipStaRt, rtncGTW[[n]])
        if (ipStaRt > rtncGTW[[n]][ipStaTemp])
          ipStaTemp <- ipStaTemp - 1
        ipEndTemp <- ipStaTemp + ipEnd - ipSta
        rtncGTWsub[n, ipSta:ipEnd] <- rtncGTW[[n]][ipStaTemp:ipEndTemp]
      }
    }
  }
  for (n in 1:dim(rtncGTWsub)[1])
    for (m in 2:dim(rtncGTWsub)[2])
      if (rtncGTWsub[n, m] < rtncGTWsub[n, m - 1]){
#        warning('Adjusted rt of sample ', n, ' is not increasing...')
        pl <- m - 1
        pr <- m
        while (rtncGTWsub[n, pl] > rtncGTWsub[n, pr]){
          if (pl > 1)
            pl <- pl - 1
          if (pr < dim(rtncGTWsub)[2])
            pr <- pr + 1
        }
        rtncGTWsub[n, pl:pr] <- approx(c(pl, pr), rtncGTWsub[n, c(pl, pl)], n = length(pl:pr))$y
      }

  for (n in 1:dim(rtncGTWsub)[1])
    rtncGTW[[n]][scanRangeOld[n,]] <- rtncGTWsub[n, ]

  return(list(rtncGTW = rtncGTW, peaks = peaks))
}
