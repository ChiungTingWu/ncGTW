ncGTW <- function(ncGTWinput, xcmsLargeWin, parSamp = 10, bpParam = bpparam(), ncGTWparam = NULL){

  downSample <- if (is.null(ncGTWparam$downSample)) 2 else ncGTWparam$downSample
  nor <- if (is.null(ncGTWparam$nor)) 1 else ncGTWparam$nor
  strNum <- if (is.null(ncGTWparam$strNum)) 1 else ncGTWparam$strNum
  diaNum <- if (is.null(ncGTWparam$diaNum)) 1 else ncGTWparam$diaNum

  mir <- if (is.null(ncGTWparam$mir)) TRUE else ncGTWparam$mir


  maxStpRat <- if (is.null(ncGTWparam$maxStpRat)) 0.6 else ncGTWparam$maxStpRat
  maxStp <- if (is.null(ncGTWparam$maxStp)) round(dim(ncGTWinput$rtRaw)[2] %/%
                                                    downSample * maxStpRat) else ncGTWparam$maxStp
  rangeThre <- if (is.null(ncGTWparam$rangeThre)) 1 else ncGTWparam$rangeThre


  biP <- if (is.null(ncGTWparam$biP)) TRUE else ncGTWparam$biP
  mu <- if (is.null(ncGTWparam$mu)) 0 else ncGTWparam$mu
  sigma <- if (is.null(ncGTWparam$sigma)) 1 else ncGTWparam$sigma
  weiP <- if (is.null(ncGTWparam$weiP)) 0 else ncGTWparam$weiP
  logt <- if (is.null(ncGTWparam$logt)) 0 else ncGTWparam$logt
  dia <- if (is.null(ncGTWparam$dia)) 0 else ncGTWparam$dia
  noiseVar <- if (is.null(ncGTWparam$noiseVar)) 1 else ncGTWparam$noiseVar


  groupInd <- ncGTWinput$groupInfo[1]
  cat("ncGTW is realigning group", groupInd, "...\n")

  dataOri <- ncGTWinput$profiles
  dataFil <- t(apply(dataOri, 1, gaussFilter))
  data <- dataFil[ , seq(1, dim(dataFil)[2], downSample)]
  scanRange <- ncGTWinput$rtRaw[ , seq(1, dim(dataFil)[2], downSample)]

  dataNum <- dim(data)[1]
  specLen <- dim(data)[2]

  parSect <- dataNum / parSamp
  if (parSect > floor(dataNum / parSamp)){
    parSect <- floor(dataNum / parSamp) + 1
  }

  parSpec <- array(0, c(parSamp, specLen, parSect))
  parInd <- matrix(0, parSamp, parSect)

  parNum <- matrix(0, parSect, 1) + parSamp
  if ((dataNum - (parSect - 1) * parSamp  ) > parSamp / 2){
    parNum[parSect] <- dataNum - (parSect - 1) * parSamp
  } else{
    parNum[parSect - 1] <- ceiling((dataNum - (parSect - 2) * parSamp) / 2)
    parNum[parSect] <- floor((dataNum - (parSect - 2) * parSamp) / 2)
  }

  tempCount <- 0
  for (n in 1:parSect){
    parSpec[1:parNum[n], , n] <- data[(tempCount + 1):(tempCount + parNum[n]), ]
    parInd[1:parNum[n], n] <- (tempCount + 1):(tempCount + parNum[n])
    tempCount <- tempCount + parNum[n]
  }

  parInfos <- vector('list', parSect)
  for (n in 1:parSect){
    parInfos[[n]]$num <- n
    parInfos[[n]]$parSpec <- parSpec[1:parNum[n], , n]
    parInfos[[n]]$parNum <- parNum[n]
    parInfos[[n]]$parInd <- parInd[ , n]
  }

  ncGTW1stOutput <- bplapply(parInfos, ncGTW1stLayer, parSect, xcmsLargeWin, groupInd, scanRange,
                            mir, strNum, diaNum, noiseVar, maxStp, dia, logt, nor,
                            mu, sigma, biP, weiP, rangeThre, BPPARAM = bpParam)

  parPath <- replicate(parSect, list(vector('list', parSamp)))
  parWarped <- array(0, c(parSamp, specLen, parSect))
  parStat <- array(0, c(3, 3, parSect))

  for (n in 1:parSect){
    parPath[[n]] <- ncGTW1stOutput[[n]]$parPath
    parWarped[1:parNum[n] , , n] <- ncGTW1stOutput[[n]]$parWarped
    parStat[ , , n] <- ncGTW1stOutput[[n]]$parStat
  }

  superSample <- matrix(0, parSect, specLen)


  if (parSect == 1){
    warning('Only one layer!!!')
    path2 <- list(cbind(0:specLen, 0:specLen, 1:(specLen+1), 1:(specLen+1)))
    warpedAll <- parWarped[ , , 1]
    smoM2 <- matrix(0, 3, 6)
  }

  for (n in 1:parSect){
    maxInd <- which.max(rowSums(parWarped[ , , n], 2))
    minInd <- which.min(rowSums(parWarped[ , , n], 2))
    superSample[n, ] <- colMeans(parWarped[setdiff(1:parNum[n], c(maxInd, minInd)), , n, drop = FALSE])
    # super_sample(n, :) = super_sample(n, :)/max(super_sample(n, :));
  }


  if (maxStp * 1.5 > specLen){
    maxStp <- specLen - 1
  } else{
    maxStp <- round(maxStp * 1.5)
  }


  if (parSect > 1){
    # fprintf('%-80s%20s\n', 'Starting SAMA top layer.....', datestr(datetime()));
    parSuperWarped <- matrix(0, parSect, specLen)
    gtwPrep <- buildMultiParaValidmap(superSample, FALSE, 1, 1, biP);
    ref <- gtwPrep$ref
    tst <- gtwPrep$tst
    validMap <- gtwPrep$validMap
    sampNum <- dim(tst)[1]

    param <- initGtwParam(validMap, noiseVar, maxStp, 10^-32, dia, logt, nor)
    gtwInfo <- buildGTWgraph(ref, tst, validMap, param, mu, sigma, biP, weiP)
    eeBet <- gtwInfo$ee[(gtwInfo$nEdgeGrid * gtwInfo$nPix + 1) : dim(gtwInfo$ee)[1], ]

    param <- initGtwParam(validMap, noiseVar, maxStp, 0, dia, logt, nor)
    gtwInfo0 <- buildGTWgraph(ref, tst, validMap, param, mu, sigma, biP, weiP)
    cut0 <- graphCut(gtwInfo0$ss, gtwInfo0$ee)
    costMin <- cut0[length(cut0)]
    s1BrokenNum <- sum((cut0[eeBet[ ,1]] + cut0[eeBet[ ,2]]) == 1)
    s1BrokenNumFinal <- s1BrokenNum

    param <- initGtwParam(validMap, noiseVar, maxStp, 100000000000, dia, logt, nor)
    gtwInfoMerge <- buildGTWgraph(ref, tst, validMap, param, mu, sigma, biP, weiP, 2)
    cutMerge <- graphCut(gtwInfoMerge$ss, gtwInfoMerge$ee)
    costMax <- cutMerge[length(cutMerge)]
    smoTemp <- (costMax - costMin) / s1BrokenNum
    s1Final <- smoTemp

    # fprintf('%-80s%20s\n', ...
    #     sprintf('Obtained the range of the total cost.'), datestr(datetime()));


    smoM2 <- matrix(0, 3, 6) + 999
    smoM2[1, 1] <- smoTemp
    tempPath <- vector('list', dim(smoM2)[1])


    for (smo in 1:dim(smoM2)[1]){
      param <- initGtwParam(validMap, noiseVar, maxStp, smoM2[smo, 1], dia, logt, nor)
      # fprintf('%-80s%20s\n', ...
      # sprintf('Built the structure of GTW 1st stage.'), datestr(datetime()));
      gtwInfo <- buildGTWgraph(ref, tst, validMap, param, mu, sigma, biP, weiP)
      #     fprintf('%-80s%20s\n', ...
      #         sprintf('Converted to the 1st maximum flow.'), datestr(datetime()));
      cut <- graphCut(gtwInfo$ss, gtwInfo$ee)
      path1 <- label2path(cut, gtwInfo)
      #     fprintf('%-80s%20s\n', ...
      #         sprintf('Got the path.'), datestr(datetime()));

      sampNum2 <- dim(superSample)[1]
      #    zero_line = zeros([size(samp_num2,1)-1,spec_len]);
      validMap2 <- matrix(0, sampNum2, sampNum2)
      for (jj in 1:sampNum2)
        for (ii in seq(jj + 1, sampNum2, length = max(0, sampNum2 - jj)))
          validMap2[ii, jj] <- 1

      maxStp2 <- maxStp
      tst2 <- matrix(0, dim(superSample)[1], dim(superSample)[2])
      ref2 <- tst2

      smoothnessV2 <- 0.00001 / noiseVar
      param2 <- initGtwParam(validMap2, noiseVar, maxStp2, smoothnessV2, dia, logt, nor)
      gtwInfo2 <- buildGTWgraph(ref2, tst2, validMap2, param2, mu, sigma, biP, weiP, 1, path1, gtwPrep$pairMap)
      cut2 <- graphCut(gtwInfo2$ss, gtwInfo2$ee)
      s2BrokenNum <- cut2[length(cut2)] / smoothnessV2
      s2BrokenNumFinal <- s2BrokenNum


      s2Final <- 2 * (maxStp2 - 1) * sampNum2 / s2BrokenNum
      smoothnessV2 <- s2Final / noiseVar

      #     fprintf('%-80s%20s\n', ...
      #         sprintf('Obtained the range of the total cost.'), datestr(datetime()));

      #     mpfv_v = zeros(100, 1);

      smo2 <- 50

      smoothnessV2 <- s2Final / noiseVar * smo2 / 50
      param2 <- initGtwParam(validMap2, noiseVar, maxStp2, smoothnessV2, dia, logt, nor)
      gtwInfo2 <- buildGTWgraph( ref2, tst2, validMap2, param2,mu, sigma, biP, weiP, 1, path1, gtwPrep$pairMap)
      # fprintf('%-80s%20s\n', ...
      # sprintf('Converted to the 2nd maximum flow.'), datestr(datetime()));

      cut2 <- graphCut(gtwInfo2$ss, gtwInfo2$ee)

      path2 <- label2path(cut2, gtwInfo2)
      # fprintf('%-80s%20s\n', ...
      # sprintf('Solved the 2nd maximum flow.'), datestr(datetime()));

      tempPath[[smo]] <- path2
      temp2path <- pathCombine(parPath, path2, parInd)

      statResult <- smoTest(xcmsLargeWin, groupInd, data, scanRange, 1:dataNum, temp2path)
      tempRange <- statResult[3, 2]
      # tempVar = statResult[3, 3]
      # tempRange = 999
      # tempVar = 999

      brokenNum1 <- sum((cut[eeBet[ ,1]] + cut[eeBet[ ,2]]) == 1)
      smoM2[smo, 2] <- brokenNum1
      smoM2[smo, 3] <- tempRange
      # smoM2[smo, 4] <- tempVar

      if (tempRange <= rangeThre)
        break

      if (smo == 1){
        smoM2[2, 1] <- (costMax - cut[length(cut)]) / brokenNum1 + smoTemp
        tempCost <- cut[length(cut)] - smoTemp * brokenNum1
        smoM2[3, 1] <- (tempCost - costMin) / (s1BrokenNum - brokenNum1)
      }

      # mpfv_v(smo2, 1) = mpfv(par_super_warped, super_sample);
      #  end
    }

    bestSmo <- which.min(smoM2[ , 3])  #+ smoM2[ , 4])
    path2 <- tempPath[[bestSmo]]

    for (m in 1:parSect)
      parSuperWarped[m, ] <- warpCurve(superSample[m, ], path2[[m]])


    warpedAll <- matrix(0, dim(data)[1], dim(data)[2])
    tempCount <- 0
    for (n in 1:parSect){
      for (m in 1:parNum[n])
        warpedAll[tempCount + m, ] <- warpCurve(parWarped[m, , n], path2[[n]])
      tempCount <- tempCount + parNum[n]
    }
  }
  return(list(parPath = parPath, path2 = path2, parStat = parStat, parInd = parInd, smoM2 = smoM2,
              data = data, scanRange= scanRange, warpedAll = warpedAll, downSample = downSample))
}



ncGTW1stLayer <- function(parInfo, parSect, xcmsLargeWin, groupInd, scanRange,
                          mir, strNum, diaNum, noiseVar, maxStp, dia, logt, nor,
                          mu, sigma, biP, weiP, rangeThre){
  suppressPackageStartupMessages({require(ncGTW)})
  n <- parInfo$num
  # cat(format(Sys.time()), 'Zero', n, 'of', parSect, '\n')
  parspec <- parInfo$parSpec
  parnum <- parInfo$parNum
  parind <- parInfo$parInd
  # cat(format(Sys.time()), 'First', n, 'of', parSect, '\n')
  gtwPrep <- buildMultiParaValidmap(parspec[1:parnum, ], mir, strNum, diaNum, biP)
  ref <- gtwPrep$ref
  tst <- gtwPrep$tst
  validMap <- gtwPrep$validMap

  param <- initGtwParam(validMap, noiseVar, maxStp, 10^-32, dia, logt, nor)
  gtwInfo <- buildGTWgraph(ref, tst, validMap, param, mu, sigma, biP, weiP)
  eeBet <- gtwInfo$ee[(gtwInfo$nEdgeGrid * gtwInfo$nPix + 1) : dim(gtwInfo$ee)[1], ]

  param <- initGtwParam(validMap, noiseVar, maxStp, 0, dia, logt, nor)
  gtwInfo0 <- buildGTWgraph(ref, tst, validMap, param, mu, sigma, biP, weiP)
  cut0 <- graphCut(gtwInfo0$ss, gtwInfo0$ee)
  costMin <- cut0[length(cut0)]
  s1BrokenNum <- sum((cut0[eeBet[ ,1]] + cut0[eeBet[ ,2]]) == 1)
  parS1BrokenNum <- s1BrokenNum #####
  # cat(format(Sys.time()), 'Second', n, 'of', parSect, '\n')
  param <- initGtwParam(validMap, noiseVar, maxStp, 10^32, dia, logt, nor)
  gtwInfoMerge <- buildGTWgraph( ref, tst, validMap, param, mu, sigma, biP, weiP, 2)
  cutMerge <- graphCut(gtwInfoMerge$ss, gtwInfoMerge$ee)
  costMax <- cutMerge[length(cutMerge)]
  smoTemp <- (costMax - costMin) / s1BrokenNum


  smoM <- matrix(999, 3, 3)
  smoM[1, 1] <- smoTemp
  tempPath <- vector('list', dim(smoM)[1])

  # cat(format(Sys.time()), 'Third', n, 'of', parSect, '\n')
  for (smo in 1:dim(smoM)[1]){
    param <- initGtwParam(validMap, noiseVar, maxStp, smoM[smo, 1], dia, logt, nor)
    cat(format(Sys.time()), 'Built the structure of ncGTW stage 1, section', n, 'of', parSect, '\n')

    gtwInfo <- buildGTWgraph(ref, tst, validMap, param, mu, sigma, biP, weiP)
    cat(format(Sys.time()), 'Converted to the 1st maximum flow, section', n, 'of', parSect, '\n')

    cut <- graphCut(gtwInfo$ss, gtwInfo$ee)
    cat(format(Sys.time()), 'Solved the 1st maximum flow, section', n, 'of', parSect, '\n')


    path1 <- label2path(cut, gtwInfo)

    #fprintf('%-80s%20s\n', ...
    #        sprintf('Got the path, section %d of %d.', n, par_sect), datestr(datetime()));

    #if n == par_sect
    #fprintf('%-80s%20s\n', ...
    #        sprintf('Initializing GTW 2nd stage, section %d of %d.', n, par_sect), datestr(datetime()));
    #end

    # zero_line = zeros([size(par_t_M,1)-1,spec_len]);
    validMap2 <- matrix(0, parnum, parnum)
    for (jj in 1:parnum)
      for (ii in seq((jj + 1), parnum, length = max(0, parnum - jj)))
        validMap2[ii, jj] <- 1

    maxStp2 <- maxStp
    tst2 <- matrix(0, dim(parspec[1:parnum, ])[1], dim(parspec[1:parnum, ])[2])
    ref2 <- tst2

    smoothnessV2 <- 0.0000000001 / noiseVar
    param2 <- initGtwParam(validMap2, noiseVar, maxStp2, smoothnessV2, dia, logt, nor)
    gtwInfo2 <- buildGTWgraph(ref2, tst2, validMap2, param2, mu, sigma, biP, weiP, 1, path1, gtwPrep$pairMap)
    cut2 <- graphCut(gtwInfo2$ss, gtwInfo2$ee)

    s2BrokenNum <- cut2[length(cut2)] / smoothnessV2
    parS2BrokenNum <- s2BrokenNum #####

    parS2 <- 2 * (maxStp2 - 1) * parnum / s2BrokenNum ###############
    smoothnessV2 <- parS2 / noiseVar


    #fprintf('%-80s%20s\n', ...
    #        sprintf('Obtained the range of the total cost, section %d of %d.', n, par_sect), datestr(datetime()));

    param2 <- initGtwParam(validMap2, noiseVar, maxStp2, smoothnessV2, dia, logt, nor)
    gtwInfo2 <- buildGTWgraph( ref2, tst2, validMap2, param2,mu, sigma, biP, weiP, 1, path1, gtwPrep$pairMap)

    #fprintf('%-80s%20s\n', ...
    #        sprintf('Converted to the 2nd maximum flow, section %d of %d.', n, par_sect), datestr(datetime()));
    cut2 <- graphCut(gtwInfo2$ss, gtwInfo2$ee)

    #        if n == par_sect
    #fprintf('%-80s%20s\n', ...
    #        sprintf('Solved the 2nd maximum flow, section %d of %d.', n, par_sect), datestr(datetime()));
    #        end

    path2 <- label2path(cut2, gtwInfo2)
    tempPath[[smo]] <- path2

    statResult <- smoTest(xcmsLargeWin, groupInd, parspec[1:parnum, ],
                          scanRange, parind[1:parnum], path2)
    tempRange <- statResult[3, 2]

    brokenNum1 <- sum((cut[eeBet[ ,1]] + cut[eeBet[ ,2]]) == 1)
    smoM[smo, 2] <- brokenNum1
    smoM[smo, 3] <- tempRange

    if (tempRange <= rangeThre)
      if (smo == 3)
        break

    if (smo == 1){
      smoM[2, 1] <- (costMax - cut[length(cut)]) / brokenNum1 + smoTemp
      tempCost <- cut[length(cut)] - smoTemp * brokenNum1
      smoM[3, 1] <- (tempCost - costMin) / (s1BrokenNum - brokenNum1)
    }
  }

  bestSmo <- which(smoM[ , 3] < rangeThre)
  if (length(bestSmo) > 1){
    bestSmo <- bestSmo[which.max(smoM[bestSmo, 1])]
  } else if (length(bestSmo) == 0){
    bestSmo <- which.min(smoM[ , 3])
  }

  # best_s = 2
  #[~, best_s] = min(smo_M(:,3));

  path2 <- tempPath[[bestSmo]]


  tempWarped <- parspec * 0
  for (m in 1:parnum)
    tempWarped[m, ] <- warpCurve(parspec[m, ], path2[[m]])

  parwarped <- tempWarped

  tempPath2 <- vector('list', dim(parspec)[1])
  tempPath2[1:parnum] <- path2

  parpath <- tempPath2
  parstat <- smoM

  return(list(parPath = parpath, parWarped = parwarped, parStat = parstat,
              parS1BrokenNum = parS1BrokenNum, parS2BrokenNum = parS2BrokenNum, parS2 = parS2))

}

