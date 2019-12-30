#' Run ncGTW alignment
#'
#' This function applies ncGTW alignment to the input feature.
#' @param ncGTWinput A \code{\link{ncGTWinput}} object.
#' @param xcmsLargeWin A \code{\link[xcms]{xcmsSet-class}} object.
#' @param parSamp Decide how many samples are in each group when considering
#'   parallel computing, and the default is 10.
#' @param k1Num Decide how many different k1 will be tested in stage 1. The
#'   default is 3.
#' @param k2Num Decide how many different k2 will be tested in stage 2. The
#'   default is 1.
#' @param bpParam A object of \pkg{BiocParallel} to control parallel processing,
#'   and can be created by
#'   \code{\link[BiocParallel:SerialParam-class]{SerialParam}},
#'   \code{\link[BiocParallel:MulticoreParam-class]{MulticoreParam}}, or
#'   \code{\link[BiocParallel:SnowParam-class]{SnowParam}}.
#' @param ncGTWparam A \code{\link{ncGTWparam}} object.
#' @details This function realign the input feature with ncGTW alignment
#' function with given m/z and RT range.
#' @return A \code{\link{ncGTWoutput}} object.
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
#' ncGTWparam <- new("ncGTWparam")
#'
#' # run ncGTW alignment
#' ncGTWoutputs <- vector('list', length(ncGTWinputs))
#' for (n in seq_along(ncGTWinputs))
#'     ncGTWoutputs[[n]] <- ncGTWalign(ncGTWinputs[[n]], xcmsLargeWin, 5,
#'         ncGTWparam = ncGTWparam)
#' }
#' @export

ncGTWalign <- function(ncGTWinput, xcmsLargeWin, parSamp=10, k1Num=3, k2Num=1,
    bpParam=BiocParallel::SnowParam(workers=1), ncGTWparam=NULL){

    if (!is(xcmsLargeWin, 'xcmsSet'))
        stop('xcmsLargeWin should be a "xcmsSet" object.')
    if (!is(ncGTWinput, 'ncGTWinput'))
        stop('ncGTWinput should be a "ncGTWinput" object.')
    if (parSamp <= 1 || parSamp %% 1 != 0)
        stop('parSamp should be an integer which is larger than 1.')
    if (k1Num <= 0 || k1Num %% 1 != 0)
        stop('k1Num should be an integer which is larger than 0.')
    if (k2Num <= 0 || k2Num %% 1 != 0)
        stop('k2Num should be an integer which is larger than 0.')

    ncGTWparam <- initHidparam(ncGTWparam, ncGTWinput, k1Num, k2Num)
    groupInd <- ncGTWinput@groupInfo[1]
    scanRangeOld <- ncGTWinput@rtRaw

    ncGTWpar <- ncGTWsplit(ncGTWinput, ncGTWparam, parSamp)
    parSamp <- ncGTWpar$parSamp
    parSect <- ncGTWpar$parSect
    parNum <- ncGTWpar$parNum
    specLen <- ncGTWpar$specLen
    scanRange <- ncGTWpar$scanRange

    message("ncGTW is realigning group ", groupInd, "...")

    ## ncGTW 1st layer alignment
    ncGTW1stOutput <- BiocParallel::bplapply(ncGTWpar$parInfos, ncGTW1stLayer,
        xcmsLargeWin, scanRange, scanRangeOld, ncGTWparam, BPPARAM=bpParam)

    parPath <- replicate(parSect, list(vector('list', parSamp)))
    parWarped <- array(0, c(parSamp, specLen, parSect))
    parStat <- array(0, c(k1Num * k2Num, 5, parSect))

    for (n in seq_len(parSect)){
        parPath[[n]] <- ncGTW1stOutput[[n]]$parPath
        parWarped[seq_len(parNum[n]) , , n] <- ncGTW1stOutput[[n]]$parWarped
        parStat[ , , n] <- ncGTW1stOutput[[n]]$parStat
    }

    ## Stop if only 1 layer
    if (parSect == 1){
        warning('Only one layer!!!')
        warpedAll <- parWarped[ , , 1]
        smoM2 <- matrix(0, 3, 6)

        return(new("ncGTWoutput", alignData=ncGTWpar$data, scanRange=scanRange,
            ncGTWpath=parPath[[1]], downSample=downSample))
    }

    ## Build super samples for 2nd layer
    superSample <- matrix(0, parSect, specLen)
    for (n in seq_len(parSect)){
        maxInd <- which.max(rowSums(parWarped[ , , n], 2))
        minInd <- which.min(rowSums(parWarped[ , , n], 2))
        superSample[n, ] <- colMeans(parWarped[setdiff(seq_len(parNum[n]),
            c(maxInd, minInd)), , n, drop=FALSE])
    }

    ## ncGTW 2nd layer alignment
    ncGTWparam2 <- ncGTWparam
    if (ncGTWparam2$maxStp * 1.5 > specLen){
        ncGTWparam2$maxStp <- specLen - 1
    } else{
        ncGTWparam2$maxStp <- round(ncGTWparam$maxStp * 1.5)
    }
    ncGTWparam2$mir <- FALSE

    allInfo <- list(groupInd=groupInd, num=0, all=0, parSpec=superSample,
        parNum=parSect, parInd=seq_len(parSect))

    ncGTWsmoRes2 <- ncGTWalignSmo(allInfo, xcmsLargeWin, scanRange,
        scanRangeOld, ncGTWparam2, TRUE, ncGTWpar$data, parPath,
        ncGTWpar$parInd)

    smoM2 <- ncGTWsmoRes2$smoM
    minRange <- min(smoM2[ , 3])
    bestSmo <- which(smoM2[ , 3] == minRange)
    bestSmo <- bestSmo[which.min(smoM2[bestSmo, 1])]

    path2 <- ncGTWsmoRes2$tempPath[[bestSmo]]

    ## Warp samples
    parSuperWarped <- matrix(0, parSect, specLen)
    for (m in seq_len(parSect))
        parSuperWarped[m, ] <- warpCurve(superSample[m, ], path2[[m]])

    warpedAll <- matrix(0, dim(ncGTWpar$data)[1], dim(ncGTWpar$data)[2])
    tempCount <- 0
    for (n in seq_len(parSect)){
        for (m in seq_len(parNum[n]))
            warpedAll[tempCount + m, ] <-
                warpCurve(parWarped[m, , n], path2[[n]])
        tempCount <- tempCount + parNum[n]
    }

    path <- pathCombine(parPath, path2, ncGTWpar$parInd)

    return(new("ncGTWoutput", alignData=ncGTWpar$data, scanRange=scanRange,
        ncGTWpath=path, downSample=ncGTWparam$downSample))
}


initHidparam <- function(ncGTWparam, ncGTWinput, k1Num, k2Num){

    ncGTWparam2 <- list()

    ncGTWparam2$downSample <- if (is.null(ncGTWparam@downSample)) 2 else
        ncGTWparam@downSample

    ncGTWparam2$nor <- if (is.null(ncGTWparam@nor)) 1 else ncGTWparam@nor

    ncGTWparam2$strNum <- if (is.null(ncGTWparam@strNum)) 1 else
        ncGTWparam@strNum

    ncGTWparam2$diaNum <- if (is.null(ncGTWparam@diaNum)) 1 else
        ncGTWparam@diaNum

    ncGTWparam2$mir <- TRUE

    ncGTWparam2$stpRat <- if (is.null(ncGTWparam@stpRat)) 0.6 else
        ncGTWparam@stpRat

    ncGTWparam2$maxStp <- if (is.nan(ncGTWparam@maxStp))
        round(dim(ncGTWinput@rtRaw)[2] %/% ncGTWparam2$downSample *
            ncGTWparam2$stpRat) else ncGTWparam@maxStp

    ncGTWparam2$rangeThre <- 1
    ncGTWparam2$biP <- TRUE
    ncGTWparam2$mu <- 0
    ncGTWparam2$sigma <- 1
    ncGTWparam2$weiP <- 0
    ncGTWparam2$logt <- 0
    ncGTWparam2$dia <- 0
    ncGTWparam2$noiseVar <- 1
    ncGTWparam2$k1Num <- k1Num
    ncGTWparam2$k2Num <- k2Num

    return(ncGTWparam2)
}


ncGTWsplit <- function(ncGTWinput, ncGTWparam, parSamp){

    downSample <- ncGTWparam$downSample

    dataOri <- ncGTWinput@profiles
    dataFil <- t(apply(dataOri, 1, gaussFilter))

    dataNum <- dim(dataFil)[1]
    specLen <- dim(dataFil)[2]

    data <- matrix(0, dataNum, specLen %/% downSample)
    scanRange <- matrix(0, dataNum, specLen %/% downSample)
    for (i in seq_len(dataNum)){
        for (j in seq_len(specLen %/% downSample)){
            subInd <- ((j - 1) * downSample + 1):(j * downSample)
            maxInd <- which.max(dataFil[i, subInd])
            data[i, j] <- dataFil[i, subInd[maxInd]]
            scanRange[i, j] <- ncGTWinput@rtRaw[i, subInd[maxInd]]
        }
    }
    data <- cbind(dataFil[ , 1], data)
    scanRange <- cbind(ncGTWinput@rtRaw[ , 1], scanRange)
    if (any(data[ , dim(data)[2]] != dataFil[ , specLen])){
        data <- cbind(data, dataFil[ , specLen])
        scanRange <- cbind(scanRange, ncGTWinput@rtRaw[ , specLen])
    }
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
    for (n in seq_len(parSect)){
        parSpec[seq_len(parNum[n]), , n] <-
            data[(tempCount + 1):(tempCount + parNum[n]), ]
        parInd[seq_len(parNum[n]), n] <- (tempCount + 1):(tempCount + parNum[n])
        tempCount <- tempCount + parNum[n]
    }

    parInfos <- vector('list', parSect)
    for (n in seq_len(parSect)){
        parInfos[[n]]$groupInd <- ncGTWinput@groupInfo['index']
        parInfos[[n]]$num <- n
        parInfos[[n]]$all <- parSect
        parInfos[[n]]$parSpec <- parSpec[seq_len(parNum[n]), , n]
        parInfos[[n]]$parNum <- parNum[n]
        parInfos[[n]]$parInd <- parInd[ , n]
    }
    return(list(parInfos=parInfos, scanRange=scanRange, parSamp=parSamp,
        specLen=specLen, parSect=parSect, parNum=parNum, data=data,
        parInd=parInd))
}


ncGTW1stLayer <- function(parInfo, xcmsLargeWin, scanRange, scanRangeOld,
    ncGTWparam){

    ncGTWsmoRes <- ncGTWalignSmo(parInfo, xcmsLargeWin, scanRange, scanRangeOld,
        ncGTWparam)

    warpRes <- chooseSmo(ncGTWsmoRes$smoM, ncGTWsmoRes$tempPath, parInfo,
        ncGTWparam)

    return(list(parPath=warpRes$parPath, parWarped=warpRes$parWarped,
        parStat=warpRes$parStat, parS1BrokenNum=ncGTWsmoRes$parS1BrokenNum,
        parS2BrokenNum=ncGTWsmoRes$parS2BrokenNum, parS2=ncGTWsmoRes$parS2))
}

ncGTWalignSmo <- function(parInfo, xcmsLargeWin, scanRange,
    scanRangeOld, ncGTWparam, second=FALSE, data=NULL, parPath=NULL,
    parInd=NULL){

    groupInd <- parInfo$groupInd
    n <- parInfo$num
    parSect <- parInfo$all
    parnum <- parInfo$parNum
    parind <- parInfo$parInd[seq_len(parnum)]
    parspec <- parInfo$parSpec[seq_len(parnum), ]

    mir <- ncGTWparam$mir
    strNum <- ncGTWparam$strNum
    diaNum <- ncGTWparam$diaNum
    biP <- ncGTWparam$biP
    noiseVar <- ncGTWparam$noiseVar
    downSample <- ncGTWparam$downSample
    rangeThre <- ncGTWparam$rangeThre
    maxStp <- ncGTWparam$maxStp
    k1Num <- ncGTWparam$k1Num
    k2Num <- ncGTWparam$k2Num

    gtwPrep <- buildMultiParaValidmap(parspec, mir, strNum, diaNum, biP)
    eeBet <- betEdge(gtwPrep, 10^-32, ncGTWparam)

    gtwRes0 <- gtwCut(gtwPrep, 0, ncGTWparam)
    cut0 <- gtwRes0$cut
    costMin <- cut0[length(gtwRes0$cut)]
    s1BrokenNum <- sum((cut0[eeBet[ ,1]] + cut0[eeBet[ ,2]]) == 1)
    parS1BrokenNum <- s1BrokenNum #####

    gtwResMerge <- gtwCut(gtwPrep, 10^32, ncGTWparam, 2)
    cutMerge <- gtwResMerge$cut
    costMax <- cutMerge[length(cutMerge)]
    smoTemp <- (costMax - costMin) / s1BrokenNum

    BNumDTW1 <- matrix(-1, 3, 3)
    BNumDTW1[1, 1] <- 0
    BNumDTW1[2, 1] <- s1BrokenNum
    BNumDTW1[3, 1] <- costMin
    BNumDTW1[1, 3] <- Inf
    BNumDTW1[2, 3] <- 0
    BNumDTW1[3, 3] <- costMax

    smoM <- array(999, c(k1Num, k2Num, 5))
    smoM[1, ,1] <- smoTemp
    BNumDTW1[1, 2] <- smoTemp
    tempPath <- vector('list', k1Num * k2Num)

    for (smo in seq_len(dim(smoM)[1])){
        gtwRes1 <- gtwCut(gtwPrep, smoM[smo, 1,1], ncGTWparam)
        cut1 <- gtwRes1$cut
        gtwInfo1 <- gtwRes1$gtwInfo
        path1 <- label2path(gtwRes1$cut, gtwRes1$gtwInfo)
        if (n != 0){
            message(format(Sys.time()),
                ' Solved the 1st maximum flow, section ', n, ' of ', parSect,
                ', bottom layer.')
        } else{
            message(format(Sys.time()),
                ' Solved the 1st maximum flow, top layer.')
        }
        brokenNum1 <- sum((cut1[eeBet[ ,1]] + cut1[eeBet[ ,2]]) == 1)
        tempInd <- (smo - 2^floor(log2(smo)) + 1) * 2
        BNumDTW1[2, tempInd] <- brokenNum1
        BNumDTW1[3, tempInd] <- cut1[length(cut1)] - smoM[smo, 1,1] * brokenNum1
        smoM[smo, ,2] <- brokenNum1

        validMap2 <- matrix(0, parnum, parnum)
        for (jj in seq_len(parnum))
            for (ii in seq((jj + 1), parnum, length = max(0, parnum - jj)))
                validMap2[ii, jj] <- 1
        tst2 <- matrix(0, dim(parspec)[1], dim(parspec)[2])
        ref2 <- tst2
        smoothnessV2 <- 0.0000000001 / noiseVar

        gtwPrep2 <- list(validMap=validMap2, tst=tst2, ref=ref2)
        gtwRes2 <- gtwCut(gtwPrep2, smoothnessV2, ncGTWparam, 1, path1,
            gtwPrep$pairMap)
        eeBet2 <- gtwRes2$eeBet
        s2BrokenNum <- gtwRes2$cut[length(gtwRes2$cut)] / smoothnessV2
        parS2BrokenNum <- s2BrokenNum #####
        parS2 <- 2 * (maxStp - 1) * parnum / s2BrokenNum ###############
        smoothnessV2 <- parS2 / noiseVar

        BNumDTW2 <- matrix(-1, 3, 3)
        BNumDTW2[1, 1] <- 0
        BNumDTW2[2, 1] <- s2BrokenNum
        BNumDTW2[3, 1] <- 0
        BNumDTW2[1, 3] <- Inf
        BNumDTW2[2, 3] <- 0
        BNumDTW2[3, 3] <- 2 * (maxStp - 1) * parnum
        smoM[smo, 1, 4] <- smoothnessV2
        BNumDTW2[1, 2] <- smoothnessV2

        for (smo2 in seq_len(dim(smoM)[2])){

            gtwRes2 <- gtwCut(gtwPrep2, smoM[smo, smo2,4], ncGTWparam, 1, path1,
                gtwPrep$pairMap)
            cut2 <- gtwRes2$cut
            gtwInfo2 <- gtwRes2$gtwInfo
            path2 <- label2path(cut2, gtwInfo2)
            tempPath[[(smo-1) * k2Num + smo2]] <- path2
            if (n != 0){
                message(format(Sys.time()),
                    ' Solved the 2nd maximum flow, section ', n, ' of ', parSect,
                    ', bottom layer.')
            } else{
            message(format(Sys.time()),
                ' Solved the 2nd maximum flow, top layer.')
            }

            brokenNum2 <- sum((cut2[eeBet2[ ,1]] + cut2[eeBet2[ ,2]]) == 1)
            tempInd <- (smo2 - 2^floor(log2(smo2)) + 1) * 2
            BNumDTW2[2, tempInd] <- brokenNum2
            BNumDTW2[3, tempInd] <- cut2[length(cut2)] -
                smoM[smo, smo2, 4] * brokenNum2
            smoM[smo, smo2, 5] <- brokenNum2

            if (second){
                temp2path <- pathCombine(parPath, path2, parInd)
                statResult <- smoTest(xcmsLargeWin, groupInd, data, scanRange,
                    seq_len(dim(data)[1]), temp2path, downSample, scanRangeOld)
            } else{
                statResult <- smoTest(xcmsLargeWin, groupInd, parspec, scanRange,
                    parind, path2, downSample, scanRangeOld)
            }

            tempRange <- statResult[3, 2]


            smoM[smo, smo2, 3] <- tempRange

            if (tempRange <= rangeThre)
                break

            if ((as.integer(log2(smo2+1)) - log2(smo2+1)) == 0 &&
                smo2 < dim(smoM)[2]){

                tBNumDTW2 <- matrix(-1, 3, smo2*2 + 3)
                for (tsmo2 in seq_len(dim(tBNumDTW2)[2])){
                    if (as.integer(tsmo2/2) != tsmo2/2){
                        tBNumDTW2[, tsmo2] = BNumDTW2[,(tsmo2+1)/2]
                    } else{
                        tBNumDTW2[1, tsmo2] =
                            (BNumDTW2[3, tsmo2/2 + 1] - BNumDTW2[3, tsmo2/2]) /
                            (BNumDTW2[2, tsmo2/2] - BNumDTW2[2, tsmo2/2 + 1])
                        if (smo2 + tsmo2/2 <= dim(smoM)[2])
                            smoM[smo, smo2 + tsmo2/2, 4] = tBNumDTW2[1, tsmo2]
                    }
                }
                BNumDTW2 <- tBNumDTW2
            }
        }
        if (tempRange <= rangeThre)
            break

        if ((as.integer(log2(smo+1)) - log2(smo+1)) == 0 && smo < dim(smoM)[1]){
            tBNumDTW1 <- matrix(-1, 3, smo*2 + 3)
            for (tsmo in seq_len(dim(tBNumDTW1)[2])){
                if (as.integer(tsmo/2) != tsmo/2){
                    tBNumDTW1[, tsmo] = BNumDTW1[,(tsmo+1)/2]
                } else{
                    tBNumDTW1[1, tsmo] =
                        (BNumDTW1[3, tsmo/2 + 1] - BNumDTW1[3, tsmo/2]) /
                        (BNumDTW1[2, tsmo/2] - BNumDTW1[2, tsmo/2 + 1])
                    if (smo + tsmo/2 <= dim(smoM)[1])
                        smoM[smo + tsmo/2, ,1] = tBNumDTW1[1, tsmo]
                }
            }
            BNumDTW1 <- tBNumDTW1
        }
    }
    return(list(smoM=matrix(aperm(smoM, c(2,1,3)), k1Num * k2Num, 5),
        tempPath=tempPath, parS1BrokenNum = parS1BrokenNum,
        parS2BrokenNum = parS2BrokenNum, parS2 = parS2))
}


chooseSmo <- function(smoM, tempPath, parInfo, ncGTWparam, bestSmo=NULL){

    rangeThre <- ncGTWparam$rangeThre

    parnum <- parInfo$parNum
    parspec <- parInfo$parSpec[seq_len(parnum), ]

    if (is.null(bestSmo)){
        bestSmo <- which(smoM[ , 3] < rangeThre)##
        if (length(bestSmo) > 1){
            bestSmo <- bestSmo[which.max(smoM[bestSmo, 1] / smoM[bestSmo, 4])]
        } else if (length(bestSmo) == 0){
            bestSmo <- which.min(smoM[ , 3])
        }
    }

    path2 <- tempPath[[bestSmo]]##

    tempWarped <- parspec * 0
    for (m in seq_len(parnum))
        tempWarped[m, ] <- warpCurve(parspec[m, ], path2[[m]])

    parwarped <- tempWarped##

    tempPath2 <- vector('list', parnum)
    tempPath2[seq_len(parnum)] <- path2

    parpath <- tempPath2##
    parstat <- smoM    ##

    return(list(parPath=parpath, parWarped=parwarped, parStat=parstat))
}

gtwCut <- function(gtwPrep, smooth, ncGTWparam, type=0, path=NA, pairMap=NA){

    param <- initGtwParam(gtwPrep$validMap, ncGTWparam$noiseVar,
        ncGTWparam$maxStp, smooth, ncGTWparam$dia, ncGTWparam$logt,
        ncGTWparam$nor)

    gtwInfo <- buildGTWgraph(gtwPrep$ref, gtwPrep$tst, gtwPrep$validMap, param,
        ncGTWparam$mu, ncGTWparam$sigma, ncGTWparam$biP, ncGTWparam$weiP, type,
        path, pairMap)

    eeBet <- gtwInfo$ee[gtwInfo$ee[, 3] == smooth & gtwInfo$ee[, 4] == smooth, ]

    cut <- graphCut(gtwInfo$ss, gtwInfo$ee)

    return(list(gtwInfo = gtwInfo, eeBet = eeBet, cut = cut))
}



betEdge <- function(gtwPrep, smooth, ncGTWparam){

    param <- initGtwParam(gtwPrep$validMap, ncGTWparam$noiseVar,
        ncGTWparam$maxStp, smooth, ncGTWparam$dia, ncGTWparam$logt,
        ncGTWparam$nor)

    gtwInfo <- buildGTWgraph(gtwPrep$ref, gtwPrep$tst, gtwPrep$validMap, param,
        ncGTWparam$mu, ncGTWparam$sigma, ncGTWparam$biP, ncGTWparam$weiP)

    eeBet <-
        gtwInfo$ee[(gtwInfo$nEdgeGrid * gtwInfo$nPix + 1):dim(gtwInfo$ee)[1], ]

    return(eeBet)
}

#' @useDynLib ncGTW
#' @importFrom Rcpp sourceCpp
NULL
