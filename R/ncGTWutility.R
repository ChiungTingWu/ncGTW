#' Compute average pairwise correlation and overlapping area
#'
#' This function computes average pairwise correlation and overlapping area of
#' each sample pair.
#' @param ncGTWinput A list in which each element is a \code{\link{ncGTWinput}}
#'   object.
#' @param sampleRt A list of the same length as the sample number in which each
#'   element is a vector corresponding to the sample raw/adjusted RT.
#' @details This function computes the pairwise correlation and overlapping area
#' of each sample pair from the input feature, and then takes average.
#' @return A list in which the first element is average pairwise correlation,
#' and the second one is average overlapping area.
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
#'
#' # load the sample profiles
#' ncGTWinputs <- loadProfile(file, excluGroups)
#'
#' XCMSCor <- matrix(0, length(ncGTWinputs), 1)
#' XCMSOl <- matrix(0, length(ncGTWinputs), 1)
#' for (n in seq_along(ncGTWinputs)){
#'     XCMSmean <- meanCorOl(ncGTWinputs[[n]],
#'         slot(xcmsLargeWin, 'rt')$corrected)
#'     XCMSCor[n] <- XCMSmean$cor
#'     XCMSOl[n] <- XCMSmean$ol
#' }
#' @export

meanCorOl <- function(ncGTWinput, sampleRt){
    if (!is(ncGTWinput, 'ncGTWinput'))
        stop('ncGTWoutput should be a "ncGTWoutput" object.')

    samNum <- dim(ncGTWinput@rtRaw)[1]

    if (length(sampleRt) != samNum)
        stop('sampleRt should be a list with length as same as sample number.')

    pointNum <- dim(ncGTWinput@rtRaw)[2]
    profiles <- ncGTWinput@profiles
    rtRange <- matrix(0, samNum, pointNum)
    for (n in seq_len(samNum)){
        profiles[n, ] <- gaussFilter(profiles[n, ])
        rtRange[n, ] <- sampleRt[[n]][ncGTWinput@rtRaw[n, ]]
    }
    proInter <- matrix(0, samNum, pointNum * 10)
    interX <- seq(max(rtRange[ , 1]), min(rtRange[ , pointNum]),
                    length.out=pointNum * 10)
    for (n in seq_len(samNum))
        proInter[n, ] <- approx(rtRange[n, ], profiles[n, ], interX,
                                yleft=NA, yright=NA)$y
    corM <- cor(t(proInter))
    olM <- matrix(0, samNum, samNum)
    for (i in seq_len(samNum))
        for (j in i:samNum){
            olM[i, j] <- sum(pmin(proInter[i, ], proInter[j, ])) /
                min(sum(proInter[i, ]), sum(proInter[j, ]))
            olM[j, i] <- olM[i, j]
        }
    return(list(cor=mean(corM), ol=mean(olM)))
}


#' Compare CV
#'
#' This function calculates the coefficient of variation of each feature.
#' @param XCMSresFilled A \code{\link[xcms]{xcmsSet-class}} object.
#' @param na.rm Omit the samples in which the feature is not detected, and the
#'   default is FALSE.
#' @details This function calculates the coefficient of variation of each
#' feature across all the samples. If a sample is detected with more than one
#' peaks in the feature, the function will pick the one with the highest
#' intensity value.
#' @return A vector of the same length as the row number of the \code{group}
#' slot in \code{XCMSresFilled}, in which each element is the CV.
#' @examples
#' # obtain data
#' data('xcmsExamples')
#' xcmsLargeWin <- xcmsExamples$xcmsLargeWin
#'
#' cv <- compCV(xcmsLargeWin)
#' @export

compCV <- function(XCMSresFilled, na.rm = FALSE){
    if (!is(XCMSresFilled, 'xcmsSet'))
        stop('XCMSresFilled should be a "xcmsSet" object.')

    groupNum <- dim(XCMSresFilled@groups)[1]
    sampleNum <- max(XCMSresFilled@peaks[, 'sample'])

    XCMSpeaks <- matrix(0, groupNum, sampleNum)
    XCMScv <- matrix(0, groupNum, 1)

    for (n in seq_len(groupNum)){
        XCMSgroupPeaks <- XCMSresFilled@peaks[XCMSresFilled@groupidx[[n]], ]

        XCMSonePeak <- matrix(0, sampleNum, dim(XCMSresFilled@peaks)[2])
        colnames(XCMSonePeak) <- colnames(XCMSgroupPeaks)

        for (m in seq_len(sampleNum)){
            XCMSpeakInd <- which(XCMSgroupPeaks[, 'sample'] == m)
            XCMSpeakInd <- XCMSpeakInd[which.max(XCMSgroupPeaks[XCMSpeakInd,
                                                                'into'])]
            if (length(XCMSpeakInd) == 0){
                XCMSonePeak[m, ] <- NA
            } else{
                XCMSonePeak[m, ] <- XCMSgroupPeaks[XCMSpeakInd, ]
            }
        }
        XCMSpeaks[n, ] <- XCMSonePeak[ , 'into']
        XCMScv[n] <- sd(XCMSpeaks[n, ], na.rm = na.rm)/mean(XCMSpeaks[n, ],
                                                            na.rm = na.rm)
    }
    return(XCMScv)
}




gaussFilter <- function(prof, sig=1){
    sz <- ceiling(sig * 6)    # length of gaussian filter vector
    if (sz < 2){
        sz <- 2
    }
    if (sz %% 2 != 0){
        sz <- sz + 1
    }
    x <- seq(-sz / 2, sz / 2, length=sz + 1)
    filterVec <- exp(-x ^ 2 / (2 * sig ^ 2))
    filterVec <- filterVec / sum (filterVec) # normalize
    filtered <- convolve(prof, filterVec, type='open')
    return(filtered[(sz / 2 + 1):(length(filtered) - sz / 2)])
}

rt2scan <- function(rt, rtAll)
    return(which.min(abs(rtAll - rt)))

smoTest <- function(xcmsLargeWin, groupInd, dataSub, scanRange,
                    sampleInd, path2, downSample, scanRangeOld){

    peaks <- xcmsLargeWin@peaks
    groupidx <- xcmsLargeWin@groupidx
    rtXCMS <- xcmsLargeWin@rt$corrected
    rtRaw <- xcmsLargeWin@rt$raw

    prePeaks <- findUniPeak(peaks, groupInd, groupidx, sampleInd=sampleInd)

    if (length(prePeaks) == 0 || length(unique(prePeaks[ , 'sample'])) == 1)
        return(matrix(-1, 3, 3))

    prePeakInd <- prePeaks[ , 'sample']
    prePeakMed <- prePeaks[, 'rt']

    sampleCount <- table(prePeakInd)
    groupNum <- max(sampleCount)
    groupSam <- as.numeric(names(which(sampleCount == groupNum)))
    maxNum <- 0
    maxRange <- 0
    maxInd <- 1
    for (ind in seq_len(length(groupSam))){
        samPeaks <- prePeaks[prePeaks[ , 'sample'] == groupSam[ind], 'rt']
        tempNum <- length(samPeaks)
        if (tempNum < maxNum )
            next
        tempRange <- max(samPeaks) - min(samPeaks)
        if (tempNum > maxNum ){
            maxNum <- tempNum
            maxRange <- tempRange
            maxInd <- groupSam[ind]
        } else if (tempRange > maxRange){
            maxRange <- tempRange
            maxInd <- groupSam[ind]
        }
    }
    groupSam <- maxInd

    if (groupNum == 1){
        kmeansPreInd <- kmeans(prePeakMed, 1)
    } else{
        kmeansPreInd <- kmeans(prePeakMed, prePeakMed[prePeakInd == groupSam])
    }

    oriPeakGroup <- vector('list', groupNum)
    XCMSPeakGroup <- vector('list', groupNum)
    ncGTWPeakGroup <- vector('list', groupNum)

    for (n in seq_len(groupNum))
        XCMSPeakGroup[[n]] <- prePeakMed[kmeansPreInd$cluster == n]

    ncGTWPeakMed <- prePeakMed * 0
    oriPeakMed <- prePeakMed * 0

    for (n in seq_len(length(ncGTWPeakMed))){
        samInd <- prePeaks[n, 'sample']
        samSubInd <- which(sampleInd == prePeaks[n, 'sample'])

        indDif <- abs(scanRange[samInd, ] - rt2scan(prePeakMed[n],
                                                    rtXCMS[[samInd]]))
        minIndDif <- min(indDif)
        medInd <- which(indDif == minIndDif)
        medInd <- medInd[which.max(dataSub[samSubInd, medInd])]
        fRange <- round(3 / mean(diff(rtRaw[[samInd]][scanRange[samInd,]])))
        if (medInd - fRange < 1){
            staInd <- 1
        } else {
            staInd <- medInd - fRange
        }
        if (medInd + fRange > dim(dataSub)[2]){
            endInd <- dim(dataSub)[2]
        } else{
            endInd <- medInd + fRange
        }
        apexRange <- staInd:endInd
        apexInd <- apexRange[which.max(dataSub[samSubInd, apexRange])]
        oriPeakMed[n] <- rtRaw[[samInd]][scanRange[samInd, apexInd]]

        samPath <- path2[[samSubInd]]


        scanSubncGTW <-
            round(mean(samPath[which(samPath[ , 2] == apexInd), 1]))
        scanncGTW <- (scanSubncGTW - 1) * downSample
        if (scanncGTW > dim(scanRangeOld)[2])
            scanncGTW <- dim(scanRangeOld)[2]
        ncGTWPeakMed[n] <- scanSubncGTW
    }
    oriPeakRt <- cbind(oriPeakMed, prePeaks[, 'rtmin'] - prePeakMed +
        oriPeakMed, prePeaks[, 'rtmax'] - prePeakMed + oriPeakMed)

    ncGTWPeakRt <- cbind(ncGTWPeakMed, ncGTWPeakMed)

    kmeansOriInd <-
        kmeans(oriPeakRt, oriPeakRt[prePeakInd == groupSam, , drop=FALSE])
    kmeansncGTWInd <-
        kmeans(ncGTWPeakRt, unique(ncGTWPeakRt[prePeakInd == groupSam, ,
            drop=FALSE]))

    for (n in seq_len(groupNum)){
        oriPeakGroup[[n]] <- oriPeakMed[kmeansOriInd$cluster == n]
        ncGTWPeakGroup[[n]] <- ncGTWPeakMed[kmeansncGTWInd$cluster == n]
    }

    statResult <- matrix(0, 3, 2)
    statResult[1, 1] <-
        sum(vapply(oriPeakGroup, var, vector("double", 1)), na.rm=TRUE)
    statResult[2, 1] <-
        sum(vapply(XCMSPeakGroup, var, vector("double", 1)), na.rm=TRUE)
    statResult[3, 1] <-
        sum(vapply(ncGTWPeakGroup, var, vector("double", 1)), na.rm=TRUE)
    statResult[1, 2] <- max(vapply(oriPeakGroup,
        function(x) range(x)[2] - range(x)[1], vector("double", 1)))
    statResult[2, 2] <- max(vapply(XCMSPeakGroup,
        function(x) range(x)[2] - range(x)[1], vector("double", 1)))
    statResult[3, 2] <- max(vapply(ncGTWPeakGroup,
        function(x) range(x)[2] - range(x)[1], vector("double", 1)))

    return(statResult)
}
