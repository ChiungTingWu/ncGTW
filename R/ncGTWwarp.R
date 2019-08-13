#' Adjust retention time
#'
#' This function produces the new warping functions (RT lists) with the
#' realignment result.
#' @param xcmsLargeWin A \code{\link[xcms]{xcmsSet-class}} object.
#' @param ncGTWinput A \code{\link{ncGTWinput}} object.
#' @param ncGTWoutput A \code{\link{ncGTWoutput}} object.
#' @param ppm Should be set as same as the one when performing the peak
#'   detection function in \code{xcms}.
#' @details This function produces the new warping functions (RT lists) with the
#' realignment result.
#' @return A \code{\link{ncGTWwarp}} object.
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
#'
#' # adjust RT with the realignment results from ncGTW
#' ncGTWres <- xcmsLargeWin
#' ncGTWRt <- vector('list', length(ncGTWinputs))
#' for (n in seq_along(ncGTWinputs)){
#'     adjustRes <- adjustRT(ncGTWres, ncGTWinputs[[n]], ncGTWoutputs[[n]], ppm)
#'     xcms::peaks(ncGTWres) <- ncGTWpeaks(adjustRes)
#'     ncGTWRt[[n]] <- rtncGTW(adjustRes)
#' }
#'
#' # apply the adjusted RT to a xcmsSet object
#' xcms::groups(ncGTWres) <- excluGroups[ , 2:9]
#' xcms::groupidx(ncGTWres) <- xcms::groupidx(xcmsLargeWin)[excluGroups[ , 1]]
#' rtCor <- vector('list', length(xcms::filepaths(ncGTWres)))
#' for (n in seq_along(file)){
#'     rtCor[[n]] <- vector('list', dim(excluGroups)[1])
#'     for (m in seq_len(dim(excluGroups)[1]))
#'         rtCor[[n]][[m]] <- ncGTWRt[[m]][[n]]
#' }
#' slot(ncGTWres, 'rt')$corrected <- rtCor
#' }
#' @export

adjustRT <- function(xcmsLargeWin, ncGTWinput, ncGTWoutput, ppm){
    if (!is(xcmsLargeWin, 'xcmsSet'))
        stop('xcmsLargeWin should be a "xcmsSet" object.')
    if (!is(ncGTWinput, 'ncGTWinput'))
        stop('ncGTWinput should be a "ncGTWinput" object.')
    if (!is(ncGTWoutput, 'ncGTWoutput'))
        stop('ncGTWoutput should be a "ncGTWoutput" object.')
    if (ppm <= 0)
        stop('ppm should be larger than 0.')

    peaks <- xcmsLargeWin@peaks
    groupInfo <- ncGTWinput@groupInfo
    groupInd <- groupInfo['index']
    rtXCMS <- xcmsLargeWin@rt$corrected
    rtRaw <- xcmsLargeWin@rt$raw
    scanRangeOld <- ncGTWinput@rtRaw
    dataOri <- ncGTWinput@profiles
    data <- ncGTWoutput@alignData
    path <- ncGTWoutput@ncGTWpath
    scanRange <- ncGTWoutput@scanRange
    downSample <- ncGTWoutput@downSample

    ## Find unique peaks
    XCMSPeaks <- findUniPeak(peaks, groupInd, xcmsLargeWin@groupidx, ppm)

    ## Find the maximum number of peaks within a sample
    sampleCount <- table(XCMSPeaks[ , 'sample'])
    groupNum <- max(sampleCount)

    ## Find corresponding scan and RT between XCMS and ncGTW
    ncGTWpeaks <- XCMSPeaks
    ncGTWpeaks[ , c('rt', 'rtmin', 'rtmax')] <- 0
    peakScanRt <- matrix(0, dim(ncGTWpeaks)[1], 11)
    colnames(peakScanRt) <-
        c('scan_XCMS', 'scan_ncGTW', 'subscan_ncGTW', 'rt_ncGTW',
            'scanmin_XCMS', 'scanmin_ncGTW', 'rtmin_ncGTW',
            'scanmax_XCMS', 'scanmax_ncGTW', 'rtmax_ncGTW', 'sample')

    warpOrder <-
        sort(XCMSPeaks[ , 'maxo'], index.return = TRUE, decreasing = TRUE)$ix
    for (nn in seq_len(dim(ncGTWpeaks)[1])){
        n <- warpOrder[nn]
        samInd <- XCMSPeaks[n, 'sample']
        samPath <- path[[samInd]]
        samRtXCMS <- rtXCMS[[samInd]][scanRange[samInd, ]]
        samRtXCMSOld <- rtXCMS[[samInd]][scanRangeOld[samInd, ]]
        samRtRaw <- rtRaw[[samInd]][scanRange[samInd, ]]
        samRtRawOld <- rtRaw[[samInd]][scanRangeOld[samInd, ]]
        meanDif <- mean(diff(samRtXCMSOld))

        ### Find the true apex points
        fRange <- round(1 / meanDif)
        scanXCMS <- rt2scan(XCMSPeaks[n, 'rt'], samRtXCMSOld)
        scanXCMSrange <- max(1, (scanXCMS - fRange)) :
            min(length(samRtXCMSOld), (scanXCMS + fRange))

        rE <- scanXCMSrange[length(scanXCMSrange)]
        for (rangeInd in seq_len(fRange)){
            ind <- scanXCMS + rangeInd
            if (ind > length(samRtXCMSOld))
                break
            if (samRtXCMSOld[ind] - samRtXCMSOld[ind - 1] > 2 * meanDif){
                rE <- ind - 1
                break
            }
        }
        rS <- scanXCMSrange[1]
        for (rangeInd in seq_len(fRange)){
            ind <- scanXCMS - rangeInd
            if (ind < 1)
                break
            if (samRtXCMSOld[ind + 1] - samRtXCMSOld[ind] > 2 * meanDif){
                rS <- ind + 1
                break
            }
        }
        scanXCMSrange <- rS:rE

        scanXCMS <- scanXCMSrange[which.max(dataOri[samInd, scanXCMSrange])]
        XCMSrt <- samRtXCMSOld[scanXCMS]
        scanSubXCMS <- rt2scan(XCMSrt, samRtXCMS)

        ### Match the apex point
        scanSubncGTW <-
            round(mean(samPath[which(samPath[ ,2] == scanSubXCMS), 1]))
        scanncGTW <- (scanSubncGTW - 1) * downSample
        if (scanncGTW > dim(scanRangeOld)[2])
            scanncGTW <- dim(scanRangeOld)[2]
        if (scanncGTW < 1)
            scanncGTW <- 1
        ncGTWrt <- rtRaw[[samInd]][scanRangeOld[samInd, scanncGTW]]

        peakScanRt[n, 'rt_ncGTW'] <- ncGTWrt
        peakScanRt[n, 'subscan_ncGTW'] <- scanSubncGTW
        peakScanRt[n, 'scan_XCMS'] <-
            which(scanRangeOld[samInd, ] == scanRange[samInd, scanSubXCMS])
        peakScanRt[n, 'scan_ncGTW'] <- scanncGTW

        ### Avoid overlapping apex points
        samPeaks <-
            cbind(peakScanRt[peakScanRt[ , 'sample'] == samInd, , drop=FALSE],
                peakInd=which(peakScanRt[ , 'sample'] == samInd))

        for (m in seq_len(dim(samPeaks)[1])){
            tempInd <- samPeaks[m, 'peakInd']
            apexDif <- XCMSPeaks[n, 'rt'] - XCMSPeaks[tempInd, 'rt']
            tempSign <- apexDif / abs(apexDif)

            if (apexDif > 0){
                widthDif <- XCMSPeaks[n, 'rt'] - XCMSPeaks[n, 'rtmin'] +
                    XCMSPeaks[tempInd, 'rtmax'] - XCMSPeaks[tempInd, 'rt']
            } else{
                widthDif <- XCMSPeaks[n, 'rtmax'] - XCMSPeaks[n, 'rt'] +
                    XCMSPeaks[tempInd, 'rt'] - XCMSPeaks[tempInd, 'rtmin']
            }
            peakRtDif <- min(abs(apexDif), widthDif)

            if (tempSign * (peakScanRt[n, 'rt_ncGTW'] - samPeaks[m, 'rt_ncGTW'])
                < peakRtDif){
                peakScanRt[n, 'scan_ncGTW'] <-
                    rt2scan(samPeaks[m, 'rt_ncGTW'] +
                                tempSign * peakRtDif, samRtXCMSOld)
                peakScanRt[n, 'rt_ncGTW'] <-
                    samRtRawOld[peakScanRt[n, 'scan_ncGTW']]
            }
        }

        ### Match the start and end points
        peakScanRt[n, 'rtmin_ncGTW'] <-
            XCMSPeaks[n, 'rtmin'] - XCMSPeaks[n, 'rt'] +
            peakScanRt[n, 'rt_ncGTW']
        peakScanRt[n, 'rtmax_ncGTW'] <-
            XCMSPeaks[n, 'rtmax'] - XCMSPeaks[n, 'rt'] +
            peakScanRt[n, 'rt_ncGTW']
        peakScanRt[n, 'scanmin_ncGTW'] <-
            rt2scan(peakScanRt[n, 'rtmin_ncGTW'], samRtRawOld)
        peakScanRt[n, 'scanmax_ncGTW'] <-
            rt2scan(peakScanRt[n, 'rtmax_ncGTW'], samRtRawOld)
        peakScanRt[n, 'scanmin_XCMS'] <- peakScanRt[n, 'scanmin_ncGTW'] -
            peakScanRt[n, 'scan_ncGTW'] + peakScanRt[n, 'scan_XCMS']
        peakScanRt[n, 'scanmax_XCMS'] <- peakScanRt[n, 'scanmax_ncGTW'] -
            peakScanRt[n, 'scan_ncGTW'] + peakScanRt[n, 'scan_XCMS']
        peakScanRt[n, 'sample'] <- XCMSPeaks[n, 'sample']

        ncGTWpeaks[n, c('rt', 'rtmin', 'rtmax')] <-
            peakScanRt[n, c('rt_ncGTW', 'rtmin_ncGTW', 'rtmax_ncGTW')]
        peaks[ncGTWpeaks[n, 'peakInd'], c('rt', 'rtmin', 'rtmax')] <-
            ncGTWpeaks[n, c('rt', 'rtmin', 'rtmax')]

    }

    ## Find a sample as the center for kmeans
    groupSam <- as.numeric(names(which(sampleCount == groupNum)))
    maxNum <- 0
    maxRange <- 0
    maxInd <- 1
    for (ind in seq_len(length(groupSam))){
        samPeaks <- XCMSPeaks[XCMSPeaks[ , 'sample'] == groupSam[ind], 'rt']
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


    ## kmeans grouping
    cenInd <- which(peakScanRt[ , 'sample'] == groupSam)
    kcenter <- peakScanRt[cenInd, 'rt_ncGTW', drop=FALSE]

    intOrder <- sort(XCMSPeaks[cenInd, 'maxo'],
                    index.return=TRUE, decreasing=TRUE)$ix[seq_len(groupNum)]
    kcenter <- kcenter[intOrder, , drop=FALSE]
    kcenter <- unique(kcenter)
    groupNum <- length(kcenter)


    if (length(kcenter) == 1){
        kmeansncGTW <- kmeans(peakScanRt[ , 'rt_ncGTW', drop=FALSE], 1)
    } else{
        kmeansncGTW <- kmeans(peakScanRt[ , 'rt_ncGTW', drop=FALSE], kcenter)
    }

    ## Put peaks into their groups
    peakGroup <- vector('list', groupNum)
    tempPeaks <- matrix(NA, max(peaks[ , 'sample']), dim(peakScanRt)[2] + 2)
    colnames(tempPeaks) <- c(colnames(peakScanRt), 'detect', 'maxo')
    colnames(tempPeaks)[c(1, 5, 8)] <-
        c('scan_Raw', 'scanmin_Raw', 'scanmax_Raw')
    for (n in seq_len(groupNum))
        peakGroup[[n]] <- tempPeaks

    for (n in seq_len(dim(peakScanRt)[1])){
        gInd <- kmeansncGTW$cluster[n]
        sInd <- peakScanRt[n, 'sample']

        if (is.na(peakGroup[[gInd]][sInd, 1]) ||
            (XCMSPeaks[n, 'maxo'] > peakGroup[[gInd]][sInd, 'maxo'])){

            peakGroup[[gInd]][sInd, ] <-
                c(peakScanRt[n, ], 1, XCMSPeaks[n, 'maxo'])
        }
    }

    ## Remove groups with too few samples
    rmInd <- matrix(TRUE, length(peakGroup), 1)
    for (n in seq_len(groupNum)){
        tempPeaks <- peakGroup[[n]]
        #        if (sum(is.na(tempPeaks[,1])) > dim(tempPeaks)[1] * (3 / 4) )
        if (sum(!is.na(tempPeaks[,1])) <
            sum(groupInfo[9:length(groupInfo)]) * (1 / 4) )
            rmInd[n] <- FALSE
    }
    peakGroup <- peakGroup[rmInd]
    groupNum <- length(peakGroup)
    peakGroupMed <- tempPeaks[seq_len(groupNum), , drop = FALSE]
    for (n in seq_len(groupNum))
        for (m in seq_len(dim(peakGroupMed)[2]))
            peakGroupMed[n, m] <- median(peakGroup[[n]][ , m], TRUE)

    ## Fill in missing peaks for each group
    tempLabel <- c('scan_ncGTW', 'rt_ncGTW', 'subscan_ncGTW','scanmin_ncGTW',
                    'rtmin_ncGTW', 'scanmax_ncGTW', 'rtmax_ncGTW')
    for (n in seq_len(groupNum)){
        tempPeaks <- peakGroup[[n]]
        for (m in seq_len(dim(tempPeaks)[1])){
            if (!is.na(tempPeaks[m, 1]))
                next

            tempPeaks[m , tempLabel] <- peakGroupMed[n, tempLabel]
            samPath <- path[[m]]
            apexInd <- samPath[which(samPath[ , 1] ==
                round(peakGroupMed[n, 'subscan_ncGTW'])), 2]

            apexInd <- apexInd[which.max(data[m, apexInd])]

            tempPeaks[m, 'scan_Raw'] <-
                which(scanRangeOld[m, ] == scanRange[m, apexInd])
            tempPeaks[m, 'scanmin_Raw'] <- peakGroupMed[n, 'scanmin_ncGTW'] +
                tempPeaks[m, 'scan_Raw'] - peakGroupMed[n, 'scan_ncGTW']
            tempPeaks[m, 'scanmax_Raw'] <- peakGroupMed[n, 'scanmax_ncGTW']+
                tempPeaks[m, 'scan_Raw'] - peakGroupMed[n, 'scan_ncGTW']
            tempPeaks[m, 'sample'] <- m
        }
        peakGroup[[n]] <- tempPeaks
    }


    ##  Warp RT of peaks of each samples
    rtncGTW <- rtXCMS
    rtncGTWsub <- scanRangeOld * NA
    rtRawSub <- scanRangeOld * NA
    for (n in seq_len(length(rtRaw)))
        rtRawSub[n, ] <- rtRaw[[n]][scanRangeOld[n, ]]

    tempShift <- 5
    while (any(t(diff(t(rtRawSub))) > 5)){
        for (n in seq_along(rtRaw)){
            rtRawSub[n, ] <- rtRaw[[n]][scanRangeOld[n, ] - tempShift]
        }
        tempShift <- tempShift + 5
    }

    for (n in seq_len(groupNum)){
        tempPeaks <- peakGroup[[n]]
        for (m in seq_len(dim(tempPeaks)[1])){
            rawSta <- tempPeaks[m, 'scanmin_Raw']
            rawEnd <- tempPeaks[m, 'scanmax_Raw']
            ncGTWSta <- tempPeaks[m, 'scanmin_ncGTW']
            ncGTWEnd <- tempPeaks[m, 'scanmax_ncGTW']

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
            if (sum(overlapFlag) == 0)
                next

            ### Solve overlapping peaks
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

            while (rtncGTWsub[m, olSta] > rtncGTWsub[m, olEnd]){
                olStaOld <- olSta
                olEndOld <- olEnd
                if ((olSta - 1 > 0) && !is.na(rtncGTWsub[m, olSta - 1]))
                    olSta <- olSta - 1

                if ((olEnd + 1 < dim(rtncGTWsub)[2]) &&
                    !is.na(rtncGTWsub[m, olEnd + 1]))
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
            rtncGTWsub[m, olSta:olEnd] <- approx(c(olSta, olEnd),
                rtncGTWsub[m, c(olSta, olEnd)], n = length(olSta:olEnd))$y

        }
    }

    ## Filling non-peak part of each sample
    for (n in seq_len(dim(rtncGTWsub)[1])){
        while (!is.na(which(is.na(rtncGTWsub[n, ]))[1])){
            ipSta <- which(is.na(rtncGTWsub[n, ]))[1]
            ipEnd <- which(!is.na(rtncGTWsub[n, ipSta:dim(rtncGTWsub)[2]]))[1]
            if (is.na(ipEnd)){
                ipEnd <- dim(rtncGTWsub)[2]
            } else {
                ipEnd <- ipEnd + ipSta - 2
            }

            if (ipSta == 1){
                ipStaRt <- rtncGTW[[n]][max(scanRangeOld[n, ipSta] - 1, 1)]
            } else{
                ipStaRt <- rtncGTWsub[n, ipSta - 1]
            }

            if (ipEnd == dim(rtncGTWsub)[2]){
                tempEnd <- min(length(rtncGTW[[n]]), scanRangeOld[n, ipEnd] + 1)
                ipEndRt <- rtncGTW[[n]][tempEnd]
            } else{
                ipEndRt <- rtncGTWsub[n, ipEnd + 1]
            }

            if (ipStaRt > ipEndRt){
                tempDif <- median(diff(rtncGTW[[n]]))
                if (ipSta == 1){
                    ipStaRt <- ipEndRt - (ipEnd - ipSta) * tempDif
                } else{
                    ipEndRt <- ipStaRt + (ipEnd - ipSta) * tempDif
                }
            }
            ipRt <- approx(c(ipSta - 1, ipEnd + 1), c(ipStaRt, ipEndRt),
                            n = length((ipSta - 1):(ipEnd + 1)))$y
            rtncGTWsub[n, ipSta:ipEnd] <- ipRt[2:(length(ipRt) - 1)]
        }
    }

    ## Fix decreasing RT
    for (n in seq_len(dim(rtncGTWsub)[1])){
        for (m in 2:dim(rtncGTWsub)[2]){
            if (rtncGTWsub[n, m] >= rtncGTWsub[n, m - 1])
                next

            pl <- m - 1
            pr <- m
            while (rtncGTWsub[n, pl] > rtncGTWsub[n, pr]){
                if (pl > 1)
                    pl <- pl - 1
                if (pr < dim(rtncGTWsub)[2])
                    pr <- pr + 1
            }
            rtncGTWsub[n, pl:pr] <- approx(c(pl, pr),
                rtncGTWsub[n, c(pl, pl)], n = length(pl:pr))$y
        }
    }
    for (n in seq_len(dim(rtncGTWsub)[1]))
        rtncGTW[[n]][scanRangeOld[n,]] <- rtncGTWsub[n, ]

    return(new("ncGTWwarp", rtncGTW=rtncGTW, ncGTWpeaks=peaks))
}


findUniPeak <- function(peaks, groupInd, groupidx, ppm=1e6, sampleInd=NULL){
    ppm <- ppm / 1e6

    Peaks <- cbind(peaks[groupidx[[groupInd]], ],
        peakInd=groupidx[[groupInd]])
    if (!is.null(sampleInd))
        Peaks <- Peaks[is.element(Peaks[ , 'sample'], sampleInd), , drop=FALSE]

    if (dim(Peaks)[1] == 0)
        return(Peaks)

    tempInd <- sort(Peaks[ , 'maxo'], index.return = TRUE, decreasing = TRUE)$ix
    Peaks <- Peaks[tempInd, ]
    Peaks <- Peaks[!duplicated(Peaks[,c('rt', 'sample')]), , drop = FALSE]
    tempInd <- sort(Peaks[ , 'peakInd'], index.return = TRUE)$ix
    Peaks <- Peaks[tempInd, ]

    rmInd <- matrix(TRUE, dim(Peaks)[1], 1)

    for (n in seq_len(dim(Peaks)[1] - 1)){
        if (rmInd[n] == FALSE)
            next
        if (abs(Peaks[n, 'rtmax'] - Peaks[n, 'rtmin']) > 30){
            rmInd[n] <- FALSE
            next
        }
        samPeakInd <-
            which(Peaks[ , 'sample'] == Peaks[n, 'sample'] & rmInd == TRUE)
        samPeakInd <- setdiff(samPeakInd, n)
        samPeaks <- cbind(Peaks[samPeakInd, , drop=FALSE])

        for (m in seq_len(dim(samPeaks)[1])){
            if (abs(Peaks[n, 'mz'] - samPeaks[m, 'mz']) / Peaks[n, 'mz'] >= ppm)
                next
            if (abs(Peaks[n, 'rt'] - samPeaks[m, 'rt']) >= 3)
                next

            if (Peaks[n, 'maxo'] < samPeaks[m, 'maxo']){
                rmInd[n] <- FALSE
                break
            } else{
                rmInd[samPeakInd[m]] <- FALSE
            }
        }
    }
    return(Peaks[rmInd, , drop = FALSE])
}

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
    ia <- matrix(is.element(ee[ ,c(1,2)], ct), dim(ee)[1], 2)
    dtmp[ia] <- 1
    isCutEE <- dtmp[ ,1] != dtmp[ , 2]

    # cut to mapping pattern in primal graph
    resPath <- vector('list', nPix)
    for (nn in seq_len(nPix)){
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
        ia <- is.element(dEdgeIntSS[ , 2], idxSrcCut)
        # get corresponding edges, src -> node
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
    p0 <- path[ , c(1,2)]
    idxValid <- (p0[ , 1] >= 1) & (p0[ , 1] <= nTps) & (p0[ , 2] >= 1) &
        (p0[ , 2] <= nTps)
    p0 <- p0[idxValid, ]

    for (tt in seq_len(dim(p0)[1])){
        pRef <- p0[tt, 2]
        pTst <- p0[tt, 1]
        warped[pTst] <- max(warped[pTst], curve[pRef])
    }
    return(warped)
}


pathCombine <- function(parPath, path2, parInd){
    dataNum <- max(parInd)
    temp2path <- vector('list', dataNum)

    for (tempInd in seq_len(dataNum)){
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

        newPath <- cbind(newPath, rbind(newPath[2:dim(newPath)[1], ],
                                        tempPath2[dim(tempPath2)[1], 3:4]))
        temp2path[[tempInd]] <- newPath
    }
    return(temp2path)
}
