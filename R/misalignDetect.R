#' Detect misaligned peak groups in xcmsSet object of XCMS
#'
#' This function detects the misaligned peak groups with two
#' \code{\link[xcms]{xcmsSet-class}} object with two different values of
#' \code{bw} parameter in \code{\link[xcms:group.density]{group}}.
#' @param xcmsLargeWin A \code{\link[xcms]{xcmsSet-class}} object with a larger
#'   \code{bw}, usually the maximum expected retension time drift.
#' @param xcmsSmallWin A \code{\link[xcms]{xcmsSet-class}} object with a smaller
#'   \code{bw}, usually the resolution of the retension time.
#' @param ppm Should be set as same as the one when performing the peak
#'   detection function in \code{xcms}.
#' @param qThre The threshould of the p-value after multiple test correction.
#'   The default is 0.05.
#' @param maxRtWin The threshould of the maximum retension time range. This is
#'   for filtering out some bad groups. The default is 50 (seconds).
#' @details This function includes two major steps to determine a peak group is
#' misaligned or not. The first step calculates the p-value of each peak group
#' in xcmsSmallWin, and find the corresponding peak group in xcmsLargeWin. The
#' second step is to find the exclusive peak groups (the groups with no
#' overlapping samples) with adjsted p-values smaller than \code{qThre}.
#' @return A matrix with all detected misaligned peak groups. The column names
#' are the same as \code{group} slot in \code{\link[xcms]{xcmsSet-class}}, but
#' the first column is the group index.
#' @examples
#' # obtain data
#' data('xcmsExamples')
#' xcmsLargeWin <- xcmsExamples$xcmsLargeWin
#' xcmsSmallWin <- xcmsExamples$xcmsSmallWin
#' ppm <- xcmsExamples$ppm
#'
#' # detect misaligned features
#' excluGroups <- misalignDetect(xcmsLargeWin, xcmsSmallWin, ppm)
#' @export

misalignDetect <- function(xcmsLargeWin, xcmsSmallWin, ppm, qThre=0.05,
    maxRtWin=50){
    groupInfo <- splitGroupPval(xcmsLargeWin, xcmsSmallWin)
    excluInfo <- exclusiveGroups(groupInfo, ppm, qThre)
    excluGroups <- excluInfo$excluGroups

    rtWin <-  excluGroups[ , 'rtmax'] - excluGroups[ , 'rtmin']
    return(excluGroups[rtWin < maxRtWin, , drop = FALSE])
}




splitGroupPval <- function(xcmsLargeWin, xcmsSmallWin) {

    largeWinGroupidx <- xcmsLargeWin@groupidx
    smallWinGroupidx <- xcmsSmallWin@groupidx

    groupLargeSmallIndex <- matrix(0, nrow(xcmsLargeWin@peaks), 2)

    for (i in seq_along(largeWinGroupidx)) {
        for (j in seq_along(largeWinGroupidx[[i]])){
            groupLargeSmallIndex[largeWinGroupidx[[i]][j], 1] <- i
        }
    }
    for (i in seq_along(smallWinGroupidx)) {
        for (j in seq_along(smallWinGroupidx[[i]])){
            groupLargeSmallIndex[smallWinGroupidx[[i]][j], 2] <- i
        }
    }

    matchTwoGroup <- matrix(0, length(largeWinGroupidx),
                            length(smallWinGroupidx))
    for (n in seq_len(nrow(xcmsLargeWin@peaks))) {
        if (groupLargeSmallIndex[n,1] != 0 &
            groupLargeSmallIndex[n,2] != 0) {
            matchTwoGroup[groupLargeSmallIndex[n,1],
            groupLargeSmallIndex[n,2]] =
            matchTwoGroup[groupLargeSmallIndex[n,1],
            groupLargeSmallIndex[n,2]] + 1
        }
    }

    matchGroupCount <- apply(matchTwoGroup, 1, function(x) sum(x != 0))
    matchGroupIdx <- apply(matchTwoGroup, 1, function(x) which(x != 0))

    match2More <- which(matchGroupCount>1)

    sampleNum <- length(xcmsLargeWin@filepaths)


    largeWinSampleidx <- vector('list', length(largeWinGroupidx))
    for (n in seq_along(largeWinGroupidx)){
        largeWinSampleidx[[n]] <-
            xcmsLargeWin@peaks[largeWinGroupidx[[n]], 'sample']
    }

    smallWinSampleidx <- vector('list', length(smallWinGroupidx))
    for (n in seq_along(smallWinGroupidx)){
        smallWinSampleidx[[n]] <-
            xcmsSmallWin@peaks[smallWinGroupidx[[n]], 'sample']
    }

    pvalues <- vector('integer', length(smallWinGroupidx)) - 1
    pvalueIdx <- vector('integer', 0)

    for (n in seq_along(smallWinGroupidx)){
        if (sum(matchTwoGroup[,n]) > 0){
            pvalues[n] <- peakGroupPvalOrder(smallWinSampleidx[[n]], sampleNum)
            pvalueIdx <- c(pvalueIdx, n)
        }
    }


    tempPval <- pvalues[pvalueIdx]
    tempQval <- p.adjust(tempPval, 'fdr')
    qvalues <- vector('integer', length(smallWinGroupidx)) + 1

    qvalues[pvalueIdx] <- tempQval

    pqValues <- matrix(0, length(largeWinGroupidx), 2)

    for (n in seq_along(largeWinGroupidx)){
        tempIdx <- matchGroupIdx[[n]]
        if (length(tempIdx) == 0){
            pqValues[n,c(1,2)] <- 2
        } else{
            pqValues[n, 1] <- min(pvalues[tempIdx])
            pqValues[n, 2] <- min(qvalues[tempIdx])
        }
    }
    groupInfo <-
        list(largeWinGroups = xcmsLargeWin@groups,
            smallWinGroups = xcmsSmallWin@groups,
            largeWinGroupidx = largeWinGroupidx,
            smallWinGroupidx = smallWinGroupidx,
            largeWinSampleidx = largeWinSampleidx,
            smallWinSampleidx = smallWinSampleidx,
            matchGroupIdx = matchGroupIdx, pvalueIdx = pvalueIdx,
            pvalues = pvalues, qvalues = qvalues, pqValues = pqValues)

    return(groupInfo)
}


exclusiveGroups <- function(groupInfo, ppm, qThre){
    excluGroupsNum <- vector('integer', length(groupInfo$largeWinGroupidx))
    excluLargeSmall <- vector('list', length(groupInfo$largeWinGroupidx))
    excluPval <- vector('list', length(groupInfo$largeWinGroupidx))

    groups <- groupInfo$smallWinGroups

    for (n in seq_along(groupInfo$largeWinGroupidx)){
        tempInd <- groupInfo$matchGroupIdx[[n]]
        tempInd <- tempInd[groupInfo$qvalues[tempInd] < qThre]

        if (length(tempInd) > 1){
            tempGroups <- groupInfo$smallWinSampleidx[tempInd]
            tempGroupNum <- length(tempGroups)

            for (i in seq_len(tempGroupNum - 1)){
                for (j in seq.int(i+1, tempGroupNum)){
                    if (length(intersect(tempGroups[[i]],
                        tempGroups[[j]])) <= 0){
                        if (abs(groups[tempInd[i], 'mzmed'] -
                                groups[tempInd[j], 'mzmed']) /
                            min(groups[tempInd[i], 'mzmed'],
                                groups[tempInd[j], 'mzmed']) *
                            1000000 < ppm){
                            excluGroupsNum[n] <- excluGroupsNum[n] + 1
                            ttest <- t.test(tempGroups[[i]], tempGroups[[j]])
                            excluPval[[n]] <- c(excluPval[[n]], ttest$p.value)
                            excluLargeSmall[[n]] <-
                                cbind(excluLargeSmall[[n]],
                                rbind(tempInd[i], tempInd[j]))


                        }

                    }
                }
            }
        }
    }
    excluGroups <-
        cbind(which(excluGroupsNum != 0),
            groupInfo$largeWinGroups[which(excluGroupsNum != 0), ,drop=FALSE])
    colnames(excluGroups)[1] <- 'index'

    excluInfo <-
        list(excluGroups = excluGroups, excluPval = excluPval,
            excluLargeSmall = excluLargeSmall)

    return(excluInfo)
}
















peakGroupPvalOrder <- function(peakIdx, sampleNum) {
    peakIdx <- unique(sort(peakIdx))
    peakIdxDif <- peakIdx[length(peakIdx)] - peakIdx[1]

    hypValue <- 0
    for (n in (length(peakIdx)-1) : peakIdxDif)
    {
        hypValue <- hypValue + (sampleNum - n) * choose(n-1, length(peakIdx)-2)
    }
    pVal <- hypValue / choose(sampleNum, length(peakIdx))
    if (length(peakIdx) == 1 | is.nan(pVal)){
        pVal <- 1
    }
    return(pVal)
}


peakGroupPval <- function(peakIdx, sampleNum, testNum) {
    peakIdx <- unique(sort(peakIdx))
    peakIdxDif <- peakIdx[seq.int(2, length(peakIdx))] -
        peakIdx[seq_len(length(peakIdx) - 1)]
    testVal <- mean(peakIdxDif)

    peakNum <- length(peakIdx)
    nullVal <- vector('numeric', testNum)

    for (n in seq_len(testNum)){
        nullIdx <- sort(sample(seq_len(sampleNum), peakNum))
        nullIdxDif <- nullIdx[seq.int(2, length(nullIdx))] -
            nullIdx[seq_len(length(nullIdx) - 1)]
        nullVal[n] <- mean(nullIdxDif)
    }

    pVal <- sum(nullVal <= testVal) / length(nullVal)

    if (length(peakIdx) == 1){
        pVal <- 1
    }
    return(pVal)
}
