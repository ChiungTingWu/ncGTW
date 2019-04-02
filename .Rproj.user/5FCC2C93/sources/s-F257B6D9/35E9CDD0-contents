peakGroupPvalOrder <- function(peakIdx, sampleNum) {
    peakIdx <- unique(sort(peakIdx))
    peakIdxDif <- peakIdx[length(peakIdx)] - peakIdx[1]

    hypValue <- 0
    for (n in (length(peakIdx)-1) : peakIdxDif)
    {
        hypValue <- hypValue + (sampleNum - n) * choose(n-1, length(peakIdx)-2)
    }
    pVal <- hypValue / choose(sampleNum, length(peakIdx))
    if (length(peakIdx) == 1){
        pVal <- 1
    }
    return(pVal)
}


peakGroupPval <- function(peakIdx, sampleNum, testNum) {
    peakIdx <- unique(sort(peakIdx))
    peakIdxDif <- peakIdx[2:length(peakIdx)] - peakIdx[1:(length(peakIdx) - 1)]
    testVal <- mean(peakIdxDif)

    peakNum <- length(peakIdx)
    nullVal <- vector('numeric', testNum)

    for (n in 1:testNum){
        nullIdx <- sort(sample(1:sampleNum, peakNum))
        nullIdxDif <- nullIdx[2:length(nullIdx)] - nullIdx[1:(length(nullIdx) - 1)]
        nullVal[n] <- mean(nullIdxDif)
    }

    pVal <- sum(nullVal <= testVal) / length(nullVal)

    if (length(peakIdx) == 1){
        pVal <- 1
    }

    return(pVal)
}

splitGroupPval <- function(xcmsLargeWin, xcmsSmallWin) {

    largeWinGroupidx <- xcmsLargeWin@groupidx
    smallWinGroupidx <- xcmsSmallWin@groupidx

    groupLargeSmallIndex <- matrix(0, nrow(xcmsLargeWin@peaks), 2)

    for (i in 1:length(largeWinGroupidx)) {
        for (j in 1:length(largeWinGroupidx[[i]])){
            groupLargeSmallIndex[largeWinGroupidx[[i]][j], 1] <- i
        }
    }
    for (i in 1:length(smallWinGroupidx)) {
        for (j in 1:length(smallWinGroupidx[[i]])){
            groupLargeSmallIndex[smallWinGroupidx[[i]][j], 2] <- i
        }
    }

    matchTwoGroup <- matrix(0, length(largeWinGroupidx), length(smallWinGroupidx))
    for (n in 1:nrow(xcmsLargeWin@peaks)) {
        if (groupLargeSmallIndex[n,1] != 0 & groupLargeSmallIndex[n,2] != 0) {
            matchTwoGroup[groupLargeSmallIndex[n,1], groupLargeSmallIndex[n,2]] =
                matchTwoGroup[groupLargeSmallIndex[n,1], groupLargeSmallIndex[n,2]] + 1
        }
    }

    matchGroupCount <- apply(matchTwoGroup, 1, function(x) sum(x != 0))
    matchGroupIdx <- apply(matchTwoGroup, 1, function(x) which(x != 0))

    match2More <- which(matchGroupCount>1)

    sampleNum <- length(xcmsLargeWin@filepaths)


    largeWinSampleidx <- vector('list', length(largeWinGroupidx))
    for (n in 1:length(largeWinGroupidx)){
        largeWinSampleidx[[n]] <- xcmsLargeWin@peaks[largeWinGroupidx[[n]], 'sample']
    }

    smallWinSampleidx <- vector('list', length(smallWinGroupidx))
    for (n in 1:length(smallWinGroupidx)){
        smallWinSampleidx[[n]] <- xcmsSmallWin@peaks[smallWinGroupidx[[n]], 'sample']
    }

    pvalues <- vector('integer', length(smallWinGroupidx)) - 1
    pvalueIdx <- vector('integer', 0)

    for (n in 1:length(smallWinGroupidx)){
        if (sum(matchTwoGroup[,n]) > 0){
            #pvalues[n] <- peakGroupPval(smallWinSampleidx[[n]], sampleNum, testNum)
            pvalues[n] <- peakGroupPvalOrder(smallWinSampleidx[[n]], sampleNum)
            pvalueIdx <- c(pvalueIdx, n)
        }
    }

#    peakMore1Idx <- which(lapply(smallWinSampleidx, length) > 1)

#    tempPval <- pvalues[peakMore1Idx]
    tempPval <- pvalues[pvalueIdx]
    tempQval <- p.adjust(tempPval, 'fdr')
    qvalues <- vector('integer', length(smallWinGroupidx)) + 1
#    qvalues[peakMore1Idx] <- tempQval
    qvalues[pvalueIdx] <- tempQval

    pqValues <- matrix(0, length(largeWinGroupidx), 2)

    for (n in 1:length(largeWinGroupidx)){
        tempIdx <- matchGroupIdx[[n]]
        if (length(tempIdx) == 0){
            pqValues[n,1:2] <- 2
        } else{
            pqValues[n, 1] <- min(pvalues[tempIdx])
            pqValues[n, 2] <- min(qvalues[tempIdx])
        }
    }
    groupInfo <- list(largeWinGroups = xcmsLargeWin@groups, smallWinGroups = xcmsSmallWin@groups,
                      largeWinGroupidx = largeWinGroupidx, smallWinGroupidx = smallWinGroupidx,
                      largeWinSampleidx = largeWinSampleidx, smallWinSampleidx = smallWinSampleidx,
                      matchGroupIdx = matchGroupIdx, pvalueIdx = pvalueIdx,
                      pvalues = pvalues, qvalues = qvalues, pqValues = pqValues)

    return(groupInfo)
}

exclusiveGroups <- function(groupInfo, ppm){
    excluGroupsNum <- vector('integer', length(groupInfo$largeWinGroupidx))
    excluLargeSmall <- vector('list', length(groupInfo$largeWinGroupidx))
    excluPval <- vector('list', length(groupInfo$largeWinGroupidx))

    groups <- groupInfo$smallWinGroups

    for (n in 1:length(groupInfo$largeWinGroupidx)){
        tempInd <- groupInfo$matchGroupIdx[[n]]
        tempInd <- tempInd[groupInfo$qvalues[tempInd] < 0.05]

        if (length(tempInd) > 1){
            tempGroups <- groupInfo$smallWinSampleidx[tempInd]
            tempGroupNum <- length(tempGroups)

            for (i in 1:(tempGroupNum - 1)){
                for (j in (i+1):tempGroupNum){
                    if (length(intersect(tempGroups[[i]], tempGroups[[j]])) <= 0){
                        if (abs(groups[tempInd[i], 'mzmed'] - groups[tempInd[j], 'mzmed']) /
                            min(groups[tempInd[i], 'mzmed'], groups[tempInd[j], 'mzmed']) *
                            1000000 < ppm){
                            excluGroupsNum[n] <- excluGroupsNum[n] + 1
                            ttest <- t.test(tempGroups[[i]], tempGroups[[j]])
                            excluPval[[n]] <- c(excluPval[[n]], ttest$p.value)
                            excluLargeSmall[[n]] <- cbind(excluLargeSmall[[n]],
                                                          rbind(tempInd[i], tempInd[j]))


                        }

                    }
                }
            }
        }
    }
    excluGroups <- cbind(which(excluGroupsNum != 0),
                         groupInfo$largeWinGroups[which(excluGroupsNum != 0),])
    colnames(excluGroups)[1] <- 'index'

    excluInfo <- list(excluGroups = excluGroups, excluPval = excluPval,
                      excluLargeSmall = excluLargeSmall)

    return(excluInfo)
}

misalignDetect <- function(xcmsLargeWin, xcmsSmallWin, ppm){
    groupInfo <- splitGroupPval(xcmsLargeWin, xcmsSmallWin)
    excluInfo <- exclusiveGroups(groupInfo, ppm)
    excluGroups <- excluInfo$excluGroups
    #keepIdx <- vector('logical', dim(excluGroups)[1])
    #for (n in 1:dim(excluGroups)[1]){
    #    if (min(excluInfo$excluPval[[excluGroups[n,'index']]]) <= 0.05){
    #        keepIdx[n] = TRUE
    #    }
    #}
    #excluGroups <- excluGroups[keepIdx,]
    return(excluGroups)
}
