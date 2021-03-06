buildMultiParaValidmap <- function(samples, mir, strNum, diaNum, biP=TRUE) {
    if (!(is.logical(mir) && is.logical(biP)))
        stop("mir and biP should be 'TRUE' or 'FALSE'!!!")

    sampNum <- nrow(samples)
    sampLen <- ncol(samples)
    pairMap <- matrix(0, sampNum, sampNum)
    pairNum <- 0

    ## build tst, ref, and pairMap depending on mir & biP
    if (mir && biP){
        tst <- matrix(0, sampNum * (sampNum - 1) / 2, sampLen)
        ref <- tst

        for (jj in seq_len(sampNum)){
            for (ii in seq(jj + 1, sampNum, length=max(0, sampNum - jj))){
                pairNum <- pairNum + 1
                pairMap[ii, jj] <- pairNum
                tst[pairNum, ] <- samples[ii, ]
                ref[pairNum, ] <- samples[jj, ]
            }
        }
    } else{
        tst <- matrix(0, sampNum * (sampNum - 1), sampLen)
        ref <- tst

        for (jj in seq_len(sampNum)){
            for (ii in seq_len(sampNum)){
                if (ii == jj)
                    next
                pairNum <- pairNum + 1
                pairMap[ii, jj] <- pairNum
                tst[pairNum, ] <- samples[ii, ]
                ref[pairNum, ] <- samples[jj, ]
            }
        }
    }

    ## build validMap
    validMap <- matrix(0, pairNum, pairNum)
    #eachMap <- array(0, c(sampNum, sampNum, pairNum))
    for (jj in seq_len(sampNum)){
        for (ii in seq_len(sampNum)){
            if (pairMap[ii, jj] == 0)
                next

            ## check vertical elements
            verEnd <- min(ii + strNum, sampNum)
            for (iii in seq(ii + 1, verEnd, length=max(0, verEnd - ii))){
                if (pairMap[iii, jj] == 0)
                    next

                if (!(mir && !((ii - jj) * (iii - jj) > 0)))
                    validMap[pairMap[iii, jj], pairMap[ii, jj]] <- 1
            }

            ## check horizontal elements
            horEnd <- min(jj + strNum, sampNum)
            for (jjj in seq(jj + 1, horEnd, length=max(0, horEnd - jj))){
                if (pairMap[ii, jjj] == 0)
                    next

                if (!(mir && !((ii - jj) * (ii - jjj) > 0)))
                    validMap[pairMap[ii, jjj], pairMap[ii, jj]] <- 1
            }

            ## check cross elements
            if (!mir){
                if ((ii - jj == 1) && (jj + 2 <= sampNum))
                    validMap[pairMap[ii, jj + 2], pairMap[ii, jj]] <- 1

                if ((jj - ii == 1) && (ii + 2 <= sampNum))
                    validMap[pairMap[ii + 2, jj], pairMap[ii, jj]] <- 1
            }

            ## check diagonal elements
            diaEnd <- min(diaNum, sampNum - max(ii, jj))
            for (kkk in seq(1, diaEnd, length=max(0, diaEnd)))
                validMap[pairMap[ii + kkk, jj + kkk], pairMap[ii, jj] ] <- 1
        }
    }
    return(list(tst=tst, ref=ref, pairMap=pairMap, validMap=validMap))
}


initGtwParam <- function(validMap, noiseVar, maxStp, smo=1, diagonalPenalty=0,
    logt=0, nor=2) {

    tmpEast <- validMap * 1
    tmpSouth <- validMap * 1
    smoMap <- array(c(tmpEast, tmpSouth) * smo, c(dim(validMap), 2))

    s2 <- noiseVar  # noise variance map
    epAdd <- 0
    epMul <- 1
    winSize <- max(maxStp, 1)
    offDiagonalPenalty <- 0
    inf0 <- 1e300
    partialMatchingCost <- (0:winSize) * 1e8

    return(list(smoBase=smo, diagonalPenalty=diagonalPenalty, logt=logt,
        nor=nor, smoMap=smoMap, s2=s2, epAdd=epAdd, epMul=epMul,
        winSize=winSize, offDiagonalPenalty=offDiagonalPenalty, inf0=inf0,
        partialMatchingCost=partialMatchingCost))
}


getDistMat <- function(ref, tst, refP, tstP, biP, weiP, nor){

    nTps <- length(ref)
    refM <- replicate(nTps, ref)
    tstM <- t(replicate(nTps, tst))

    refPM <- replicate(nTps, refP)
    tstPM <- replicate(nTps, tstP)

    d0 <- (abs(refM - tstM)) ^ nor
    if ((nor != 1) && ((nor != 2)))
        warning('norm is not an integer!')

    if (sum(is.na(d0))){
        warning('There are some NA in d0.')
        d0[is.na(d0)] <- 0
    }

    stopifnot(is.logical(biP), !is.na(biP))
    if (weiP == -1) stopifnot(biP)
    switch(as.character(weiP),
        "0"={ dp <- d0 },
        "-1"={
            a0 <- (tstPM + refPM) / 2
            g0 <- (tstPM * refPM) ^ (0.5)
            a1 <- (a0 + g0) / 2
            g1 <- (a0 * g0) ^ (0.5)
            a2 <- (a1 + g1) / 2
            g2 <- (a1 * g1) ^ (0.5)
            a3 <- (a2 + g2) / 2
            dp <- d0 / a3
        },
        dp <- ifelse(biP, list((d0 / ((tstPM * refPM)^weiP))),
            list(d0 / (tstPM ^ weiP)))[[1]]
    )

    dp[is.infinite(dp)] <- 10 ^ 300
    dp[is.infinite(dp)] <- 10 ^ 300

    dM <- array(0, c(dim(d0), 2))
    dM[, , 1] <- d0
    dM[, , 2] <- dp
    return(dM)
}

buildGTWgraph <- function(ref, tst, validMap, param, mu=0, sigma=1, biP=TRUE,
    weiP=0, type=0, path=NA, pairMap=NA) {

    ## type=0: typical GTW graph
    ## type=1: ncGTW graph
    ## type=2: merge GTW graph

    smoMap <- param$smoMap
    smoBase <- param$smoBase
    s2 <- param$s2
    epAdd <- param$epAdd
    epMul <- param$epMul
    win <- param$winSize
    ofstPen <- param$offDiagonalPenalty
    diagonalPenalty <- param$diagonalPenalty
    capRev <- param$inf0
    pmCost <- param$partialMatchingCost

    logt <- param$logt
    nor <- param$nor

    if (sum(abs(diag(validMap))) != 0)
        stop('The diagonal terms of validMap should all be zero.')

    nPix <- dim(tst)[1]
    nTps <- dim(tst)[2]

    if (length(mu) != length(sigma))
        stop('The number of mu, sigma and pairs should be the same.')
    if ((length(mu) != 1) && (length(mu) != nPix))
        stop('The number of mu should be 1 or number of pairs.')

    if (logt == 0){
        refP <- 1 - pnorm(log(ref+1), mu, sigma)
        tstP <- 1 - pnorm(log(tst+1), mu, sigma)
    } else{
        refP <- 1 - pnorm(ref, mu, sigma)
        tstP <- 1 - pnorm(tst, mu, sigma)
    }


    distMask <- matrix(0, nTps, nTps)
    diagMask <- matrix(0, nTps, nTps)
    for (ii in seq_len(nTps)){
        diagMask[ii,ii] <- diagonalPenalty
        for (jj in seq_len(nTps))
            distMask[ii, jj] <- ofstPen * (abs(ii - jj) > 0)
    }

    scl0 <- 1e8  # larger than time points
    #  validMapIdx <- validMap * 0
    #  validMapIdx[validMap > 0] <- 1 : sum(validMap > 0)

    nNeibPair <- (sum(lower.tri(validMap) * validMap != 0))

    if (length(s2) == 1)
        s2 <- matrix(s2, nPix, 1)
    s2x <- s2[s2 > 0]
    s2xMean <- mean(s2x, na.rm=TRUE)
    pmCost <- pmCost / 2 / s2xMean

    ## template using coordinate
    edgeInfo  <- buildPairTemplate(nTps,win,pmCost)
    pEdge     <- (edgeInfo$pEdge)
    dEdge     <- (edgeInfo$dEdge)
    weightPos <- (edgeInfo$cPos)
    weightVal <- (edgeInfo$cVal)
    st01      <- (edgeInfo$st01)
    eType     <- (edgeInfo$eType)


    ## -- direction penalty, additive or multiplicative
    weightVal[eType == 2] <- weightVal[eType == 2] + epAdd
    weightValMul <- weightVal * 0 + 1
    weightValMul[eType == 2] <- epMul

    ## re-code coordinates to single number
    tmp <- dEdge * 4
    tmp <- cbind(tmp[ , 1] * scl0 + tmp[ , 2], tmp[ , 3] * scl0 + tmp[ , 4])
    ia <- which(!(duplicated(as.vector(tmp))))
    tmpUniq <- tmp[ia]
    tmpUniqSort <- sort(tmpUniq, index.return = TRUE)
    tmpUniq <- tmpUniqSort$x
    ia <- ia[tmpUniqSort$ix]

    nNodeGrid <- length(tmpUniq) - 2
    srcNode <- nNodeGrid * nPix + 1
    sinkNode <- nNodeGrid * nPix + 2

    mapObj <- as.list(c(sinkNode, seq_len(nNodeGrid), srcNode))
    names(mapObj) <- tmpUniq

    dEdgeInt <- apply(tmp, c(1, 2), function(x) mapObj[[as.character(x)]])
    dEdge1 <- rbind(dEdge[ ,c(1, 2)], dEdge[ , 3:4])
    nodePos <- dEdge1[ia[2:(length(ia) - 1)], ]

    ## split the edge matrix to src/sink and within grid
    ## nodes connected with src or sink
    d1 <- dEdgeInt[ , 1]
    d2 <- dEdgeInt[ , 2]
    idxSrc <- which(d1 == srcNode)
    idxSink <- which(d2 == sinkNode)

    ## extra weight for ss edges
    ssTmpWt <- matrix(0, nNodeGrid, 2)
    ssTmpWt[d2[idxSrc], 1]  <- weightVal[idxSrc]
    ssTmpWt[d1[idxSink], 2] <- weightVal[idxSink]

    ## weight from curve distance for ss edges
    wtPos1 <- cbind(weightPos[,1] + (weightPos[,2] - 1) * nTps, eType)
    ssTmpPos <- matrix(0, nNodeGrid, 2) + nTps * nTps + 1
    ssTmpPos[d2[idxSrc],  1] <- wtPos1[idxSrc,  1]
    ssTmpPos[d1[idxSink], 2] <- wtPos1[idxSink, 1]

    ## edges between nodes (except src and sink)
    idxNotSS <- d1 != srcNode & d2 != sinkNode
    eeTmp <- cbind(dEdgeInt[idxNotSS, ], wtPos1[idxNotSS, ],
        weightVal[idxNotSS], weightValMul[idxNotSS])
    nEdgeGrid <- dim(eeTmp)[1]

    ## output, for label and path mapping
    pEdgeSS <- pEdge[!idxNotSS, ]
    pEdgeEE <- pEdge[idxNotSS, ]
    dEdgeIntSS <- dEdgeInt[!idxNotSS, ]
    dEdgeIntEE <- dEdgeInt[idxNotSS, ]

    ## edges for within pairs using integer
    if (type == 2){
        ssPair <- matrix(0, nNodeGrid, 2)
        eePair <- matrix(0, nEdgeGrid, 4)
    } else if (smoBase > 0 && type == 0){
        ssPair <- matrix(0, nPix * nNodeGrid, 2)
        eePair <- matrix(0, nPix * nEdgeGrid + nNodeGrid * nNeibPair, 4)
    } else{
        ssPair <- matrix(0, nPix * nNodeGrid, 2)
        eePair <- matrix(0, nPix * nEdgeGrid, 4)
    }
    eePair[ , 4] <- capRev

    for (ii in seq_len(nPix)){
        # s2x <- s2[ii]
        if (type == 2){
            eeOfst <- 0
            ssOfst <- 0
        } else{
            eeOfst <- (ii - 1) * nEdgeGrid
            ssOfst <- (ii - 1) * nNodeGrid
        }

        # position (1,T+1) means not using distance matrix
        if (type == 1){
            d0 <- matrix(0, nTps, nTps)
            dp <- matrix(1, nTps, nTps)
            # d0 <- d0 / s2x + distMask
            # dp <- dp / s2x + distMask
        } else{
            if (dim(ref)[1] == 1){
                dM <-
                    getDistMat(ref, tst[ii, ], refP, tstP[ii, ], biP, weiP, nor)
            } else {
                dM <-
                    getDistMat(ref[ii, ], tst[ii, ], refP[ii, ], tstP[ii, ],
                        biP, weiP, nor)
            }
            d0 <- dM[ , , 1, drop=TRUE]
            dp <- dM[ , , 2, drop=TRUE]
            # d0 <- d0 / s2x + distMask + diagMask
            # dp <- dp / s2x + distMask
        }

        d0ext <- cbind(d0, matrix(0, dim(d0)[1], 1))
        dpext <- cbind(dp, matrix(0, dim(dp)[1], 1))
        # edges from src and to sink
        ssPair[(ssOfst + 1):(ssOfst + nNodeGrid), ] <-
            ssPair[(ssOfst + 1):(ssOfst + nNodeGrid), ] +
            matrix(d0ext[as.vector(ssTmpPos)] + as.vector(ssTmpWt), nNodeGrid,2)
        # edges between nodes
        if (type == 1){
            tmp <- matrix(0, dim(eeTmp)[1], 3)
            tmp[eeTmp[ , 4] == 2, ] <-
                cbind(eeTmp[eeTmp[ , 4] == 2, c(1,2)] + ssOfst,
                    (dpext[eeTmp[eeTmp[ , 4] == 2, 3]] +
                        eeTmp[eeTmp[ , 4] == 2, 5]) * eeTmp[eeTmp[, 4] == 2, 6])
            tmp[eeTmp[ , 4] != 2, ] <-
                cbind(eeTmp[eeTmp[ , 4] != 2, c(1,2)] + ssOfst,
                    (d0ext[eeTmp[eeTmp[ ,4] != 2, 3]] +
                        eeTmp[eeTmp[ , 4] != 2, 5]) * eeTmp[eeTmp[ ,4] != 2, 6])

        } else{
            tmp <- cbind(eeTmp[, c(1,2)] + ssOfst,
                (dpext[eeTmp[, 3]] + eeTmp[, 5]) * eeTmp[, 6])
        }

        if (type == 2){
            eePair[ ,c(1,2)] <- tmp[ ,c(1,2)]
            eePair[ ,3] <- eePair[ ,3] + tmp[ ,3]
        } else{
            eePair[(eeOfst + 1):(eeOfst + nEdgeGrid), c(1,2,3)] <- tmp
        }
    }

    ss <- ssPair
    ee <- eePair

    gtwInfo <-
        list(ss = ss, ee = ee, nodePos = nodePos, weightPos = weightPos,
            weightVal = weightVal, pEdge = pEdge, pEdgeSS = pEdgeSS,
            pEdgeEE = pEdgeEE, dEdge = dEdge, dEdgeInt = dEdgeInt,
            dEdgeIntSS = dEdgeIntSS, dEdgeIntEE = dEdgeIntEE,
            validMap = validMap, srcSinkPos = st01, nTps = nTps,
            nPix = nPix, nNodeGrid = nNodeGrid,
            nNodes = nNodeGrid * nPix + 2, nEdgeGrid = nEdgeGrid,
            ssTmpPos = ssTmpPos)

    if (smoBase == 0 || type == 2)
        return(gtwInfo)

    ## edges for between pairs
    if (type == 0){
        nn <- 0
        for (jj in seq_len(nPix)){
            for (ii in seq(jj + 1, nPix, length = max(0, nPix - jj))){
                if (validMap[ii, jj] != 1)
                    next
                idx1 <- nn + 1
                idx2 <- nn + nNodeGrid
                eePair[(idx1:idx2) + nPix * nEdgeGrid, 1] <- (jj - 1) *
                    nNodeGrid + (seq_len(nNodeGrid))
                eePair[(idx1:idx2) + nPix * nEdgeGrid, 2] <- (ii - 1) *
                    nNodeGrid + (seq_len(nNodeGrid))
                smo0 <- smoMap[ii, jj, 2]
                # smo0 = min(smoMap(ii,jj),smoMap(ii+1,jj))
                eePair[(idx1:idx2) + nPix * nEdgeGrid, 3:4] <- smo0
                nn <- nn + nNodeGrid
            }
        }
        gtwInfo$ee <- eePair
        return(gtwInfo)
    }


    ## type == 1
    ia <- which(!(duplicated(as.vector(dEdgeInt))))
    tdp <- dEdgeInt[ia]
    tdpSort <- sort(tdp, index.return=TRUE)
    tdp <- tdpSort$x
    ia <- ia[tdpSort$ix]

    ll <- dim(dEdgeInt)[1]
    tdP2 <- cbind(dEdge[trunc(ia / ll) * ll + ia],
        dEdge[trunc(ia / ll + 1) * ll + ia])
    GridMap <- matrix(0, 2 * nTps, 2 * nTps)

    for (ii in seq_len(dim(tdP2)[1] - 2)){
        if (round(tdP2[ii, 1] / 0.5) != tdP2[ii, 1] / 0.5){
            if (round(tdP2[ii, 2] / 0.5) != tdP2[ii, 2] / 0.5){
                GridMap[tdP2[ii, 1] * 2 - 0.5, tdP2[ii, 2] * 2 - 0.5] <- ii
            }
        }
    }
    eeSpa <- matrix(0, nTps * nTps * 4 * nNeibPair, 4)
    outIdx <- matrix(TRUE, nTps * nTps * 4 * nNeibPair, 1)

    nn <- 0
    for (jj in seq_len(nPix)){
        for (ii in seq(jj + 1, nPix, length = max(0, nPix - jj))){
            if (validMap[ii, jj] != 1)
                next
            nowIdx <- jj
            tgtIdx <- ii
            idx1 <- nn + 1
            idx2 <- nn + nTps * nTps * 4

            tmpPath <- path[[pairMap[ii, jj]]]
            gtwEdge <- matrix(0, nTps * nTps * 4, 3)
            gtwEdgePP <- matrix(0, nTps * nTps * 4, 2)
            # outIdxt = zeros(nTps*nTps*4, 1)
            gtwEdge[ ,3] <- 1
            for (kk in seq_len(dim(tmpPath)[1])){
                if ((tmpPath[kk, 1] == 0 || tmpPath[kk, 2] == 0))
                    next
                tmpSp <- nTps * 4 * (tmpPath[kk, 2] - 1)

                gtwEdge[(tmpSp + 1):(tmpSp + nTps * 2), c(1, 2)] <-
                    GridMap[ , 2 * tmpPath[kk, c(1, 2)] - 1]
                gtwEdge[(tmpSp + nTps * 2 + 1):(tmpSp + nTps * 4), c(1, 2)] <-
                    GridMap[ , 2 * tmpPath[kk, c(1, 2)]]
                gtwEdgePP[(tmpSp + 1):(tmpSp + nTps * 4), c(1, 2)] <-
                    tmpPath[kk, c(1, 2)]
            }
            outIdxt <- gtwEdge[ , 1] != 0 & gtwEdge[ , 2] != 0

            tempV1 <- log(tst[nowIdx, gtwEdgePP[outIdxt, 1]] + 1)
            tempV2 <- log(tst[tgtIdx, gtwEdgePP[outIdxt, 2]] + 1)
            tempP1 <- (1 - pnorm(tempV1, mu[nowIdx], sigma[nowIdx])) ^ weiP
            tempP2 <- (1 - pnorm(tempV2, mu[tgtIdx], sigma[tgtIdx])) ^ weiP
            gtwEdge[outIdxt, 3] <- min(tempP1, tempP2)

            gtwEdge[gtwEdge[ , 3] == 0, 3] <- 1 / capRev

            outIdx[idx1:idx2] <- outIdxt
            eeSpa[idx1:idx2, 1] <- (nowIdx - 1) * nNodeGrid + gtwEdge[ , 1]
            eeSpa[idx1:idx2, 2] <- (tgtIdx - 1) * nNodeGrid + gtwEdge[ , 2]
            smo0 <- smoMap[ii, jj, 1]
            eeSpa[idx1:idx2, 3] <- smo0 / gtwEdge[ , 3]
            eeSpa[idx1:idx2, 4] <- smo0 / gtwEdge[ , 3]
            #             eeSpa(idx1:idx2,3) = smo0
            #             eeSpa(idx1:idx2,4) = smo0
            nn <- nn + nTps * nTps * 4
        }
    }
    eeSpa <- eeSpa[outIdx, ]
    tmpSort <- sort(eeSpa[ , 2], index.return = TRUE)
    eeSpa <- eeSpa[tmpSort$ix, ]
    ee <- rbind(eePair, eeSpa)

    gtwInfo$ee <- ee
    return(gtwInfo)
}

buildPairTemplate <- function(nTps, win, pmCost) {
    # Get a template of graphs for single pair using coordinates
    # Get both primal edges and dual edges in coordinates
    # Use extra src and sink node to allow partial matching
    # Specify partial matching cost by pmCost of size 1 x winSize
    # weightPos gives pixel dependent portion of the weight based on position
    # in the distance matrix
    # If take position (1,T+1), it means not using the distance matrix
    # weightVal gives prespecified weight independent of pixel, zero means not
    # using Edge types (eType).
    # 1: diagonal/from grid to sink
    # 2: horizontal or vertical
    # 3: from src to grid
    #
    # !! this could introduce bias, need better normalization
    #
    # Inf means maximum value (usually infinite capacity)
    # 0 means there is nothing
    #
    # Pairs in the form (ref,tst). Weight position also (ref,tst)
    # For coordinates, actually this means ref in x axis and tst in y axis

    ## setup

    if (win <= 1)
        stop("winSize should be integer larger than 1\n")

    s0 <- c(0,0)
    t0 <- c(nTps + 1, nTps + 1)
    s1 <- c(nTps + 1, 0)
    t1 <- c(0, nTps + 1)
    st01 <- c(s0, t0, s1, t1)

    # diagnal, sub diagonals, extra edges to s and t, vertical and horizontal
    # edges all in two directions, the reverse is infinite, but we ONLY build
    # finite portion here
    nEdges <- nTps - 1 + (nTps - 2 + nTps - win) * (win - 1) +
        2 * (2 * win - 1) + 2 * (nTps - 1 + nTps - win + 1) * (win - 1)

    # primal and dual edges as well as weights in coordinate form
    # node1 -> node2: (x,y) for node1, (x,y) for node 2
    # (x,y) for weight matrix position, (1,T+1) means none exist position
    pEdge <- matrix(0, nEdges, 4)
    dEdge <- matrix(0, nEdges, 4)
    cPos <- cbind(matrix(1, nEdges, 1), matrix(1, nEdges, 1) * (nTps + 1))
    cVal <- matrix(0, nEdges, 1)
    eType <- matrix(0, nEdges,1)

    ## assign edges
    # loop through all nodes in primal graph, get edges to the right, top or
    # top-right save each edge, and the corresponding dual edge, as well as the
    #  position in the weighting matrix edges outside the window is not included

    # from s to nodes
    pEdge[1, ] <- c(s0,1,1)
    dEdge[1, ] <- c(1.5,0.5,0.5,1.5)
    nn <- 2
    for (w in 2:win){
        pEdge[nn, ] <- c(s0, 1, w)  # close to top left
        pEdge[nn + 1, ] <- c(s0, w, 1)  # close to bottom right
        if (w == win){
            dEdge[nn, ] <- c(0.5, w - 0.5, t1)
            dEdge[nn + 1, ] <- c(s1, w - 0.5, 0.5)
        } else {
            dEdge[nn, ] <- c(0.5, w - 0.5, 0.5, w + 0.5)
            dEdge[nn + 1, ] <- c(w + 0.5, 0.5, w - 0.5, 0.5)
        }
        eType[nn:(nn + 1)] <- 3
        cVal[nn] <- pmCost[w]
        cVal[nn + 1] <- pmCost[w]
        nn <- nn + 2
    }

    # within grid
    for (x0 in seq_len(nTps)){  # ref
        for (y0 in seq_len(nTps)){  # tst
            if (x0 == nTps && y0 == nTps)
                next

            if (y0 < (x0 - win + 1) || y0 > (x0 + win - 1))
                next

            x <- x0
            y <- y0 + 1  # top

            if (y >= (x - win + 1) && y <= (x + win - 1) && x <= nTps &&
                y <= nTps){
                pEdge[nn, ] <- c(x0, y0, x, y)
                cPos[nn, ] <- c(x0, y0)
                if (x0 == 1){
                    dEdge[nn, ] <- c(x0 + 0.25, y0 + 0.75, 0.5, y0 + 0.5)
                } else if (x0 == nTps){
                    dEdge[nn, ] <- c(nTps + 0.5, y0 + 0.5, x0 - 0.25, y0 + 0.25)
                } else{
                    dEdge[nn, ] <- c(x0 + 0.25, y0 + 0.75, x0 - 0.25, y0 + 0.25)
                }
                eType[nn] <- 2
                nn <- nn + 1
            }

            x <- x0 + 1
            y <- y0  # right

            if (y >= (x - win + 1) && y <= (x + win - 1) && x <= nTps &&
                y <= nTps){
                pEdge[nn, ] <- c(x0, y0, x, y)
                cPos[nn, ] <- c(x0, y0)
                if (y0 == 1){
                    dEdge[nn, ] <- c(x0 + 0.5, 0.5, x0 + 0.75, y0 + 0.25)
                } else if (y0 == nTps){
                    dEdge[nn, ] <- c(x0 + 0.25, y0 - 0.25, x0 + 0.5, nTps + 0.5)
                } else{
                    dEdge[nn, ] <- c(x0 + 0.25, y0 - 0.25, x0 + 0.75, y0 + 0.25)
                }
                eType[nn] <- 2
                nn <- nn + 1
            }

            x <- x0 + 1
            y <- y0 + 1   # topright (diagonal)

            if (y >= (x - win + 1) && y <= (x + win - 1) && x0 < nTps &&
                y0 < nTps){
                pEdge[nn, ] <- c(x0, y0, x, y)
                cPos[nn, ] <- c(x0, y0)
                if (y == (x + win - 1)){
                    dEdge[nn, ] <- c(x0 + 0.75, y0 + 0.25, t1)
                } else if (y == (x - win + 1)){
                    dEdge[nn, ] <- c(s1, x0 + 0.25, y0 + 0.75)
                } else{
                    dEdge[nn, ] <- c(x0 + 0.75, y0 + 0.25, x0 + 0.25, y0 + 0.75)
                }
                eType[nn] <- 1
                nn <- nn + 1
            }
        }
    }
    # from nodes to t
    pEdge[nn, ] <- c(nTps, nTps, t0)
    dEdge[nn, ] <- c(nTps + 0.5, nTps - 0.5, nTps - 0.5, nTps + 0.5)
    cPos[nn, ] <- c(nTps, nTps)
    nn <- nn + 1
    for (w in 2:win){
        pEdge[nn, ] <- c(nTps - w + 1, nTps, t0)  # close to top left
        pEdge[nn + 1, ] <- c(nTps, nTps - w + 1, t0)  # close to bottom right
        if (w == win){
            dEdge[nn, ] <- c(nTps + 1 - w + 0.5, nTps + 0.5, t1)
            dEdge[nn + 1, ] <- c(s1, nTps + 0.5, nTps + 1 - w + 0.5)
        } else {
            dEdge[nn, ] <- c(nTps + 1 - w + 0.5, nTps + 0.5, nTps + 1 - w - 0.5,
                nTps + 0.5)
            dEdge[nn + 1, ] <- c(nTps + 0.5, nTps + 1 - w - 0.5, nTps + 0.5,
                nTps + 1 - w + 0.5)
        }
        cPos[nn, ] <- c(nTps - w + 1, nTps)
        cPos[nn + 1, ] <- c(nTps, nTps - w + 1)
        cVal[nn] <- pmCost[w]
        cVal[nn + 1] <- pmCost[w]
        eType[nn : (nn + 1)] <- 1
        nn <- nn + 2
    }

    return(list(pEdge = pEdge, dEdge = dEdge, cPos = cPos, cVal = cVal,
                st01 = st01, eType = eType))
}
