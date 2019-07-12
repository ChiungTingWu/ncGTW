#' Load sample profiles for each peak group
#'
#' This function loads each raw sample profiles from the file with certain m/z
#' and RT range.
#' @param filePaths The character vector of the loading file paths.
#' @param excluGroups The output matrix of \code{\link{misalignDetect}} or
#'   \code{\link[xcms]{xcmsSet-class}}$\code{group}, in which \code{mzmin},
#'   \code{mzmax}, \code{rtmin}, and \code{rtmax} are set as the m/z and RT
#'   range for loading.
#' @param mzAdd The extra m/z range for loading (both sides), and the default is
#'   0.005.
#' @param rtAdd The extra RT range for loading (both sides), and the default is
#'   10 (seconds).
#' @param profstep The size of each m/z bin for peak integration, and the
#'   default is 0.
#' @param BPPARAM A object of \pkg{BiocParallel} to control parallel processing,
#'   and can be created by
#'   \code{\link[BiocParallel:SerialParam-class]{SerialParam}},
#'   \code{\link[BiocParallel:MulticoreParam-class]{MulticoreParam}}, or
#'   \code{\link[BiocParallel:SnowParam-class]{SnowParam}}.
#'
#' @details This function obtains the extracted ion chromatogram for each sample
#' at the givin m/z and RT range with a certain m/z bin size for integration.
#' Considering there may be missing peak by peak detection, \code{mzAdd} and
#' \code{rtAdd} are to increase the integration range.
#' @return A list of the same length as the row number of \code{excluGroups}, in
#' which each element is a \code{\link{ncGTWinput}} object.
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
#' @export

loadProfile <- function(filePaths, excluGroups, mzAdd=0.005, rtAdd=10,
    profstep=0, BPPARAM=BiocParallel::SnowParam(workers=1)){

    if (!is.character(filePaths))
        stop('the \'filePaths\' argument should be a character vector')
    if (!is.numeric(excluGroups))
        stop('the \'excluGroups\' argument should be a numeric matrix')
    cName <- c("mzmed", "mzmin", "mzmax","rtmed", "rtmin", "rtmax")
    if (!all(cName %in% colnames(excluGroups)))
        stop('the \'excluGroups\' argument should contain column names:
             "mzmed", "mzmin", "mzmax","rtmed", "rtmin", "rtmax"')


    groupNum <- dim(excluGroups)[1]
    fileNum <- length(filePaths)

    mzRange <- excluGroups[,c('mzmin', 'mzmax'), drop=FALSE]
    mzRange[,'mzmin'] <- mzRange[,'mzmin'] - mzAdd
    mzRange[,'mzmax'] <- mzRange[,'mzmax'] + mzAdd
    rtRange <- excluGroups[,c('rtmin', 'rtmax'), drop=FALSE]
    rtRange[,'rtmin'] <- rtRange[,'rtmin'] - rtAdd
    rtRange[,'rtmax'] <- rtRange[,'rtmax'] + rtAdd


    timeS <- Sys.time()
    eicList <-
        suppressWarnings(
            BiocParallel::bplapply(filePaths, loadEic, mzRange,
            rtRange, profstep, BPPARAM=BPPARAM))
    timeE <- Sys.time()
    timeDif <- timeE - timeS
    message('The loading time is ',round(timeDif,2), ' ', units(timeDif), '.')
    message('Starting to check the length of each profile...')

    groupProf <- vector("list", groupNum)
    groupRt <- vector("list", groupNum)
    if ('index' %in% colnames(excluGroups)){
        groupInd <- excluGroups[, 'index']
    } else{
        groupInd <- seq_len(nrow(excluGroups))
    }
    names(groupProf) <- paste0('group', groupInd)
    names(groupRt) <- names(groupProf)

    maxMinArr <- array(0, c(groupNum, fileNum, 2))
    maxMinMat <- matrix(0, groupNum, 2)
    maxMinMat <- matrix(0, groupNum, 2)
    maxMinInd <- array(0, c(groupNum, fileNum, 2))

    LenMat <- matrix(0, groupNum, fileNum)
    LenVec <- matrix(0, groupNum, 1)

    ncGTWinputs <- vector('list', groupNum)
    for (i in seq_len(groupNum)){
        for (j in seq_len(fileNum)){
            maxMinArr[i, j, 1] <- min(eicList[[j]][[i]][,1])
            maxMinArr[i, j, 2] <- max(eicList[[j]][[i]][,1])
        }
        maxMinMat[i, 1] <- max(maxMinArr[i, , 1, drop=FALSE])
        maxMinMat[i, 2] <- min(maxMinArr[i, , 2, drop=FALSE])

        for (j in seq_len(fileNum)){
            maxMinInd[i, j, 1] <-
                which.min(abs(eicList[[j]][[i]][, 1] - maxMinMat[i, 1]))
            maxMinInd[i, j, 2] <-
                which.min(abs(eicList[[j]][[i]][, 1] - maxMinMat[i, 2]))
            LenMat[i, j] <- maxMinInd[i, j, 2] - maxMinInd[i, j, 1] + 1
        }

        if (length(unique(LenMat[i,])) != 1){
            warning('The RT resolution may not be same for some samples...')
            LenVec[i] <- min(LenMat[i,])
        } else{
            LenVec[i] <- unique(LenMat[i,])
        }

        groupProf[[i]] <- matrix(0, fileNum, LenVec[i])
        groupRt[[i]] <- matrix(0, fileNum, LenVec[i])
        for (j in seq_len(fileNum)){
            rtProf <- eicList[[j]][[i]][maxMinInd[i, j, 1]:
                                            (maxMinInd[i, j, 1] +
                                                 LenVec[i] - 1), ]
            groupRt[[i]][j, ] <- rtProf[ , 1]
            groupProf[[i]][j, ] <- rtProf[ , 2]
        }

        groupInfo <- excluGroups[i, ]
        if (!'index' %in% names(groupInfo))
            groupInfo <- c(index=i, groupInfo)
        ncGTWinputs[[i]] <- new("ncGTWinput", groupInfo=groupInfo,
                                profiles=groupProf[[i]], rtRaw=groupRt[[i]])

    }
    return(ncGTWinputs)
}



loadEic <- function(filePath, mzRange, rtRange, profstep){
    # suppressMessages(requireNamespace("xcms", quietly = TRUE))
    # rawProf <- suppressMessages(xcmsRaw(filePath, profstep))
    # rawEic <- getEIC(rawProf, mzRange, rtRange)
    rawProf <- suppressMessages(xcmsRaw(filePath, profstep))
    rawEic <- vector('list', dim(mzRange)[1])

    for (n in seq_len(dim(mzRange)[1])){
        tempEic <- xcms::rawEIC(rawProf, mzRange[n, ], rtRange[n, ])
        rawEic[[n]] <- cbind(tempEic$scan, tempEic$intensity)
    }
    rm(rawProf)
    gc()

    # return(rawEic@eic$xcmsRaw)
    return(rawEic)
}
