#' Edited XCMS fillPeaksChromPar for feature-wise warping functions
#'
#' This function is edited from \code{fillPeaksChromPar} in
#' \code{\link[xcms]{fillPeaks.chrom-methods}} to accept feature-wise warping
#' functions.
#' @param arg A list sent from \code{\link[xcms]{fillPeaks.chrom-methods}}.
#' @details This function is for parallelly filling missing peaks with
#' feature-wise warping functions. The original function in
#' \code{\link[xcms]{fillPeaks.chrom-methods}} can only handle sample-wise
#' warping functions.
#' @return A list of sameple index vector and filled peak matrix.

fillPeaksChromPar <- function(arg) {

    suppressMessages(requireNamespace("xcms", quietly=TRUE))
    suppressMessages(requireNamespace("ncGTW", quietly=TRUE))
    message('This is edited from XCMS (fillPeaksChromPar).....')

    params <- arg$params
    myID <- arg$id
    message(arg$file)

    prof <- params$prof
    rtcor <- params$rtcor
    peakrange <- params$peakrange
    expand.mz <- params$expand.mz
    expand.rt <- params$expand.rt
    gvals <- params$gvals$gvals

    lcraw <- xcmsRaw(arg$file, profmethod=params$prof$method, profstep=0)

    if (length(params$dataCorrection) > 1) {
        ## Note: dataCorrection (as set in the xcmsSet function) is either
        ## 1 for all or for none.
        if (any(params$dataCorrection == 1))
            lcraw <- stitch(lcraw, AutoLockMass(lcraw))
    }

    if (exists("params$polarity") && length(params$polarity) > 0) {
        if (length(params$polarity) > 0) {
            ## Retain wanted polarity only
            lcraws <- split(lcraw, lcraw@polarity, DROP=TRUE)
            lcraw <- lcraws[[params$polarity]]
        }
    }


    if (!is.list(rtcor))
        rtcor <- list(rtcor)
    lcraw@profparam <- rtcor

    ## Expanding the peakrange
    incrMz <- (peakrange[, "mzmax"] - peakrange[, "mzmin"]) / 2 * (expand.mz -1)
    peakrange[, "mzmax"] <- peakrange[, "mzmax"] + incrMz
    peakrange[, "mzmin"] <- peakrange[, "mzmin"] - incrMz
    incrRt <- (peakrange[, "rtmax"] - peakrange[, "rtmin"]) / 2 * (expand.rt -1)
    peakrange[, "rtmax"] <- peakrange[, "rtmax"] + incrRt
    peakrange[, "rtmin"] <- peakrange[, "rtmin"] - incrRt

    naidx <- which(is.na(gvals[,myID]))

    newpeaks <- getPeaksncGTW(lcraw, peakrange[naidx, , drop=FALSE],
        step=prof$step, naidx)



    list(myID=myID, newpeaks=cbind(newpeaks, sample=myID))
}

#' Edited XCMS getPeaks for feature-wise warping functions
#'
#' This function is edited from \code{\link[xcms]{getPeaks-methods}} to accept
#' feature-wise warping functions.
#' @param object An \code{\link[xcms]{xcmsRaw-class}} object.
#' @param peakrange \code{matrix} with 4 required columns \code{"mzmin"},
#'   \code{"mzmax"}, \code{"rtmin"} and \code{"rtmax"}.
#' @param step \code{numeric(1)} defining the bin size for the profile matrix
#'   generation.
#' @param naidx A vector contains the sample indexes need to be filled.
#' @details This function is for parallelly filling missing peaks with
#' feature-wise warping functions. The original code function in
#' \code{\link[xcms]{getPeaks-methods}} can only handle sample-wise warping
#' functions.
#' @return A list of sameple index vector and filled peak matrix.

getPeaksncGTW <- function(object, peakrange, step=0.1, naidx) {
    suppressMessages(requireNamespace("xcms", quietly=TRUE))
    message('ncGTW fillpeaks')
    ## Here we're avoiding the profFun call.
    if (all(c("mzmin","mzmax","rtmin","rtmax") %in% colnames(peakrange)))
        peakrange <- peakrange[,c("mzmin","mzmax","rtmin","rtmax"), drop=FALSE]
    stime <- object@profparam

    pi <- profinfo(object)
    method <- pi$method
    if (missing(step))
        step <- pi$step
    if (step == 0)
        step <- 0.1
    baselevel <- pi$baselevel
    basespace <- pi$basespace
    vps <- diff(c(object@scanindex, length(object@env$mz)))

    message("method_new: ", method, " ")
    message("step: ", step)
    ## Create the profile matrix:
    pMat <- xcms:::.createProfileMatrix(mz = object@env$mz,
                                        int = object@env$intensity,
                                        valsPerSpect = vps,
                                        method = method,
                                        step = step,
                                        baselevel = baselevel,
                                        basespace = basespace,
                                        returnBreaks = TRUE,
                                        baseValue = 0,
                                        mzrange. = NULL)
    brks <- pMat$breaks
    pMat <- pMat$profMat  ## rows are masses, cols are retention times/scans.
    bin_size <- diff(brks[c(1,2)])
    bin_half <- bin_size / 2
    ## Calculate the mean mass per bin using the breaks used for the binning.
    ## Note: these define the real mass breaks as they have been used for the
    ## binning. Simply using seq(floor...) as in the original code is wrong
    ## because the mass bins are calculated wrongly. The bin size is != step,
    ## bin size is marginally smaller and, for larger mz the correct mass
    ## bin will be wrongly identified.
    mass <- brks[-length(brks)] + bin_half ## midpoint for the breaks
    mass_range <- range(mass)

    ## Prepare the result matrix.
    cnames <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into", "maxo")
    rmat <- matrix(nrow = nrow(peakrange), ncol = length(cnames))
    colnames(rmat) <- cnames

    for (i in order(peakrange[, 1])) {

        if (length(stime) == 1) {
            stime_temp <- stime[[1]]
        } else {
            message("Feature-wise warpping function:", naidx[i])
            stime_temp <- stime[[naidx[i]]]
        }

        imz <- xcms:::findRange(mass, c(peakrange[i, 1] - bin_half,
                                        peakrange[i, 2] + bin_half), TRUE)
        iret <- xcms:::findRange(stime_temp, peakrange[i, 3:4], TRUE)
        idx_imz <- imz[1]:imz[2]
        idx_iret <- iret[1]:iret[2]
        ## Extract the intensity matrix for the mz-rt range: rows are mz, cols
        ## rt values.
        ymat <- pMat[idx_imz, idx_iret, drop = FALSE]
        ## Define the maximum intensity, is one value per mz.
        ymax <- xcms:::colMax(ymat)
        iymax <- which.max(ymax)

        ## The width in rt.
        pwid <- diff(stime_temp[iret])/diff(iret)

        ## Calculate sum across rt. For each mz we get one value.
        rosm <- rowSums(ymat)
        limz <- length(idx_imz)
        if (length(rosm) != limz) { ## that happens for some reason
            warning("weighted.mean  : x and w must have the same length \n")
            rosm <- rep(1, limz)  ## fallback to mean
        }
        ## mean mz:
        rmat[i, 1] <- weighted.mean(mass[idx_imz], rosm) ## mz; its not the
        ## position of the largest intensity!
        if (is.nan(rmat[i,1]) || is.na(rmat[i,1])) ##  R2.11 :  weighted.mean()
            ## results in NA (not NaN) for zero weights
            rmat[i, 1] <- mean(peakrange[i, c(1,2)])

        rmat[i, 2:3] <- peakrange[i, c(1,2)]            ## mzmin, mzmax
        rmat[i, 4] <- stime_temp[idx_iret][iymax] ## rt
        rmat[i, 5:6] <- peakrange[i, c(3,4)]            ## rtmin, rtmax

        if (peakrange[i, 3] <  stime_temp[1] ||
            peakrange[i, 4] > stime_temp[length(stime_temp)] ||
            is.nan(pwid)) {
            warning("getPeaks: Peak  m/z:", peakrange[i, 1], "-",
                    peakrange[i, 2], ",  RT:", peakrange[i, 3], "-",
                    peakrange[i, 4], "is out of retention time range for ",
                    "this sample (", object@filepath,
                    "), using zero intensity value.\n")
            rmat[i, 7:8] <- 0
        } else {
            rmat[i, 7] <- pwid * sum(ymax)  ## into
            rmat[i, 8] <- ymax[iymax]       ## maxo
        }
    }
    invisible(rmat)
}
