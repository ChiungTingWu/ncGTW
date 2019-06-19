context("Valid ncGTW")
library(xcms)

test_that("misalignDetect detects misaligned features", {
    data('xcmsExamples')
    xcmsLargeWin <- xcmsExamples$xcmsLargeWin
    xcmsSmallWin <- xcmsExamples$xcmsSmallWin
    ppm <- xcmsExamples$ppm

    # detect misaligned features
    expect_equal(nrow(misalignDetect(xcmsLargeWin, xcmsSmallWin, ppm)), 2)
})

test_that("loadProfile loads all sample profiles", {
    data('xcmsExamples')
    xcmsLargeWin <- xcmsExamples$xcmsLargeWin
    xcmsSmallWin <- xcmsExamples$xcmsSmallWin
    ppm <- xcmsExamples$ppm

    # detect misaligned features
    excluGroups <- misalignDetect(xcmsLargeWin, xcmsSmallWin, ppm)

    # obtain the paths of the sample files
    filepath <- system.file("extdata", package = "ncGTW")
    file <- list.files(filepath, pattern="mzxml", full.names=TRUE)

    tempInd <- matrix(0, length(file), 1)
    for (n in seq_along(file)){
        tempCha <- file[n]
        tempLen <- nchar(tempCha)
        tempInd[n] <- as.numeric(substr(tempCha, regexpr("example", tempCha) + 7,
                                        tempLen - 6))
    }
    # sort the paths by data acquisition order
    file <- file[sort.int(tempInd, index.return = TRUE)$ix]

    # load the sample profiles
    ncGTWinputs <- loadProfile(file, excluGroups)

    expect_equal(length(ncGTWinputs), 2)
    expect_equal(nrow(ncGTWinputs[[1]]@profiles), length(file))
    expect_equal(nrow(ncGTWinputs[[2]]@profiles), length(file))
})

test_that("ncGTWalign realign the misaligned feature, and adjustRT produces the
          new warping functions", {
    data('xcmsExamples')
    xcmsLargeWin <- xcmsExamples$xcmsLargeWin
    xcmsSmallWin <- xcmsExamples$xcmsSmallWin
    ppm <- xcmsExamples$ppm

    # detect misaligned features
    excluGroups <- misalignDetect(xcmsLargeWin, xcmsSmallWin, ppm)

    # obtain the paths of the sample files
    filepath <- system.file("extdata", package = "ncGTW")
    file <- list.files(filepath, pattern="mzxml", full.names=TRUE)

    tempInd <- matrix(0, length(file), 1)
    for (n in seq_along(file)){
        tempCha <- file[n]
        tempLen <- nchar(tempCha)
        tempInd[n] <- as.numeric(substr(tempCha, regexpr("example", tempCha) + 7,
                                        tempLen - 6))
    }
    # sort the paths by data acquisition order
    file <- file[sort.int(tempInd, index.return = TRUE)$ix]

    # load the sample profiles
    ncGTWinputs <- loadProfile(file, excluGroups)

    # initialize the parameters of ncGTW alignment with default
    ncGTWparam <- initncGTWparam()

    # run ncGTW alignment
    ncGTWoutputs <- ncGTWalign(ncGTWinputs[[1]], xcmsLargeWin, 5,
                               ncGTWparam = ncGTWparam)
    expect_equal(length(ncGTWoutputs@ncGTWpath), length(file))

    # adjust RT with the realignment results from ncGTW
    ncGTWres <- xcmsLargeWin

    adjustRes <- adjustRT(ncGTWres, ncGTWinputs[[1]], ncGTWoutputs, ppm)

    expect_equal(nrow(adjustRes@ncGTWpeaks), nrow(peaks(xcmsLargeWin)))
    expect_equal(length(adjustRes@rtncGTW), length(file))
})


