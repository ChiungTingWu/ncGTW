#' @useDynLib ncGTW, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

data <-     rbind( c(0,  3,  3,  4,  8, 10),
                   c(8,  8, 14, 18, 24, 27),
                   c(2,  9, 17, 21, 28, 36),
                   c(3, 13, 22, 26, 34, 49),
                   c(5, 14, 23, 30, 38, 58),
                   c(1, 31,  9, 23, 11, 10),
                   c(7, 99, 55, 12, 33, 43),
                   c(4, 52, 73,  4, 66, 38),
                   c(9, 74, 21, 12, 44, 55))



filepath <- "MESA_COMBI_BIO_P2_LP_iQC171.mzXML"



line <- grep("fileName",txt,fixed=T) ## find the line that contains the file name
if (length(line) >0) { ## found it
  res <- strsplit(txt[line], ...)  ## extract the result using strsplit etc.
}

Sys.time()
test <- openMSfile(filepath)
Sys.time()
head <- header(test)


Sys.time()
msd <- openMSfile(filepath)

Sys.time()



  pks <- mzR::peaks(msd, 1:2626)

  Sys.time()
  ## Fix problem with single spectrum files (issue #66)
  if (is(pks, "matrix"))
    pks <- list(pks)
  valsPerSpect <- lengths(pks) / 2
  pks <- do.call(rbind, pks)
  hdr_ms1 <- hdr[idx_ms1, header_cols,
                 drop = FALSE]








ttttt = xcmsrawtest@env$profile


filepath <- dir(pattern = "mzXML")
# filepath <- list.files(cdfpath, recursive = FALSE, full.names = TRUE)

# filepath <- "C:/Users/cbil/Google Drive/5th project/test dataset2"

mzrange <- matrix(0, 1, 2)
mzrange[1, 1] <- 269.2253
mzrange[1, 2] <- 269.2274

rtrange <- matrix(0, 1, 2)
rtrange[1, 1] <- 322.1551
rtrange[1, 2] <- 355.3655

data <- vector('list', length(filepath))

Sys.time()
for (n in 1:length(data)){
  xcmsrawtest <- xcmsRaw(filepath[n], 0.01)
  show(Sys.time())
  raweic <- getEIC(xcmsrawtest, mzrange, rtrange)
  show(Sys.time())
  data[[n]] <- raweic@eic$xcmsRaw[[1]]
  show(Sys.time())
}

dataM <- matrix(0, 6, 101)
for (n in 1:length(data)){
  dataM[n, ] <- data[[n]][,2]
}



Sys.time()
seteic <- getEIC(xcmssettest, mzrange)
Sys.time()

Sys.time()
