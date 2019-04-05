library(xcms)
library(ncGTW)

CPWmin = 2 # UPLC peakwidths can be as narrow as this - in practice, they are not usually so short but the other parameters ensure noise peaks are not prevalent
CPWmax = 25 # widest peak appears to be about 20s but I generally add some time to this as peaks are not easy to deconvolve by eye
ppm = 15 # calculated from the raw data - see spreadsheet in zipfile
xsnthresh = 3 # generous, but minimum acceptable S/N
#snth = 10
LM=F # in practice, including LM files does not appear to really influence final result - also peakpicking does not always appear to work when included
integrate=2 # raw data of sufficient quality that it can be used directly. Would consider waves in presence of strong peak asymmetry/skew, but generally unecessary
RTerror=6 # check drift of peak tops from chromatograms - see spreadsheet in zipfile
MZerror=0.05
prefilter=c(8, 1000)
mzdiff = -0.001


temp_path <- "J:/mesa/iQC335"
tempFiles <- list.files(temp_path, pattern = "mzXML", recursive = F, full.names = T)



tempInd <- matrix(0, length(tempFiles), 1)
for (n in 1:length(tempFiles))
{
  tempCha <- tempFiles[n]
  tempLen <- nchar(tempCha)
  tempInd[n] <- as.numeric(substr(tempCha, regexpr("_iQC", tempCha) + 4, tempLen - 6))
}

tempFiles <- tempFiles[sort.int(tempInd, index.return = TRUE)$ix]


bpParam <- SnowParam(workers = 14)

print(time1 <- Sys.time())
ds <- xcmsSet(tempFiles, method="centWave", peakwidth=c(CPWmin,CPWmax), ppm=ppm,
              noise=xsnthresh, integrate=integrate,
              prefilter=prefilter, mzdiff = mzdiff, BPPARAM = bpParam)
print(time2 <- Sys.time())
print((time2 - time1))

gds <- group(ds, method="density", mzwid=MZerror, bw=RTerror, minfrac = 0.25)
agds <- retcor(gds, method="peakgroups", family="s", missing = 5, extra=5)
bw = 7 # adjust
Sys.time()
xcmsLargeWin <- group(agds, method="density", mzwid=MZerror, bw=bw, minfrac = 0.5)
Sys.time()
xcmsSmallWin <- group(agds, method="density", mzwid=MZerror, bw=0.3, minfrac = 0.00001, minsamp = 2)
Sys.time()

excluGroups <- misalignDetect(xcmsLargeWin, xcmsSmallWin, ppm)
filePaths <- xcmsLargeWin@filepaths
time3 <- Sys.time()
ncGTWinputs <- loadProfile(filePaths, excluGroups, BPPARAM = bpParam)
time4 <- Sys.time()
print(time4 - time3)

for (n in 1:length(ncGTWinputs)){
  plotGroup(ncGTWinputs[[n]], xcmsLargeWin@rt$corrected)
  Sys.sleep(5)
}

allTime <- vector('list', dim(excluGroups)[1])
ncGTWoutputs <- vector('list', dim(excluGroups)[1])

cat(format(Sys.time()), 'ncGTW starts !\n')
for (n in 1:dim(excluGroups)[1]){
  ncGTWoutputs[[n]] <- ncGTW(ncGTWinputs[[n]], xcmsLargeWin, 10)
  allTime[[n]] <- Sys.time()
  cat(format(Sys.time()), 'Group', ncGTWinputs[[n]]$groupInfo['index'], 'is realigned.\n')
}

ncGTWres <- xcmsLargeWin
ncGTWRt <- vector('list', dim(excluGroups)[1])
for (n in 1:dim(excluGroups)[1]){
  adjustRes <- adjustRT(ncGTWres, ncGTWinputs[[n]], ncGTWoutputs[[n]], ppm)
  ncGTWres@peaks <- adjustRes$peaks
  ncGTWRt[[n]] <- adjustRes$rtncGTW
  cat(format(Sys.time()), 'Group', ncGTWinputs[[n]]$groupInfo['index'], 'is adjusted.\n')
}

ncGTWres@groups <- excluGroups[ , 2:9]
ncGTWres@groupidx <- xcmsLargeWin@groupidx[excluGroups[ , 1]]
ncGTWres@rt$corrected <- vector('list', length(ncGTWres@filepaths))
for (n in 1:length(ncGTWres@filepaths)){
  ncGTWres@rt$corrected[[n]] <- vector('list', dim(excluGroups)[1])
  for (m in 1:dim(ncGTWres@groups)[1])
    ncGTWres@rt$corrected[[n]][[m]] <- ncGTWRt[[m]][[n]]
}


require(xcms)
require(ncGTW)

assignInNamespace("fillPeaksChromPar", fillPeaksChromPar, ns="xcms",
                  envir = as.environment("package:xcms"))

Sys.time()
ncGTWresFilled <- fillPeaks(ncGTWres, BPPARAM = SnowParam(workers = 8))
Sys.time()
View(ncGTWresFilled@peaks)
n = 59

ttt= groupidx(ncGTWres)

dim(ncGTWres@rt$corrected[[excluGroups[n, 1]]])
y <- matrix(0, dim(ncGTWinputs[[n]]$rtRaw)[1], dim(ncGTWinputs[[n]]$rtRaw)[2])
for (m in 1:dim(y)[1])
  y[m, ] <- ncGTWRt[[n]][[m]][ncGTWinputs[[n]]$rtRaw[m, ]]

matplot(t(y), t(ncGTWinputs[[n]]$profiles), type = 'l')

Sys.time()
tt = fillPeaks(xcmsLargeWin, BPPARAM = SnowParam(workers = 8))
Sys.time()



dataSub <- parSpec[1:parNum[n], , n]
sampleInd <- parInd[1:parNum[n], n]

parInfo = parInfos[[7]]
dataSub = parspec[1:parnum, ]
sampleInd = parind[1:parnum]

n = 2
ncGTWinput <- ncGTWinputs[[n]]
ncGTWoutput <- ncGTWoutputs[[n]]

for (n in 1:dim(excluGroups)[1]){
  ncGTWoutputs[[n]]$downSample = 2
}

