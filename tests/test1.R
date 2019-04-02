samples <-  rbind( c(0,  3,  3,  4,  8, 10),
                   c(0,  8, 14, 18, 24, 27),
                   c(0,  9, 17, 21, 28, 36),
                   c(0, 13, 22, 26, 34, 49),
                   c(0, 14, 23, 30, 38, 58) )

mir = FALSE
strNum = 1
diaNum = 1
biP = TRUE


yyy = buildMultiParaValidmap(samples, mir, strNum, diaNum, biP)
View(yyy$validMap)
yyy$eachMap

validMap <- yyy$validMap

logt = 0
nor = 1
dia = 0
noiseVar = 1
smo = 10^-32
maxStp = 2

param <- initGtwParam(validMap, noiseVar, maxStp, smo, dia, logt, nor)

mu = 0
sigma = 1

edgeInfo <- buildPairTemplate(nTps, win, pmCost)
View(edgeInfo$pEdge)
View(edgeInfo$dEdge)
View(edgeInfo$cPos)
View(edgeInfo$cVal)
View(edgeInfo$st01)
View(edgeInfo$eType)



ll = vector('list', 10)
names(ll)[1] <- 2


c(1,2,3) + (c(5,6,9) - 1) * (10)


ref = 1:10
tst = (3:12) * 2
refP = ref
tstP = tst
biP = TRUE
weiP = 0
nor = 1
dM <- getDistMat(ref, tst, refP, tstP, biP, weiP, nor)
