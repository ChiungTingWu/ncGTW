samples <-  rbind( c(0,  3,  3,  4,  8, 10),
                   c(0,  8, 14, 18, 24, 27),
                   c(0,  9, 17, 21, 28, 36),
                   c(0, 13, 22, 26, 34, 49),
                   c(0, 14, 23, 30, 38, 58) )

mir = TRUE
strNum = 1
diaNum = 1
biP = TRUE

logt = 0
nor = 1
dia = 0
noiseVar = 1
smo = 10^-32
maxStp = 2

mu = 0
sigma = 1
weiP = 0

gtwPrep <- buildMultiParaValidmap(samples, mir, strNum, diaNum, biP)
ref <- gtwPrep$ref
tst <- gtwPrep$tst
validMap <- gtwSamples$validMap

