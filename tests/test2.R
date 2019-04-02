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
validMap <- gtwPrep$validMap


param <- initGtwParam(validMap, noiseVar, maxStp, 10^-32, dia, logt, nor)
gtwInfo <- buildGTWgraph(ref, tst, validMap, param, mu, sigma, biP, weiP)
eeBet <- gtwInfo$ee[(gtwInfo$nEdgeGrid * gtwInfo$nPix + 1) : dim(gtwInfo$ee)[1], ]

param <- initGtwParam(validMap, noiseVar, maxStp, 0, dia, logt, nor)
gtwInfo0 <- buildGTWgraph(ref, tst, validMap, param, mu, sigma, biP, weiP)
test <- graphCut(gtwInfo0$ss, gtwInfo0$ee)
#s_temp_num = sum((labels0(ee_s(:,1)) + labels0(ee_s(:,2)))==1);


#param = gtw.initGtwParam(validMap,noiseVar,maxStp, 100000000000, dia, logt, nor);
#[ ss,ee,~, ~, ~ ] = gtw.buildGraph4Aosokin_GraphMerge( ref, tst, validMap, param, mu, sigma, biP, weiP);
#[cost_max, ~] = aoIBFS.graphCutMex(ss,ee);
#s_temp = (cost_max - cost_min)/s_temp_num;
