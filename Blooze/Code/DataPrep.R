
# Author: tim
###############################################################################
setwd("/home/tim/git/Blooze/Blooze")
library(HMDHFDplus)
library(reshape2)
library(LifeTable)

source("Code/ExtrapFuns.R")


Dx  <- readHMDweb("JPN","Deaths_1x1",us,pw)
Ex  <- readHMDweb("JPN","Exposures_1x1",us,pw)

DX <- acast(Dx, Age~Year, value.var = "Female")
EX <- acast(Ex, Age~Year, value.var = "Female")

MX <- suppressWarnings(extrap(DX,EX,130))

qx <- LT(ages = 0:130, Mx = MX[,"2010"],mxsmooth = FALSE, axsmooth = FALSE, axmethod = "midpoint", sex = "female")$qx

save(qx,file = "Data/qxf.Rdata")
#save(MX, file = "Data/JPNMXflong.Rdata")
#matplot(0:130, MX,log = "y",type='l')
#abline(v=90)
