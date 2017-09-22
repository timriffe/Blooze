
# TR: make some R code to do the primary simulations.
setwd("/home/tim/git/Blooze/Blooze")
library(LifeTable)
repetitions <- 15000
burnin      <- 150
obsperiod   <- 40


qx <- local(get(load("Data/qxf.Rdata")))


SimProp <- function(qx, radix = 1e5, burnin = length(qx)+1, obsperiod = 40){
	Na  <- length(qx)
	Ny  <- burnin + obsperiod
	Age <- 1:Na - 1
	Pop <- matrix(0, nrow = Na, ncol = burnin + obsperiod, dimnames = list(Age, 1:Ny-1))
#	sum(Start[Age >= 100]) / sum(Start)
	Pop[, 1] <- qx2lx(qx) * radix
	Pop[1, ] <- radix
	for (i in 2:ncol(Pop)){
		Dxlamb      <- Pop[-Na, i - 1] * qx[-Na]
		Dx          <- rpois(Na - 1, Dxlamb)
		ind         <- Dx >  Pop[-Na, i - 1]
		if (any(ind)){
			Dx[ind]     <- Pop[-Na, i - 1][ind]
		}
		Pop[2:Na,i] <- Pop[-Na, i - 1] - Dx
		#Pop[Pop[,i] < 0,i] <- 0
	}
	Pop <- Pop[,-c(1:burnin)]
	colSums(Pop[Age >= 100, ]) /
	colSums(Pop)
}

highwater <- function(props, n = 5){
	max(as.vector(stats::filter(props, rep(1 / n, n), sides = 2)),na.rm=TRUE)
}
lowwater <- function(props, n = 5){
	min(as.vector(stats::filter(props, rep(1 / n, n), sides = 2)),na.rm=TRUE)
}
highlowwater <- function(props, n = 5){
	mavg <- as.vector(stats::filter(props, rep(1 / n, n), sides = 2))
	range(mavg,na.rm=TRUE)
}
#props <- SimProp(qx)
#plot(props)
#abline(h = asymptote)
#abline(h = highwater(props))
lx <- qx2lx(qx) 
age <- 0:130
asymptote <- sum(lx[age >= 100]) / sum(lx)

repetitions <- 15000
probs <- c(0,.01,.025,.05,.95,.975,.99,1)


# I herewish declare that the plural of radix shall not be
# radices, but rather radishes.
radishes <- round(exp(log(10)*seq(log10(100),log10(1e6),by=0.25) ))

Maxsq <- matrix(0,nrow = length(probs),ncol = length(radishes), dimnames=list(probs,radishes))

# takes a good while to calculate
for (i in 1:length(radishes)){
	cat(radishes[i],"\n")
    simsradix <- replicate(repetitions, 
			highwater(
					SimProp(qx=qx,radix=radishes[i])
	                     )
		             )
	maxs       <- quantile(simsradix, probs = probs)
	Maxsq[, i] <- maxs
}

save(Maxsq, file = "Data/JPNf2010qs.Rdata")

pdf("Figures/SimJPNf2010.pdf")
par(mai=c(1,2,.1,.1))
matplot(radishes, t(Maxsq), type = 'l', lty = 1, col = gray(.2), log = "x",
		xlab = "radix size", ylab = "",xlim = c(57, 1e6), las = 1)
abline(h = asymptote, col = "red")
text(100, Maxsq[, 1], probs, cex = .8, pos = 2)
text(1e6, asymptote * 1.1, "Asympotote", pos = 2)
segments(1e6 * .9, asymptote * 1.1, 1.2e6, asymptote * 1.01, col = "red")
text(2, mean(Maxsq[, 1]), "highest prop\ncentenarian\nfrom 5yr\nmoving avg\nin 40 yrs",xpd=TRUE)
text(240,.0039, "Quantile contours\nfrom 15000 simulations\nper radix",pos=4)
dev.off()