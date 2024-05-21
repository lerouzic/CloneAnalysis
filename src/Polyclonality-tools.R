
############################## Constants ###############################
max.nbclones <- 7

col.names <-  c("G","R","Y","B")
colors    <-  c(G="green", R="darkred", Y="yellow", B="turquoise")
obs.names <-  apply(expand.grid(lapply(col.names, function(x) c("", x))), 1, paste, collapse="") 

######################## Frequency manipulation ########################
#
# conventions:
# p  : color frequencies (vector of size 4, sum = 1).      Names = col.names
# f  : observation frequencies (vector of size 16, sum=1). Names = obs.names
# obs: observation counts (vector of size 16).             Names = obs.names

logit    <- function(p) log(p / (1-p))
invlogit <- function(u) 1/(1+exp(-u))

logitp2p <- function(lpG, lpR, lpY) {
	pG <- invlogit(lpG)
	pR <- (1-pG)*invlogit(lpR)
	pY <- (1-pG)*(1-pR)*invlogit(lpY)
	pB <- 1-pG - pR - pY
	c(G=unname(pG), R=unname(pR), Y=unname(pY), B=unname(pB))
}

obs2p <- function(obs) { 
	# By a maximum likelihood procedure, assuming independent observations
	stopifnot(length(obs) == length(obs.names))
	mll <- function(lpG, lpR, lpY) {
		f <- setNames(rep(1, length(obs)), nm=obs.names)
		p <- logitp2p(lpG, lpR, lpY)
		f[which(names(f) == "")] <- 0 # Assuming observations without colors are absent
		f <- f * ifelse(grepl("G", names(f)), p["G"], 1-p["G"])
		f <- f * ifelse(grepl("R", names(f)), p["R"], 1-p["R"])
		f <- f * ifelse(grepl("Y", names(f)), p["Y"], 1-p["Y"])
		f <- f * ifelse(grepl("B", names(f)), p["B"], 1-p["B"])
		
		# Fixing numerical errors
		f[f < 0] <- 0
		f <- f/sum(f)
		
		-dmultinom(obs, prob=f, log=TRUE)
	}
	mm <- try(stats4::mle(mll, start=c(lpG=0, lpR=0, lpY=0)))
	if (class(mm) == "try-error") return(NA)
	list(lp.estimates = mm@coef, p.estimates=logitp2p(mm@coef[1], mm@coef[2], mm@coef[3]))
}

############################ Model fitting #############################

cache.gridcat <- lapply(seq_len(max.nbclones), function(cl) apply(expand.grid(replicate(cl, col.names, simplify=FALSE)), 1, function(x) paste0(unique(x[order(match(x, col.names))]), collapse="")))

# Nomenclature:
#
# exact  : exactly nbclones clones
# random : average meanclones clones (Poisson distribution)

freq.model.exact <- function(p, nbclones) {
	ref.nn <- cache.gridcat[[nbclones]]
	prbs <- setNames(apply(expand.grid(replicate(nbclones, p, simplify=FALSE)), 1, prod), nm=ref.nn)
	sapply(obs.names, function(nn) sum(prbs[ref.nn == nn]))
}

freq.model.random <- function(p, meanclones) {
	distclone <- dpois(seq_len(max.nbclones), lambda=meanclones)
	distclone <- distclone / sum(distclone)
	colSums(do.call(rbind, lapply(seq_len(max.nbclones), freq.model.exact, p=p))*distclone)
}

mllik.exact <- function(x, obs, p)  
	-dmultinom(obs, prob=freq.model.exact(p, nbclones=round(x)), log=TRUE)

mllik.random <- function(x, obs, p) # x can be a vector
	-dmultinom(obs, prob=freq.model.random(p, meanclones=x), log=TRUE)

maxlik.exact.manual <- function(obs, p=obs2p(obs)$p.estimates, CI=0.95) {
	nbcl <- seq_len(max.nbclones)
	mll  <- sapply(nbcl, mllik.exact, obs=obs, p=p)
	sup  <- mll - min(mll) + qchisq(CI, 1)/2
	list(
		x   = nbcl, 
		mll = mll,
		estim.x = nbcl[which.min(mll)],
		CI      = range(nbcl[sup < 0])
	)
}

maxlik.random.manual <- function(obs, p=obs2p(obs)$p.estimates, CI=0.95, n=20) {
	nbcl <- seq(0, max.nbclones, length=n+1)[-1]
	mll  <- sapply(nbcl, mllik.random, obs=obs, p=p)
	sup  <- mll - min(mll) - qchisq(CI, 1)/2
	list(
		x   = nbcl, 
		mll = mll,
		estim.x = nbcl[which.min(mll)],
		CI      = range(nbcl[sup < 0])
	)
}


plot.distcomp <- function(frq1, frq2, ...) {
	barplot(rbind(frq1, frq2), beside=TRUE)
	
	
}

plot.distcomp.model <- function(obs, p=estim.freq(obs)$p.estimates, exact=FALSE, ...) {
	fm <- freq.model(p, max.lik.model(obs=obs, p=p, exact=exact))
	fo <- obs/sum(obs)
	
	plot.distcomp(fo, fm, ...)
}
