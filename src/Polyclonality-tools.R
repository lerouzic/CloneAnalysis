library(stats4)
library(gtools)
library(multicool)


############################## Constants ###############################
max.nbclones <- 18

col.names <-  c("G","R","Y","B")
colors    <-  c(G="springgreen2", R="firebrick2", Y="gold2", B="cyan3")
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

cache.cat <- lapply(1:max.nbclones, function(cl) { 
	# Pre-calculated and stored in a variable to save time.
		pm <- gtools::combinations(4, cl, rev(col.names), set=FALSE, repeats.allowed=TRUE)
		list(
			comb = pm, 
			choose = apply(pm, 1, multicool::multinom), 
			names = apply(pm, 1, function(x) paste0(unique(x[order(match(x, col.names))]), collapse=""))
		)
	})


predict.f.clones <- function(nbclones, p) {
	pf <- matrix(p[cache.cat[[nbclones]]$comb], ncol=nbclones)
	full.f <- apply(pf, 1, prod) * cache.cat[[nbclones]]$choose
	f <- tapply(full.f, INDEX=cache.cat[[nbclones]]$names, sum)
	
	# Numerical errors are common
	f[f<0] <- 0
	f <- f/sum(f)
	
	ans <- setNames(rep(0, length(obs.names)), nm=obs.names)
	ans[names(f)] <- f
	ans
}

predict.f <- function(lmuC, lpG, lpR, lpY) {
	p <- logitp2p(lpG, lpR, lpY)
	
	distclone <- dpois(seq_len(max.nbclones), lambda=exp(lmuC))
	distclone <- distclone / sum(distclone)
	
	colSums(do.call(rbind, lapply(seq_len(max.nbclones), predict.f.clones, p=p))*distclone)
}

maxlik.mle <- function(obs, CI=0.95) {
	mll <- function(lmuC, lpG, lpR, lpY) {
		if (exp(lmuC) > max.nbclones) return (mll(log(max.nbclones), lpG, lpR, lpY) + 1000*(lmuC-log(max.nbclones)))
		if (exp(lmuC) < 1)            return (mll(log(1), lpG, lpR, lpY) + 1/exp(lmuC)^2)
		pred.f <- predict.f(lmuC, lpG, lpR, lpY)
		-dmultinom(obs, prob=pred.f, log=TRUE)
	}
	
	fit <- try(stats4::mle(mll, start=c(lmuC=1, lpG=0, lpR=0, lpY=-2)))
	
	if (class(fit) == "try-error") return(list(coef=NA, x=NA))
	
	fit.ci <- confint(fit, level=CI)
		
	list(
		coef     = fit@coef, 
		x        = c(muC=unname(exp(fit@coef["lmuC"])), logitp2p(fit@coef["lpG"], fit@coef["lpR"], fit@coef["lpY"])),
		CI       = fit.ci,
		CI.x     = rbind(muC=unname(exp(fit.ci["lmuC",])), cbind(logitp2p(fit.ci["lpG",1], fit.ci["lpR",1], fit.ci["lpY",1]), logitp2p(fit.ci["lpG",2], fit.ci["lpR",2], fit.ci["lpY",2])))
	) 
}

plot.distcomp <- function(frq1, frq2=NULL, ...) {
	barplot(if(is.null(frq2)) frq1 else rbind(frq1, frq2), beside=TRUE, ...)
	
	
}

plot.distcomp.model <- function(obs, p=estim.freq(obs)$p.estimates, exact=FALSE, ...) {
	fm <- freq.model(p, max.lik.model(obs=obs, p=p, exact=exact))
	fo <- obs/sum(obs)
	
	plot.distcomp(fo, fm, ...)
}
