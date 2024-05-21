library(readxl)

lvl <- apply(expand.grid(list(c("","G"), c("","R"), c("","Y"), c("","B"))), 1, paste, collapse="") 
max.nbclones <- 7

logit    <- function(p) log(p / (1-p))
invlogit <- function(u) 1/(1+exp(-u))

logitp2freq <- function(lpG, lpR, lpY) {
	pG <- invlogit(lpG)
	pR <- (1-pG)*invlogit(lpR)
	pY <- (1-pG)*(1-pR)*invlogit(lpY)
	pB <- 1-pG - pR - pY
	c(G=unname(pG), R=unname(pR), Y=unname(pY), B=unname(pB))
}

estim.freq <- function(tab) {
	mll <- function(lpG, lpR, lpY) {
		prob <- setNames(rep(1, length(tab)), nm=names(tab))
		fp <- logitp2freq(lpG, lpR, lpY)
		prob[which(names(prob) == "")] <- 0 # Assuming observations without colors are absent
		prob <- prob * ifelse(grepl("G", names(prob)), fp["G"], 1-fp["G"])
		prob <- prob * ifelse(grepl("R", names(prob)), fp["R"], 1-fp["R"])
		prob <- prob * ifelse(grepl("Y", names(prob)), fp["Y"], 1-fp["Y"])
		prob <- prob * ifelse(grepl("B", names(prob)), fp["B"], 1-fp["B"])
		
		# Fixing numerical errors
		prob[prob < 0] <- 0
		prob <- prob/sum(prob)
		
		-dmultinom(tab, prob=prob, log=TRUE)
	}
	mm <- try(stats4::mle(mll, start=c(lpG=0, lpR=0, lpY=0)))
	if (class(mm) == "try-error") return(NA)
	list(lp.estimates = mm@coef, p.estimates=logitp2freq(mm@coef[1], mm@coef[2], mm@coef[3]))
}

cache.gridcat <- lapply(seq_len(max.nbclones), function(cl) apply(expand.grid(replicate(cl, c("G","R","Y","B"), simplify=FALSE)), 1, function(x) paste0(unique(x[order(match(x, c("G","R","Y","B")))]), collapse="")))

freq.model.exact <- function(p, nbclones) {
	ref.nn <- cache.gridcat[[nbclones]]
	prbs <- setNames(apply(expand.grid(replicate(nbclones, p, simplify=FALSE)), 1, prod), nm=ref.nn)
	sapply(lvl, function(nn) sum(prbs[ref.nn == nn]))
}

freq.model <- function(p, meanclone, exact=FALSE) {
	if (exact) {
		return(freq.model.exact(p, round(meanclone)))
	}
	distclone <- dpois(seq_len(max.nbclones), lambda=meanclone)
	distclone <- distclone / sum(distclone)
	
	colSums(do.call(rbind, lapply(seq_len(max.nbclones), freq.model.exact, p=p))*distclone)
}

explore.lik.model <- function(obs, p=estim.freq(obs)$p.estimates, exact=FALSE, n=20) {
	nbcl <- if(exact) seq_len(max.nbclones) else seq(0, max.nbclones, length=n+1)[-1]
	llik <- sapply(nbcl, function(cl) dmultinom(obs, prob=freq.model(p, meanclone=cl, exact=exact), log=TRUE))
	list(nbcl=nbcl, loglik=llik)
}

max.lik.model <- function(obs, p=estim.freq(obs)$p.estimates, exact=FALSE, n=20) {
	ll <- explore.lik.model(obs=obs, p=p, exact=exact, n=n)
	return(ll$nbcl[which.max(ll$loglik)])
}

plot.distcomp <- function(frq1, frq2, ...) {
	barplot(rbind(frq1, frq2), beside=TRUE)
	
	
}

plot.distcomp.model <- function(obs, p=estim.freq(obs)$p.estimates, exact=FALSE, ...) {
	fm <- freq.model(p, max.lik.model(obs=obs, p=p, exact=exact))
	fo <- obs/sum(obs)
	
	plot.distcomp(fo, fm, ...)
}

cc <- as.data.frame(suppressWarnings(read_excel("../data/FlybowPolyclonality.xlsx")))

cc$'Time HS' <- factor(cc$'Time HS')
cc$Days <- factor(ifelse(cc$'Days post-HS' < 21, "short", "long"))
cc$'Tumor' <- factor(cc$'Tumor')

# Filtering
cc <- cc[cc$'GFP'+cc$'mCherry'+cc$'mCitrine'+cc$'mTurquoise2'+cc$'GFP in tumors'+cc$'mCherry in tumors'+cc$'mCitrine in tumors'+cc$'mTurquoise2 in tumors' > 0,]
cc <- cc[cc$'Time HS' %in% c("30 min", "1 hr", "2 hr"),]


summ <- by(cc, list(cc$'Time HS', cc$'Days', cc$'Tumor'), FUN=function(ccc) {
			G <- ifelse (ccc$Tumor == "F", ccc$'GFP'       > 0, ccc$'GFP in tumors'         > 0)
			R <- ifelse (ccc$Tumor == "F", ccc$'mCherry'   > 0, ccc$'mCherry in tumors'     > 0)
			Y <- ifelse (ccc$Tumor == "F", ccc$'mCitrine'  > 0, ccc$'mCitrine in tumors'    > 0)
			B <- ifelse (ccc$Tumor == "F", ccc$'mTurquoise2'>0, ccc$'mTurquoise2 in tumors' > 0)
			pap <- factor(paste0(ifelse(G, "G", ""), ifelse(R, "R", ""), ifelse(Y, "Y", ""), ifelse(B, "B", "")), levels=lvl)
			table(pap)
		})
names(summ) <- apply(expand.grid(dimnames(summ)), 1, paste, collapse=".")
summ <- summ[!sapply(summ, is.null)]

pdf("../results/poly1.pdf", width=14, height=6)
	par(mar=c(7,4,1,1))
	barplot(do.call(cbind, lapply(summ, function(x) sapply(setNames(nm=c("G","R","Y","B")), function(c) sum(x[grepl(c,names(x))])))), beside=TRUE, las=2, col=c("green","red","yellow","blue"))
dev.off()

pdf("../results/poly2.pdf", width=14, height=6)
	par(mar=c(7,4,1,1))
	bb <- barplot(do.call(cbind, lapply(summ, function(x) sapply(setNames(nm=1:4), function(c) sum(x[nchar(names(x))==c])/sum(x)))), beside=FALSE, las=2)
	text(x=bb[10], y=c(0.2, 0.5, 0.8, 0.9), as.character(1:4), col="white")
dev.off()

pdf("../results/poly3.pdf", width=8, height=8)
	layout(rbind(c(1,1),c(2,3)))
	par(mar=c(7,4,1,1))
	
	v.tim <- c('30 min'=0.5, '1 hr'=1, '2 hr'=2)
	summ <- summ[order(v.tim[sapply(strsplit(names(summ), ".", fixed=TRUE), "[", 1)])]
	
	ff <- lapply(summ, estim.freq)
	
	ss <- strsplit(names(summ), split=".", fixed=TRUE)
	ss.tum <- sapply(ss, "[", 3) == "T"
	ss.tim <- v.tim[sapply(ss, "[", 1)]
	ss.lon <- sapply(ss, "[", 2) == "long"
		
	bb <- barplot(sapply(ff, function(x) if (is.na(x[1])) rep(NA, 4) else x$p.estimates), col=c("green","red","yellow","blue"), las=2)
	
	
	plot(NULL, xlim=c(0,2), ylim=c(0,1), xlab="Time", ylab="Frequency", main= "Non-tumorous clones")
	lines(v.tim, sapply(ff[!ss.tum & !ss.lon], function(x) x$p.estimates["G"]), col="green")
	lines(v.tim, sapply(ff[!ss.tum & !ss.lon], function(x) x$p.estimates["R"]), col="red")
	lines(v.tim, sapply(ff[!ss.tum & !ss.lon], function(x) x$p.estimates["Y"]), col="yellow")
	lines(v.tim, sapply(ff[!ss.tum & !ss.lon], function(x) x$p.estimates["B"]), col="blue")
	lines(v.tim, sapply(ff[!ss.tum & ss.lon],  function(x) x$p.estimates["G"]), col="green",  lty=2)
	lines(v.tim, sapply(ff[!ss.tum & ss.lon],  function(x) x$p.estimates["R"]), col="red",    lty=2)
	lines(v.tim, sapply(ff[!ss.tum & ss.lon],  function(x) x$p.estimates["Y"]), col="yellow", lty=2)
	lines(v.tim, sapply(ff[!ss.tum & ss.lon],  function(x) x$p.estimates["B"]), col="blue",   lty=2)
	
	legend("topleft", lty=c(1,2), legend=c("short", "long"))

	plot(NULL, xlim=c(0,2), ylim=c(0,1), xlab="Time", ylab="Frequency", main= "Tumors")
	lines(v.tim, sapply(ff[ss.tum & !ss.lon], function(x) x$p.estimates["G"]), col="green")
	lines(v.tim, sapply(ff[ss.tum & !ss.lon], function(x) x$p.estimates["R"]), col="red")
	lines(v.tim, sapply(ff[ss.tum & !ss.lon], function(x) x$p.estimates["Y"]), col="yellow")
	lines(v.tim, sapply(ff[ss.tum & !ss.lon], function(x) x$p.estimates["B"]), col="blue")
	lines(v.tim, sapply(ff[ss.tum & ss.lon],  function(x) x$p.estimates["G"]), col="green",  lty=2)
	lines(v.tim, sapply(ff[ss.tum & ss.lon],  function(x) x$p.estimates["R"]), col="red",    lty=2)
	lines(v.tim, sapply(ff[ss.tum & ss.lon],  function(x) x$p.estimates["Y"]), col="yellow", lty=2)
	lines(v.tim, sapply(ff[ss.tum & ss.lon],  function(x) x$p.estimates["B"]), col="blue",   lty=2)

dev.off()


pdf("../results/poly4.pdf", width=8, height=4)

	layout(t(1:2))

	plot(NULL, xlim=c(0,2), ylim=c(0,max.nbclones), xlab="Time", ylab="Average number of clones/image", main= "Non-tumorous clones")
	lines(v.tim, sapply(summ[!ss.tum & !ss.lon], function(ss) max.lik.model(ss)))
	lines(v.tim, sapply(summ[!ss.tum &  ss.lon], function(ss) max.lik.model(ss)), lty=2)
	
	legend("topleft", lty=c(1,2), legend=c("short", "long"))

	plot(NULL, xlim=c(0,2), ylim=c(0,max.nbclones), xlab="Time", ylab="Average number of clones/image", main= "Tumors")
	lines(v.tim, sapply(summ[ ss.tum & !ss.lon], function(ss) max.lik.model(ss)))
	lines(v.tim, sapply(summ[ ss.tum &  ss.lon], function(ss) max.lik.model(ss)), lty=2)

dev.off()
