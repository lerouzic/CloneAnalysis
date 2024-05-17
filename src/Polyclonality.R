library(readxl)

logit    <- function(p) log(p / (1-p))
invlogit <- function(u) 1/(1+exp(-u))

logitp2freq <- function(lpG, lpR, lpY) {
	pG <- invlogit(lpG)
	pR <- (1-pG)*invlogit(lpR)
	pY <- (1-pG)*(1-pR)*invlogit(lpY)
	pB <- 1-pG - pR - pY
	c(G=pG, R=pR, Y=pY, B=pB)
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
		-dmultinom(tab, prob=prob, log=TRUE)
	}
	mm <- stats4::mle(mll, start=c(lpG=0, lpR=0, lpY=0))
	list(lp.estimates = mm@coef, p.estimates=logitp2freq(mm@coef[1], mm@coef[2], mm@coef[3]))
}

cc <- as.data.frame(suppressWarnings(read_excel("../data/FlybowPolyclonality.xlsx")))

cc$'Time HS' <- factor(cc$'Time HS')
cc$Days <- factor(ifelse(cc$'Days post-HS' <= 7, "short", "long"))
cc$'Tumor' <- factor(cc$'Tumor')

lvl <- apply(expand.grid(list(c("","G"), c("","R"), c("","Y"), c("","B"))), 1, paste, collapse="") 

summ <- by(cc, list(cc$'Time HS', cc$'Days', cc$'Tumor'), FUN=function(ccc) {
			G <- ccc$'GFP'      > 0 | ccc$'GFP in tumors'      > 0
			R <- ccc$'mCherry'  > 0 | ccc$'mCherry in tumors'  > 0
			Y <- ccc$'mCitrine' > 0 | ccc$'mCitrine in tumors' > 0
			B <- ccc$'mTurquoise2'>0| ccc$'mTurquoise2 in tumors' > 0
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

