source("../src/Polyclonality-tools.R")

library(readxl)

use.cache <- TRUE # Do not recompule likelihoods

cc <- as.data.frame(suppressWarnings(read_excel("../data/FlybowPolyclonality.xlsx")))

cc$'Time HS' <- factor(cc$'Time HS')
cc$Days <- factor(ifelse(cc$'Days post-HS' < 21, "short", "long"))
cc$'Tumor' <- factor(cc$'Tumor')

# Filtering
cc <- cc[cc$'GFP'+cc$'mCherry'+cc$'mCitrine'+cc$'mTurquoise2'+cc$'GFP in tumors'+cc$'mCherry in tumors'+cc$'mCitrine in tumors'+cc$'mTurquoise2 in tumors' > 0,]
cc <- cc[cc$'Time HS' %in% c("30 min", "1 hr", "2 hr"),]

# Summary (obs) data
summ <- by(cc, list(cc$'Time HS', cc$'Days', cc$'Tumor'), FUN=function(ccc) {
			G <- ifelse (ccc$Tumor == "F", ccc$'GFP'       > 0, ccc$'GFP in tumors'         > 0)
			R <- ifelse (ccc$Tumor == "F", ccc$'mCherry'   > 0, ccc$'mCherry in tumors'     > 0)
			Y <- ifelse (ccc$Tumor == "F", ccc$'mCitrine'  > 0, ccc$'mCitrine in tumors'    > 0)
			B <- ifelse (ccc$Tumor == "F", ccc$'mTurquoise2'>0, ccc$'mTurquoise2 in tumors' > 0)
			obs <- factor(paste0(ifelse(G, "G", ""), ifelse(R, "R", ""), ifelse(Y, "Y", ""), ifelse(B, "B", "")), levels=obs.names)
			table(obs)
		})
names(summ) <- apply(expand.grid(dimnames(summ)), 1, paste, collapse=".")
summ <- summ[!sapply(summ, is.null)]
v.tim <- c('30 min'=0.5, '1 hr'=1, '2 hr'=2)
summ <- summ[order(v.tim[sapply(strsplit(names(summ), ".", fixed=TRUE), "[", 1)])]

if (use.cache && file.exists("./lik-cache.rds")) {
	fitsumm <- readRDS("./lik-cache.rds")
} else {
	fitsumm <- parallel::mclapply(summ, function(ss) maxlik.mle(ss), mc.cores=parallel::detectCores())
	if (use.cache) 
		saveRDS(fitsumm, file="./lik-cache.rds")
}

x.shift <- 0.05

pdf("../results/poly3.pdf", width=8, height=4)
	layout(t(1:2))
	par(mar=c(5,4,2,1))
		
	ss <- strsplit(names(summ), split=".", fixed=TRUE)
	ss.tum <- sapply(ss, "[", 3) == "T"
	ss.tim <- v.tim[sapply(ss, "[", 1)]
	ss.lon <- sapply(ss, "[", 2) == "long"
			
	# Non-tumors
	plot(NULL, xlim=c(0,2.2), ylim=c(0,1), xlab="HS duration", ylab="Estimated frequency", main= "Control clones")
	for (ccc in col.names) {
		lines(v.tim, sapply(fitsumm[!ss.tum & !ss.lon], function(x) x$x[ccc]), col=colors[ccc])
		arrows(	x0=v.tim, 
				y0=sapply(fitsumm[!ss.tum & !ss.lon], function(x) x$CI.x[ccc,1]), 
				y1=sapply(fitsumm[!ss.tum & !ss.lon], function(x) x$CI.x[ccc,2]), 
				col=colors[ccc], length=0)
	}
	for (ccc in col.names) {
		lines(v.tim+x.shift, sapply(fitsumm[!ss.tum &  ss.lon], function(x) x$x[ccc]), col=colors[ccc], lty=2)
		arrows(	x0=v.tim+x.shift, 
				y0=sapply(fitsumm[!ss.tum &  ss.lon], function(x) x$CI.x[ccc,1]), 
				y1=sapply(fitsumm[!ss.tum &  ss.lon], function(x) x$CI.x[ccc,2]), 
				col=colors[ccc], length=0, lty=2)
	}
	
	legend("topright", lty=c(1,2), legend=c("short", "long"))

	# Tumors
	plot(NULL, xlim=c(0,2.2), ylim=c(0,1), xlab="HS duration", ylab="Estimated frequency", main= "Tumors")
	for (ccc in col.names) {
		lines(v.tim, sapply(fitsumm[ss.tum & !ss.lon], function(x) x$x[ccc]), col=colors[ccc])
		arrows(	x0=v.tim, 
				y0=sapply(fitsumm[ss.tum & !ss.lon], function(x) x$CI.x[ccc,1]), 
				y1=sapply(fitsumm[ss.tum & !ss.lon], function(x) x$CI.x[ccc,2]), 
				col=colors[ccc], length=0)
	}
	for (ccc in col.names) {
		lines(v.tim+x.shift, sapply(fitsumm[ss.tum &  ss.lon], function(x) x$x[ccc]), col=colors[ccc], lty=2)
		arrows(	x0=v.tim+x.shift, 
				y0=sapply(fitsumm[ss.tum &  ss.lon], function(x) x$CI.x[ccc,1]), 
				y1=sapply(fitsumm[ss.tum &  ss.lon], function(x) x$CI.x[ccc,2]), 
				col=colors[ccc], length=0, lty=2)
	}

dev.off()


pdf("../results/poly4.pdf", width=4, height=4)

	llty4 <- c(short=1, long=2)
	
	layout(1)
	
	plot(NULL, xlim=c(0,2.1), ylim=c(0,max.nbclones), xlab="Time", ylab="Average number of clones/tumor", main= "")
	lines(v.tim, sapply(fitsumm[ ss.tum & !ss.lon], function(x) x$x["muC"]), type="b", lty=llty4["short"])
	arrows(x0=v.tim, y0=sapply(fitsumm[ ss.tum & !ss.lon], function(x) x$CI.x["muC",1]), y1= sapply(fitsumm[ ss.tum & !ss.lon], function(x) x$CI.x["muC",2]), length=0)
	lines(v.tim+x.shift, sapply(fitsumm[ ss.tum &  ss.lon], function(x) x$x["muC"]), type="b", lty=llty4["long"])
	arrows(x0=v.tim+x.shift, y0=sapply(fitsumm[ ss.tum &  ss.lon], function(x) x$CI.x["muC",1]), y1= sapply(fitsumm[ ss.tum &  ss.lon], function(x) x$CI.x["muC",2]), length=0, lty=2)
	
	legend("topright", lty=llty4, legend=names(llty4))
dev.off()


pdf("../results/poly6.pdf", width=8, height=6)

	layout(matrix(1:4, byrow=TRUE, ncol=2))
	par(mar=c(4, 4, 2, 1))
	ccol6 <- c('30 min'="gray80",'1 h'="gray50",'2 h'="gray20")
	for (tt in c("F","T")) {
		for (dr in c("short","long")) {
			x <- do.call(rbind, summ[grepl(dr, names(summ)) & grepl(tt, names(summ))])
			plot.distcomp(x, main=paste0(dr, ", ", if(tt == "F") "Control" else "Tumor"), ylim=c(0,50), las=2, col=ccol6)
			if (tt=="F" && dr == "short") legend("topleft", fill=ccol6, legend=names(ccol6))
		}
	}

dev.off()

pdf("../results/poly7.pdf", width=4, height=4)

	ccol7 <- c(F="gray", T="gray20")
	layout(1)
	par(mar=c(4, 4, 2, 1))
	plot(NULL, xlim=c(0, 2), ylim=c(0,200), xlab="HS time", ylab="Number of observations")
	for (tt in c("F","T")) {
		x <- sapply(summ[grepl(tt, names(summ))], sum)
		points(v.tim, colSums(matrix(x, nrow=2)), col=ccol7[tt], type="b")
	}
	legend("topleft", lty=1, col=rep(ccol7), legend=c("Control", "Tumor"))

dev.off()
