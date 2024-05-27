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

pdf("../results/poly3.pdf", width=8, height=8)
	layout(rbind(c(1,1),c(2,3)))
	par(mar=c(7,4,1,1))
		
	ss <- strsplit(names(summ), split=".", fixed=TRUE)
	ss.tum <- sapply(ss, "[", 3) == "T"
	ss.tim <- v.tim[sapply(ss, "[", 1)]
	ss.lon <- sapply(ss, "[", 2) == "long"
		
	bb <- barplot(sapply(fitsumm, function(x) if (is.na(x[1])) rep(NA, 4) else x$x[names(colors)]), col=colors, las=2)
	
	# Non-tumors
	plot(NULL, xlim=c(0,2.2), ylim=c(0,1), xlab="Time", ylab="Frequency", main= "Non-tumorous clones")
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
	plot(NULL, xlim=c(0,2.2), ylim=c(0,1), xlab="Time", ylab="Frequency", main= "Tumors")
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


pdf("../results/poly4.pdf", width=8, height=4)

	layout(t(1:2))

	plot(NULL, xlim=c(0,2), ylim=c(0,max.nbclones), xlab="Time", ylab="Average number of clones/image", main= "Non-tumorous clones")
	lines(v.tim, sapply(fitsumm[!ss.tum & !ss.lon], function(x) x$x["muC"]), type="b")
	arrows(x0=v.tim, y0=sapply(fitsumm[!ss.tum & !ss.lon], function(x) x$CI.x["muC",1]), y1= sapply(fitsumm[!ss.tum & !ss.lon], function(x) x$CI.x["muC",2]), length=0)
	lines(v.tim+x.shift, sapply(fitsumm[!ss.tum &  ss.lon], function(x) x$x["muC"]), type="b", lty=2)
	arrows(x0=v.tim+x.shift, y0=sapply(fitsumm[!ss.tum &  ss.lon], function(x) x$CI.x["muC",1]), y1= sapply(fitsumm[!ss.tum &  ss.lon], function(x) x$CI.x["muC",2]), length=0, lty=2)
		
	legend("bottomright", lty=c(1,2), legend=c("short", "long"))

	plot(NULL, xlim=c(0,2), ylim=c(0,max.nbclones), xlab="Time", ylab="Average number of clones/image", main= "Tumors")
	lines(v.tim, sapply(fitsumm[ ss.tum & !ss.lon], function(x) x$x["muC"]), type="b")
	arrows(x0=v.tim, y0=sapply(fitsumm[ ss.tum & !ss.lon], function(x) x$CI.x["muC",1]), y1= sapply(fitsumm[ ss.tum & !ss.lon], function(x) x$CI.x["muC",2]), length=0)
	lines(v.tim+x.shift, sapply(fitsumm[ ss.tum &  ss.lon], function(x) x$x["muC"]), type="b", lty=2)
	arrows(x0=v.tim+x.shift, y0=sapply(fitsumm[ ss.tum &  ss.lon], function(x) x$CI.x["muC",1]), y1= sapply(fitsumm[ ss.tum &  ss.lon], function(x) x$CI.x["muC",2]), length=0, lty=2)

dev.off()


pdf("../results/poly5.pdf", width=8, height=10)

	layout(matrix(1:length(summ), byrow=TRUE, ncol=2))
	par(mar=c(4, 4, 2, 1))
	for (nn in names(summ)) {
		ccol <- c(obs="darkgray", pred="lightgray")
		obs  <- summ[[nn]]
		pred <- do.call(predict.f, as.list(fitsumm[[nn]]$coef))
		pval <- chisq.test(obs[-1], p=pred[-1], simulate.p.value=TRUE, B=1e4)$p.value
		stars <- if(pval > 0.05) "" else if (pval > 0.01) "*" else if (pval > 0.001) "**" else "***"
		plot.distcomp(obs, sum(obs)*pred, main=paste0(nn, if (stars != "") paste0(" (", stars, ")")), las=2, col=ccol)
		if (nn == names(summ)[1])
			legend("topright", fill=ccol, legend=c("Observed", "Predicted"))
	}

dev.off()
