source("../src/Polyclonality-tools.R")

library(readxl)

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
	
pdf("../results/poly3.pdf", width=8, height=8)
	layout(rbind(c(1,1),c(2,3)))
	par(mar=c(7,4,1,1))
	
	ff <- lapply(summ, obs2p)
	
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

	fitsumm <- lapply(summ, function(ss) maxlik.random.manual(ss))

	plot(NULL, xlim=c(0,2), ylim=c(0,max.nbclones), xlab="Time", ylab="Average number of clones/image", main= "Non-tumorous clones")
	lines(v.tim, sapply(fitsumm[!ss.tum & !ss.lon], function(x) x$estim.x), type="b")
	arrows(x0=v.tim, y0=sapply(fitsumm[!ss.tum & !ss.lon], function(x) x$CI[1]), y1= sapply(fitsumm[!ss.tum & !ss.lon], function(x) x$CI[2]), length=0)
	lines(v.tim+0.05, sapply(fitsumm[!ss.tum &  ss.lon], function(x) x$estim.x), type="b", lty=2)
	arrows(x0=v.tim+0.05, y0=sapply(fitsumm[!ss.tum &  ss.lon], function(x) x$CI[1]), y1= sapply(fitsumm[!ss.tum &  ss.lon], function(x) x$CI[2]), length=0, lty=2)
		
	legend("topleft", lty=c(1,2), legend=c("short", "long"))

	plot(NULL, xlim=c(0,2), ylim=c(0,max.nbclones), xlab="Time", ylab="Average number of clones/image", main= "Tumors")
	lines(v.tim, sapply(fitsumm[ ss.tum & !ss.lon], function(x) x$estim.x), type="b")
	arrows(x0=v.tim, y0=sapply(fitsumm[ ss.tum & !ss.lon], function(x) x$CI[1]), y1= sapply(fitsumm[ ss.tum & !ss.lon], function(x) x$CI[2]), length=0)
	lines(v.tim+0.05, sapply(fitsumm[ ss.tum &  ss.lon], function(x) x$estim.x), type="b", lty=2)
	arrows(x0=v.tim+0.05, y0=sapply(fitsumm[ ss.tum &  ss.lon], function(x) x$CI[1]), y1= sapply(fitsumm[ ss.tum &  ss.lon], function(x) x$CI[2]), length=0, lty=2)

dev.off()
