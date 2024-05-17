library(readxl)

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

