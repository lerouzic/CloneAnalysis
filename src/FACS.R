library(readxl)

p.stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))


############### FACS data #######################

dd <- as.data.frame(suppressWarnings(read_excel("../data/FACSdata.xlsx", col_types=c("text", "numeric", "numeric", "numeric"))))

dd$Genotype[dd$Genotype %in% c("Control 1", "Control 2")] <- "Control 1&2" 
dd$Genotype[dd$Genotype %in% c("control 1 bis", "Control 3", "Control 4")] <- "Control Other"
dd <- dd[dd$Genotype !=  "Control Other",]

colnames(dd)[2] <- "All" # Easier to handle
dd <- dd[-219,] # Line 219 is nonsensical

# Due to some preliminary treatment (?), counts are not all integers
dd$All <- round(dd$All)
dd$GFP <- round(dd$GFP)

sn <- strsplit(dd$Genotype, split=" ")

dd$Genotype <- factor(
	paste0(sapply(sn, "[", 1), ".", sapply(sn, "[", 2), sapply(sn, function(ss) if (length(ss) > 3) paste0(".", ss[4]) else ""))
	)
dd$Genotype <- relevel(dd$Genotype, "Control.1&2")

dd.lm1  <- lm(I(dd$GFP/dd$All) ~ dd$Genotype)
# glm1 <- glm(cbind(dd$GFP, dd$All-dd$GFP) ~ dd$Genotype, family=binomial()) # Not interesting
dd.glm2 <- glm(cbind(dd$GFP, dd$All-dd$GFP) ~ dd$Genotype, family=quasibinomial())

dd.s1 <- summary(dd.lm1)
dd.s3 <- summary(dd.glm2)

dd.ans <- data.frame(
	Estimate.1 = dd.s1$coefficients[,1], 
	Estimate.3 = dd.s3$coefficients[,1],
	Std.err.1  = dd.s1$coefficients[,2], 
	Std.err.3  = dd.s3$coefficients[,2], 
	p.1        = dd.s1$coefficients[,4], 
	p.3        = dd.s3$coefficients[,4], 
	pa.1       = p.adjust(dd.s1$coefficients[,4], "holm"), 
	pa.3       = p.adjust(dd.s3$coefficients[,4], "holm")
	)
	
rownames(dd.ans) <- gsub("dd\\$Genotype", "", rownames(dd.s1$coefficients))

dd.ans$st.1 <- p.stars(dd.ans$pa.1)
dd.ans$st.3 <- p.stars(dd.ans$pa.3)
dd.ans$st.old <- sapply(rownames(dd.ans), function(rn) {w <- which(dd$Genotype == rn); if (length(w) > 0) as.character(sn[[w[1]]][3]) else ""})

write.table(format(dd.ans, digits=4), "../results/FACS.txt", quote=FALSE, sep="\t")

##################### Motility data #######################
	
dd2 <- as.data.frame(suppressWarnings(read_excel("../data/DataMotility.xlsx", col_types=c("text", "numeric", "numeric", "numeric"))))

dd2$Genotype <- gsub(" ", ".", dd2$Genotype)
colnames(dd2)[2] <- "All"

dd2$Genotype[dd2$Genotype == "Control.1"] <- "Control"
dd2$Genotype <- factor(dd2$Genotype)
dd2$Genotype <- relevel(dd2$Genotype, "Control")

dd2.lm1  <- lm(I(dd2$GFP/dd2$All) ~ dd2$Genotype)
dd2.glm2 <- glm(cbind(dd2$GFP, dd2$All-dd2$GFP) ~ dd2$Genotype, family=quasibinomial())

dd2.s1 <- summary(dd2.lm1)
dd2.s3 <- summary(dd2.glm2)

dd2.ans <- data.frame(
	Estimate.1 = dd2.s1$coefficients[,1], 
	Estimate.3 = dd2.s3$coefficients[,1],
	Std.err.1  = dd2.s1$coefficients[,2], 
	Std.err.3  = dd2.s3$coefficients[,2], 
	p.1        = dd2.s1$coefficients[,4], 
	p.3        = dd2.s3$coefficients[,4], 
	pa.1       = p.adjust(dd2.s1$coefficients[,4], "holm"), 
	pa.3       = p.adjust(dd2.s3$coefficients[,4], "holm")
	)

rownames(dd2.ans) <- gsub("dd2\\$Genotype", "", rownames(dd2.s1$coefficients))

dd2.ans$st.1 <- p.stars(dd2.ans$pa.1)
dd2.ans$st.3 <- p.stars(dd2.ans$pa.3)

write.table(format(dd2.ans, digits=4), "../results/Motility.txt", quote=FALSE, sep="\t")
