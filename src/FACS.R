library(readxl)

dd <- as.data.frame(read_excel("../data/FACSdata.xlsx", col_types=c("text", "numeric", "numeric", "numeric")))

dd$Genotype[dd$Genotype %in% c("Control 1", "Control 2")] <- "Control 1&2" 
dd$Genotype[dd$Genotype %in% c("control 1 bis", "Control 3", "Control 4")] <- "Control Other"

dd <- dd[-219,] # Line 219 is nonsensical

colnames(dd)[2] <- "All"
sn <- strsplit(dd$Genotype, split=" ")

dd$Genotype <- factor(paste(sapply(sn, "[", 1), sapply(sn, "[", 2), sep="."))
dd$Genotype <- relevel(dd$Genotype, "Control.1&2")

lm1  <- lm(I(dd$GFP/dd$All) ~ dd$Genotype)
glm1 <- glm(cbind(dd$GFP, dd$All-dd$GFP) ~ dd$Genotype, family=binomial())
glm2 <- glm(cbind(dd$GFP, dd$All-dd$GFP) ~ dd$Genotype, family=quasibinomial())

s1 <- summary(lm1)
s2 <- summary(glm1)
s3 <- summary(glm2)

ans <- data.frame(
	Estimate.1 = s1$coefficients[,1], 
	Estimate.2 = s2$coefficients[,1],
	Estimate.3 = s3$coefficients[,1],
	Std.err.1  = s1$coefficients[,2], 
	Std.err.2  = s2$coefficients[,2], 
	Std.err.3  = s3$coefficients[,2], 
	p.1        = s1$coefficients[,4], 
	p.2        = s2$coefficients[,4], 
	p.3        = s3$coefficients[,4], 
	pa.1       = p.adjust(s1$coefficients[,4], "holm"), 
	pa.2       = p.adjust(s2$coefficients[,4], "holm"), 
	pa.3       = p.adjust(s3$coefficients[,4], "holm")
	)
rownames(ans) <- gsub("dd\\$Genotype", "", rownames(s2$coefficients))

p.stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))

ans$st.1 <- p.stars(ans$pa.1)
ans$st.2 <- p.stars(ans$pa.2)
ans$st.3 <- p.stars(ans$pa.3)

write.table(ans, "../results/FACS.txt", quote=FALSE, sep="\t")

	
