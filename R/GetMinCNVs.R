GetMinCNVs <- function(NCases=10000, NControl=25000, NTests=800, nCNVs = c(0,1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30))
{
	Total <- sum(NCases, NControl)

	MultTest <- 0.05/NTests

	df <- data.frame(State=rep(0, Total), CN=rep(2, Total), stringsAsFactors=F)
	df$State[1:NCases] <- 1 # Cases

	Res <- sapply(1:length(nCNVs), function(i) # Adding CNVs at cases
	{
		df$CN[1:nCNVs[i]] <- 1
		res2 <- sapply(1:length(nCNVs), function(i2) # Adding CNVs at control
		{		
			df$CN[(NCases+1):(NCases+nCNVs[i2])] <- 1
			df$CN <- factor(df$CN, levels=c(2,1))
			mylogit <- glm(State ~ CN, data = df, family = "binomial")
			OR <- exp(coef(mylogit))[2]
			pvalue <- summary(mylogit)$coefficients[2,4]
		
			if(pvalue < MultTest){ Sig <- TRUE }else{ Sig <- FALSE}
	
			res <- data.frame(Cases=nCNVs[i], Control=nCNVs[i2], pvalue=pvalue, OR=OR, Sig=Sig)
			cat(xtabs(~State+CN, df), "\n")
			
			cat(nCNVs[i], nCNVs[i2], pvalue, OR, "\n")
			return(res)
		})
		df$CN =rep(2, Total)

		res3 <- MatrixOrList2df(res2)
		return(res3)
	})
	Res2 <- MatrixOrList2df(Res)	
	Res3 <- Res2[order(Res2$Control, Res2$Cases),]

	tmp <- subset(Res3, Sig)
	CNVs <- tmp$Cases + tmp$Control
	Min <- min(CNVs)
	cat("Minimum number of CNVs is:", Min,"\n")
	cat("NCases:", NCases, ", NControl:",NControl,"NTests:",NTests, "\n")

	return(Res3)
}
	
