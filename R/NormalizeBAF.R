NormalizeBAF <- function(BAF)
{
	Centers <- c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1)
	names(Centers) <- c("AAAA", "AAAB", "AAB", "AB", "ABB", "ABBB", "BBBB")
	
	tmp <- sapply(Centers, function(X){ abs(BAF - X) })
	tmp2 <- apply(tmp, 1, function(X){ names(X)[min(X) == X][1] })
	tmp3 <- data.frame(Class=tmp2, BAF=BAF, Indx=1:length(tmp2), stringsAsFactors=F)

	res <- sapply(unique(tmp3$Class), function(C)
	{
		V <- tmp3$BAF[tmp3$Class %in% C]
		if(length(V) > 10) # smooth.spline needs at least 4 numbers using more than 10
		{
			V.norm <- NormalizeLRR(V, Centers[C]) # Normalize spline
		}else{ V.norm <- V }
		return(V.norm)
	})
	
	BAF_Indx <- sapply(unique(tmp3$Class), function(C){ tmp3$Indx[tmp3$Class %in% C] })
	BAF.order <- unlist(BAF_Indx)
	BAF.norm <- unlist(res)
	Res <- data.frame(BAF.order=BAF.order, BAF.norm=BAF.norm, stringsAsFactors=F)
	Res2 <- Res[order(Res$BAF.order),]
	return(Res2$BAF.norm)
}
