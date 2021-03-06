##' CompressCNVs: unknown
##'
##' Unknown - description?
##' @title CompressCNVs
##' @param tmp4: Unknown.
##' @param OverlapCutoff: Unknown.
##' @param Cores: Number of cores used, default = 1.
##' @return Hotspots.
##' @author Marcelo Bertalan, Louise K. Hoeffding.
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples Unknown

CompressCNVs <- function(tmp4, OverlapCutoff, Cores)
{
	library(plyr)
	library(parallel)

	CHRs <- unique(tmp4$Chr)

	tmp5 <- mclapply(CHRs, mc.cores=Cores, mc.preschedule = FALSE, function(X) #tmp5 <- sapply(CHRs, function(X)
	{

		tmp <- subset(tmp4, Chr %in%  X)
		rownames(tmp) <- tmp$NewName
		cat("Chr\t",X, "\t", nrow(tmp), "\r")
		if(nrow(tmp) > 1)
		{
			tmp2 <- FindHighFreqCNVs(tmp,OverlapCutoff)
			colnames(tmp2) <- tmp$NewName
			rownames(tmp2) <- tmp$NewName
			
			tmp3 <- apply(tmp2, 1, function(Z) { names(Z[Z]) })
			
			# Check if all regions match with itself.
			if(sum(unlist(lapply(tmp3, length)) == 0) > 0){ stop("The region does not match with itself. Please check your start, stop and length of your region\n") }
				
			tmp5 <- lapply(tmp3, function(L)
			{
				NewStart <- sort(as.numeric(tmp[L,]$Start))[1]
				NewStop <- sort(as.numeric(tmp[L,]$Stop), decreasing=T)[1]
				#NewStart <- as.integer(quantile(tmp[L,]$Start, probs=1:100/100)[20]) # 20%
				#NewStop <- as.integer(quantile(tmp[L,]$Stop, , probs=1:100/100)[80]) # 85%
				NewChr <- X
				NewLength <- NewStop - NewStart
				NewCount <- sum(as.numeric(tmp[L,]$Count))
				NewSource <- unique(tmp$Source)
				NewName <- paste(NewChr, NewStart, NewStop, NewLength, sep="_", collapse="")
				Newdf <- data.frame("Chr"=NewChr,"Start"=NewStart, "Stop"=NewStop, "Length"=NewLength, "Count"=NewCount, Source=NewSource, "NewName"=NewName, stringsAsFactors=F)
			})
			df3 <- ldply(tmp5, data.frame)
			df3 <- df3[, !names(df3) %in% ".id"]
			df3 <- df3[!duplicated(df3[,c("Start","Stop","Chr")]),]
			return(df3)
		}
		else
		{
			tmp$NewName <- rownames(tmp) #paste(rownames(tmp), tmp$Source, sep="_", collapse="")
			return(tmp)
		}
		cat("Chr\t", X, "\tdone.\r")
	})
	#save(tmp5, file="tmp5.RData")
	tmp <- MatrixOrList2df(tmp5)
	tmp <- tmp[, !names(tmp) %in% ".id"]
	return(tmp)
}
