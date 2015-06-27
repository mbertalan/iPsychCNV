MergeCNVs <- function(df, MaxNumSNPs=50)
{
	tmp2 <- subset(df, CN != 2)
	Test2 <- sapply(unique(tmp2$ID), function(IDs) # Sample Loop	
	{
		#cat("Sample: ", IDs, "\n")
		tmp3 <- subset(tmp2, ID %in% IDs) 	 
		Test <- sapply(unique(tmp3$Chr), function(CHRs) # Chr loop
		{
			#cat("Chr: ", CHRs, "\n")
			tmp4 <- subset(tmp3, Chr %in% CHRs) 
			if(nrow(tmp4) > 1)
			{		
				# If the next CNV has distance shorter than MaxNumSNPs.
				M2 <- sapply(1:(nrow(tmp4[,c("StartIndx","StopIndx")])-1), function(X){ if((as.numeric(tmp4[(X+1),1]) - as.numeric(tmp4[(X),2])) < MaxNumSNPs){ c(X, (X+1)) } })
				library(igraph)
				g2 <- graph(unlist(M2), directed=FALSE)
				k <- clusters(g2)
				Groups <- k$membership
				names(Groups) <- 1:length(Groups)

				Size <- k$csize
				names(Size) <- 1:length(unique(k$membership))

				Groups <- Groups[Groups %in% names(Size[Size > 1])]
				Merge <- sapply(unique(Groups), function(X){ paste(names(Groups[Groups == X]), sep="", collapse="_") })
				names(Merge) <- Merge
			

				CheckingMerge <- sapply(Merge, function(X)
				{
					#cat("2-) Merging: ", X, "\n")
					# Getting the rows that should be merged.
					Indx <- as.numeric(unlist(strsplit(X, "_")))
					Indx <- unique(as.numeric(Indx))
					#cat("3-) Indx: ", Indx, "\n")
					
					if(length(Indx) > 1)
					{
						# CN 1 and 3 should not be merged. If more than 1 CN, select the most common.
						tmp <- subset(tmp4[Indx,], CN == as.numeric(names(sort(table(tmp4[Indx,"CN"]), decreasing=T)[1])))
						Indx <- as.numeric(row.names(tmp))
						
						#cat("4-) Great ! Join ", Indx, "\n")
						NewStart <- sort(tmp4[Indx,"Start"])[1]
						NewStop <- sort(tmp4[Indx,"Stop"], decreasing=TRUE)[1]
						NewStartIndx <- sort(tmp4[Indx,"StartIndx"])[1]
						NewStopIndx <- sort(tmp4[Indx,"StopIndx"], decreasing=TRUE)[1]
						NewLength <-  NewStop - NewStart
						NewNumSNPs <- NewStopIndx - NewStartIndx
						df <- data.frame(Start=NewStart, Stop=NewStop, StartIndx=NewStartIndx, StopIndx=NewStopIndx, Length=NewLength, NumSNPs=NewNumSNPs)
						# Original DF
						df2 <- tmp4[Indx[1], !colnames(tmp4) %in% colnames(df)] # Original data without the modified merged info.
						df3 <- cbind(df, df2)
						return(df3)
					}
					#cat("5-) Done ! ", Indx, "\n")
				})
				Merged <- MatrixOrList2df(CheckingMerge)
				# Collecting the CNVs that are merged
				Indx <- sapply(rownames(Merged), function(X){ unlist(strsplit(X, "_")) }) 
				Indx <- as.numeric(unlist(Indx))
				df <- tmp4[Indx*-1,]
				Merged <- rbind(Merged, df)
				# If there was nothing to merge (all CNVs are too far way from each other.
				if(nrow(Merged) == 0)
				{
					Merged <- tmp4
				}
			}
			else
			{ 
				#cat("Northing to merge ! \n")
				Merged <- tmp4
			}							
			return(Merged)	
		})
		Test2 <- MatrixOrList2df(Test)
		return(Test2)
	})
	Test2 <- MatrixOrList2df(Test2)
	Test2 <- Test2[order(Test2$StartIndx),]
	return(Test2)
}
