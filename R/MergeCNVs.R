MergeCNVs <- function(df, MaxNumSNPs=4)
{
	MaxNumSNPs <- 2 # Have to fix this
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
				M <- sapply(tmp4$StopIndx, function(Stop){ sapply(tmp4$StartIndx, function(Start) { abs(Start - Stop ) < MaxNumSNPs || abs(Stop - Start ) < MaxNumSNPs }) })
				diag(M) <- FALSE
				rownames(M) <- 1:nrow(M)
				colnames(M) <- 1:nrow(M)

				tmp5 <- apply(M, 1, function(Z) { names(Z[Z]) })
				if(length(tmp5) != 0) 
				{
					Merge <- sapply(1:length(tmp5), function(X){ paste( sort(as.numeric(c(names(tmp5[X]), tmp5[[X]]))), sep="", collapse="_") })
					MergeSize <- sapply(unique(Merge), function(X){ length(unlist(strsplit(X, "_"))) })
					MergeSize <- MergeSize[MergeSize > 1] # remove CNVs that will not be merged
				
					Mergetmp <- sapply(names(MergeSize), function(X){ unlist(strsplit(X, "_")) })

					# Checking for multiple merges				
					Rep <- table(unlist(strsplit(names(MergeSize), "_")))
					if(sum(Rep > 1) > 0)
					{
						library(igraph)
						g2 <- graph(c(sort(as.numeric(as.vector(Mergetmp)))), directed=FALSE)
						k <- clusters(g2)
						Groups <- k$membership
						names(Groups) <- 1:length(Groups)

						Size <- k$csize
						names(Size) <- 1:length(unique(k$membership))

						Groups <- Groups[Groups %in% names(Size[Size > 1])]
						Merge <- sapply(unique(Groups), function(X){ paste(names(Groups[Groups == X]), sep="", collapse="_") })
						names(Merge) <- Merge
					}
					else
					{
						Merge <- names(MergeSize)
					}

					names(Merge) <- Merge
					#cat("1-) Merge1: ", Merge, " Names Merge1:", names(Merge), "\n")

					CheckingMerge <- sapply(Merge, function(X)
					{
						#cat("2-) Merging: ", X, "\n")
						Indx <- as.numeric(unlist(strsplit(X, "_")))
						Indx <- unique(as.numeric(Indx))
						#cat("3-) Indx: ", Indx, "\n")
						if(length(unique(tmp4[Indx,]$CN)) == 1 & length(Indx) > 1)
						{
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
				}
				else
				{ 
					#cat("Northing to merge ! \n")
					Merged <- tmp4
				}							
				return(Merged)
			}
			else
			{
				return(tmp4) # If there is only one CNV
			}
		})
		Test2 <- MatrixOrList2df(Test)
		return(Test2)
	})
	Test2 <- MatrixOrList2df(Test2)
	return(Test2)
}
