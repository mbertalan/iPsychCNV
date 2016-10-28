##' MockData: Unknown. 
##'
##' Specifically designed to handle noisy data from amplified DNA on phenylketonuria (PKU) cards. The function is a pipeline using many subfunctions.
##' @title MockData
##' @param Heterozygosity, default = 10.
##' @return List with BAF frequency for each copy-number state.
##' @author Marcelo Bertalan. 
##' @source \url{http://biopsych.dk/iPsychCNV}
##' @export
##' @examples 
##' BAF_list <- MakeBAF(Heterozygosity = 10)

MakeBAF <- function(Heterozygosity = 10)
{
	BAFs <- seq(from=0, to=1, by=0.01) # 101
	BAF_Basic <- rep(0.00001, 101)
	names(BAF_Basic) <- BAFs

	##
	BAFs <- seq(from=0, to=1, by=0.01)
	# 0 or 1
	Y = (100 - Heterozygosity)/36

	BAF_Normal <- rep(0, 101)
	BAF_Normal[c(1:2)] <- BAF_Normal[c(1:2)] + (Y *8)  # 40%
	BAF_Normal[c(3:4)] <- BAF_Normal[c(3:4)] + Y # 5%
	
	BAF_Normal[c(98:99)] <- BAF_Normal[c(98:99)] + Y # 5%
	BAF_Normal[c(100:101)] <- BAF_Normal[c(100:101)] + (Y*8) # 40%
		
	BAF_Normal[c(49:52)] <- BAF_Normal[c(49:52)] + (Heterozygosity/5)
	BAF_Normal[c(50:51)] <- BAF_Normal[c(50:51)] + (Heterozygosity/10)
	

	# BAF Del prob
	BAF_Del <- BAF_Basic
	BAF_Del[1:2] <- BAF_Del[1:2] + 1
	BAF_Del[100:101] <- BAF_Del[100:101] + 1

	# BAF CN=0
	BAF_CN0 <- BAF_Basic
	
	# BAF Dup prob
	BAF_Dup <-  BAF_Basic
	BAF_Dup[1:2] <- BAF_Dup[1:2] + (Y*9)
	BAF_Dup[100:101] <- BAF_Dup[100:101] + (Y*9)
	#BAF_Dup[30:35] <- BAF_Dup[30:35] + 0.05 
	BAF_Dup[32:33] <- BAF_Dup[32:33] + Heterozygosity/4
	#BAF_Dup[65:70] <- BAF_Dup[65:70] + 0.05
	BAF_Dup[67:68] <- BAF_Dup[67:68] + Heterozygosity/4

	# BAF CN=4
	BAF_CN4 <- BAF_Basic
	BAF_CN4[1:2] <- BAF_CN4[1:2] + (Y*9)
	BAF_CN4[100:101] <- BAF_CN4[100:101] + (Y*9)
	BAF_CN4[24:26] <- BAF_CN4[25:26] + Heterozygosity/6
	BAF_CN4[50:51] <- BAF_CN4[50:51] + Heterozygosity/6
	BAF_CN4[74:76] <- BAF_CN4[75:76] + Heterozygosity/6
	
	tmp <- list(BAF_CN0=BAF_CN0, BAF_CN4=BAF_CN4, BAF_Normal=BAF_Normal, BAF_Del=BAF_Del, BAF_Dup=BAF_Dup, BAFs=BAFs)
	return(tmp)
}
