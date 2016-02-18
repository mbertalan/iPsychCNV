##' Verify position - checks if a genomic position is valid
##'
##' @title Verify position
##' @param Pos Position of the loci to plot in the form chr21:1050000-1350000
##' @param argument A character string, default being position, which is
##' added to the beining of the error message. E.g position/highlight
##' argument is not valid.
##' @return An error message if the postion is not valid
##' @author Johan Hilge Thygesen
##' @examples
##' VerifyPos("chr21:1050000-135")
##' @export
VerifyPos <- function(Pos, argument = "Position"){
    regmatch <- regexpr("chr[0-9X]{1,2}:[0-9]+-[0-9]+", Pos) # Regexp test of position
    if(attr(regmatch, "match.length") == nchar(Pos)) { # Split postion
        chr <- sub("chr", "", unlist(strsplit(Pos,":"))[1])
        Pos <- unlist(strsplit(Pos,":"))[2]
        x.start <- as.numeric(unlist(strsplit(Pos,"-"))[1])
        x.stop <- as.numeric(unlist(strsplit(Pos,"-"))[2])
        ## Test if start is smaller than stop
        if(x.stop<=x.start) {  
            stop(paste(argument,"stop is <= start /n/n"))
        }
    }else{
        stop(paste("the",argument,"position argument is not valid see --help"))
    }
    return(c(chr,x.start,x.stop))
}
