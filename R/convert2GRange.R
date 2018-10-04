# convert2GRange
#' @title convert output of BigLD/GPART to "GRangesList" object
#' @name convert2GRange
#' @aliases convert2GRange
#' @description \code{convert2GRange} convert a BigLD or GPART output to a data of \code{GRangesList} object.
#' @usage convert2GRange(blockresult)
#' @param blockresult BigLD or GPART output
# <output>
#' @return \code{GRangeList} object including the BigLD or GPART output
#'
#' @examples
#'
#' testBigLD <- BigLD(geno=geno[,1:200], SNPinfo=SNPinfo[1:200,])
#' testBigLD_grange <- convert2GRange(testBigLD)
#'
#' @author Sun-Ah Kim <sunny03@snu.ac.kr>, Yun Joo Yoo <yyoo@snu.ac.kr>
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @export
convert2GRange = function(blockresult){
  BigLDcolname = c("chr", "start.index", "end.index", "start.rsID",
                   "end.rsID", "start.bp", "end.bp")
  GPARTcolname = c("chr", "start.index", "end.index", "start.rsID",
                   "end.rsID", "start.bp", "end.bp", "blocksize", "Name")
  if(ncol(blockresult) == 7 & all(colnames(blockresult)[seq_len(7)] == BigLDcolname)){
    blocktype <- "bigld"
  }else if(ncol(blockresult) ==9 & all(colnames(blockresult) == GPARTcolname)){
    blocktype <- "gpart"
  }else{
    stop("The input data is not an output of BigLD or GPART functions")
  }
  if(blocktype == "bigld"){
    rowRanges <- GenomicRanges::GRanges(seqnames = blockresult[,1],
                                        IRanges::IRanges(blockresult[,6], width=(blockresult[,7]-blockresult[,6]+1)),
                                        start_rsID = as.character(blockresult[,4]),
                                        end_rsID = as.character(blockresult[,5]))
  }else if(blocktype == "gpart"){
    rowRanges <- GenomicRanges::GRanges(seqnames = blockresult[,1],
                                        IRanges::IRanges(blockresult[,6], width=(blockresult[,7]-blockresult[,6]+1)),
                                        start_rsID = as.character(blockresult[,4]),
                                        end_rsID = as.character(blockresult[,5]),
                                        blocksize = blockresult[,8],
                                        blockname = as.character(blockresult[,9]))
  }
  return(rowRanges)
}
