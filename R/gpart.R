#' @useDynLib gpart, .registration = TRUE
#'

VCFtogeno <- function(vcf) {
  .Call(`_gpart_VCFtogeno`, vcf)
}

pairCubeX <- function(b1, b2) {
  .Call(`_gpart_pairCubeX`, b1, b2)
}

matCubeX <- function(geno) {
  .Call(`_gpart_matCubeX`, geno)
}

matCubeX2 <- function(geno1, geno2) {
  .Call(`_gpart_matCubeX2`, geno1, geno2)
}

estiD <- function(pA, pB, n11, n12, n21, n22, n1212) {
  .Call(`_gpart_estiD`, pA, pB, n11, n12, n21, n22, n1212)
}

CIDp <- function(b1, b2) {
  .Call(`_gpart_CIDp`, b1, b2)
}

genoDp <- function(geno, strLD = TRUE, lower = 0.7, upper = 0.98) {
  .Call(`_gpart_genoDp`, geno, strLD, lower, upper)
}

genoCubeDp <- function(geno) {
  .Call(`_gpart_genoCubeDp`, geno)
}

genoDp2 <- function(geno1, geno2, strLD = TRUE, lower = 0.7, upper = 0.98) {
  .Call(`_gpart_genoDp2`, geno1, geno2, strLD, lower, upper)
}

calc_lnlike <- function(known11, known12, known21, known22, center_ct_d, freq11, freq12, freq21, freq22, half_hethet_share, freq11_incr) {
  .Call(`_gpart_calc_lnlike`, known11, known12, known21, known22, center_ct_d, freq11, freq12, freq21, freq22, half_hethet_share, freq11_incr)
}

cubic_real_roots <- function(coef_a, coef_b, coef_c, solutions) {
  .Call(`_gpart_cubic_real_roots`, coef_a, coef_b, coef_c, solutions)
}

calc_lnlike_quantile <- function(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, quantile) {
  .Call(`_gpart_calc_lnlike_quantile`, known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, quantile)
}

em_phase_hethet <- function(known11, known12, known21, known22, center_ct, onside_sol_ct_ptr) {
  .Call(`_gpart_em_phase_hethet`, known11, known12, known21, known22, center_ct, onside_sol_ct_ptr)
}

haploview_blocks_classify <- function(counts, lowci_max, lowci_min, recomb_highci, strong_highci, strong_lowci, strong_lowci_outer, is_x, recomb_fast_ln_thresh) {
  .Call(`_gpart_haploview_blocks_classify`, counts, lowci_max, lowci_min, recomb_highci, strong_highci, strong_lowci, strong_lowci_outer, is_x, recomb_fast_ln_thresh)
}

CIDp_strLD <- function(b1, b2, lower, upper) {
  .Call(`_gpart_CIDp_strLD`, b1, b2, lower, upper)
}


################# CLQD subFUNCTIONS ################
# subfunctions : 1.newSplitCliques 2. heuristicCLQ
#1
newSplitCliques <- function(cliques.bp, gapdist)
{
  nowlist <- lapply(cliques.bp, sort)
  fixlist <- NULL
  repeat {
    need.split <- which(lapply(nowlist, function(x) max(diff(x)) > gapdist) == TRUE)
    need.fix <- which(lapply(nowlist, function(x) max(diff(x)) > gapdist) == FALSE)
    addlist <- nowlist[need.fix]
    fixlist <- c(fixlist, addlist)
    if (length(need.split) == 0) {
      break
    }
    nowlist <- nowlist[need.split]
    nowlength <- length(nowlist)
    newlist <- as.list(rep(NA, nowlength))
    for (i in seq_len(nowlength)) {
      gap <- diff(nowlist[[i]])
      frontpart <- nowlist[[i]][seq_len(min(which(gap > gapdist)))]
      restpart <- nowlist[[i]][-(seq_len(min(which(gap > gapdist))))]
      nowlist[[i]] <- frontpart
      newlist[[i]] <- restpart
    }
    addlist <- nowlist[vapply(nowlist, function(x) length(x) > 1, TRUE)]
    fixlist <- c(fixlist, addlist)
    nowlist <- newlist[vapply(newlist, function(x) length(x) > 1, TRUE)]
  }
  return(fixlist)
}
#2
heuristicCLQ <- function(subOCM, hrstParam)
{
  degs <- apply(subOCM, 1, sum)
  top5pctDegVs <- which(degs>=quantile(degs, 0.95))
  top5pctDegVs <- top5pctDegVs[order(degs[top5pctDegVs],decreasing=TRUE)]
  heuristicBinbasket <- NULL
  total_density <- NULL
  v_dens_all <- c()
  find_opt <- FALSE
  for(v in top5pctDegVs){
    v_nbds <- as.integer(which(subOCM[v,]!=0))
    candibin <- c(v, v_nbds)
    v_density <- sum(subOCM[candibin,candibin])/(length(candibin)^2 - length(candibin))
    v_dens <- c(v_density, candibin)
    if(v_density>=0.95) {
      return(candibin)
      find_opt <- TRUE
      break
    }else{
      new_v_density <- 1
      while((v_density+0.0001)<new_v_density){
        newdeg <- apply(subOCM[v_nbds,], 1, sum)
        v_nbds <- v_nbds[-(which(newdeg==min(newdeg)))]
        candibin <- c(v, v_nbds)
        new_v_density <- sum(subOCM[candibin,candibin])/(length(candibin)^2 - length(candibin))
        if(new_v_density>=0.95){
          return(candibin)
        }else{
          v_dens_all <- c(v_dens_all, new_v_density)
          heuristicBinbasket<- c(heuristicBinbasket, list(candibin))
        }
        if(length(v_nbds) < (2*hrstParam/3)) break
      }#end while
    }#end else
  }#end for
  ## return results
  max_den_bin = heuristicBinbasket[which(v_dens_all == max(v_dens_all))]
  if(length(max_den_bin) == 1) {
    return(max_den_bin[[1]])
  }else{
    max_den_bin_leng = vapply(max_den_bin, length, 1)
    candibin = max_den_bin[[which(max_den_bin_leng == max(max_den_bin_leng))[1]]]
    return(candibin)
  }
}


################# BigLD subFUNCTIONS ##########################
# subfunctions : CLQD,
# dataPreparation , cutsequence, intervalCliqueList,
# findMaximumIndept, constructLDblock, subBigLD, appendcutByForce
dataPreparation <- function(genofile, SNPinfofile, geno, SNPinfo,
                            chrN = NULL , startbp = -Inf, endbp = Inf)
{
  #chrN : vector of chromosome
  #startbp, endbp : integer
  # filetype <- txt, the input data have header
  # if(is.null(genofile)) stop("Please input genofile (and SNPinfofile)")
  filetype <- tail(strsplit(genofile, "\\.")[[1]], 1)
  if(filetype == "txt"){
    if(is.null(SNPinfo)){
      SNPinfo <- data.table::fread(SNPinfofile, header=TRUE)
    }
    class(SNPinfo)<-"data.frame"
    if(dim(SNPinfo)[2]==4){
      SNPinfo <- SNPinfo[,c(1,3,4)] # chrN, rsID, bp
    }
    geno <- data.table::fread(genofile, header = TRUE) # rsID
    geno <- as.matrix(geno)
    # header generation
    if(is.null(chrN)) chrN <- unique(SNPinfo[,1])
    colnames(SNPinfo) <- c("chrN", "rsID", "bp")
    if(startbp != -Inf | endbp != Inf){
      if((length(chrN)>1) | (length(unique(SNPinfo[,1]))>1))
        stop("If your input data include more than one chromosome or you choose more than one chromosome,
             then do not specify the startbp and endbp!")
      SNPloca <- which(SNPinfo$bp>=startbp & SNPinfo$bp<=endbp & SNPinfo$chrN == chrN)
      SNPinfo <- SNPinfo[SNPloca,]
      geno <- geno[,SNPloca]
    }
    colnames(geno)<- SNPinfo$rsID
    if(dim(geno)[2] != dim(SNPinfo)[1]) stop("Please check the input file format!")
  }else if(filetype == "raw"){
    if(is.null(SNPinfofile) & is.null(SNPinfo)) stop("Please input map file for SNPinfofile")
    if(!is.null(SNPinfofile)){
      if((tail(strsplit(SNPinfofile, "\\.")[[1]],1)!="map")) stop("Please input map file for SNPinfofile")
      SNPinfo <- data.table::fread(SNPinfofile, header=FALSE, sep = "\t")
      class(SNPinfo)<-"data.frame"
      SNPinfo <- SNPinfo[,c(1,2,4)]
      colnames(SNPinfo) <- c("chrN", "rsID", "bp")
      geno <- data.table::fread(genofile, header = TRUE, sep = " ")
      geno <- geno[, -seq_len(6)]
      geno <- as.matrix(geno)
      if(is.null(chrN)) chrN <- unique(SNPinfo[,1])
      if(startbp != -Inf | endbp != Inf){
        if((length(chrN)>1) | (length(unique(SNPinfo[,1]))>1))
          stop("If your input data include more than one chromosome or you choose more than one chromosome,
               then do not specify the startbp and endbp!")
        SNPloca <- which(SNPinfo$bp>=startbp & SNPinfo$bp<=endbp & SNPinfo$chrN == chrN)
        SNPinfo <- SNPinfo[SNPloca,]
        geno <- geno[,SNPloca]
      }
      colnames(geno)<- SNPinfo$rsID
    }else{
      stop("Please check your input files")
    }
  }else if(filetype == "traw"){
    genodata <- data.table::fread(genofile, header = TRUE, sep = "\t")
    SNPinfo <- genodata[,c(1,2,4)]
    colnames(SNPinfo) <- c("chrN", "rsID", "bp")
    class(SNPinfo)<-"data.frame"
    geno <- t(genodata[,-seq_len(6)])
    geno <- as.matrix(geno)
    colnames(geno) <- SNPinfo$rsID
    if(is.null(chrN)) chrN <- unique(SNPinfo[,1])
    if(startbp != -Inf | endbp != Inf){
      if((length(chrN)>1) | (length(unique(SNPinfo[,1]))>1))
        stop("If your input data include more than one chromosome or you choose more than one chromosome,
             then do not specify the startbp and endbp!")
      takenIndex <- which((SNPinfo$bp >= startbp) & (SNPinfo$bp <= endbp) & (SNPinfo$chrN == chrN))
      SNPinfo <- SNPinfo[takenIndex,]
      geno <- geno[, takenIndex]
    }
  }else if(filetype == "ped"){
    if(is.null(SNPinfofile) & is.null(SNPinfo)) stop("Please input map file for SNPinfofile")
    if(!is.null(SNPinfofile)){
      # system.time({
      if((tail(strsplit(SNPinfofile, "\\.")[[1]],1)!="map")) stop("Please input map file for SNPinfofile")
      SNPinfo <-  data.table::fread(SNPinfofile, header=FALSE, sep = "\t")
      SNPinfo <- SNPinfo[,c(1,2,4)]
      colnames(SNPinfo) <- c("chrN", "rsID", "bp")
      class(SNPinfo)<-"data.frame"
      geno <- data.table::fread(genofile, header=FALSE, sep = " ", colClasses = "character")
      geno <- geno[,-seq_len(6)]
      geno <- as.matrix(geno)
      # message(dim(geno))
      if(is.null(chrN)) chrN <- unique(SNPinfo[,1])
      if(startbp != -Inf | endbp != Inf){
        if((length(chrN)>1) | (length(unique(SNPinfo[,1]))>1))
          stop("If your input data include more than one chromosome or you choose more than one chromosome,
               then do not specify the startbp and endbp!")
        SNPloca <- which(SNPinfo$bp >= startbp & SNPinfo$bp <= endbp & SNPinfo$chrN == chrN)
        SNPinfo <- SNPinfo[SNPloca,]
        genoloca1 <- min(SNPloca) * 2 - 1
        genoloca2 <- max(SNPloca) * 2
        geno <- geno[, genoloca1:genoloca2]
      }
      multiAllele<-c()
      nInd = nrow(geno)
      for(i in seq_len(ncol(geno)/2)){
        nowseq <- c(geno[, (i*2-1)], geno[, (i*2)])
        naInd = unique(c(which(geno[, (i*2-1)] == "0"), which(geno[, (i*2)] == "0")))
        uniqueAllele <- setdiff(unique(nowseq), "0")
        if(length(uniqueAllele)>2) multiAllele <- c(multiAllele, i)
        line1 <- rep(NA, nInd)
        line2 <- rep(NA, nInd)
        line1[geno[,(i*2-1)] == uniqueAllele[1]] <- 0
        line1[geno[,(i*2-1)] == uniqueAllele[2]] <- 1
        line2[geno[,(i*2)] == uniqueAllele[1]] <- 0
        line2[geno[,(i*2)] == uniqueAllele[2]] <- 1
        nline = (line1+line2)
        if(sum(nline, na.rm = TRUE) > nInd) nline = 2-nline
        geno[,(i)]<- nline
      }
      geno<-geno[,seq_len(ncol(geno)/2)]
      # message(dim(geno))
      class(geno)<-"integer"
      # remove multi-allelic
      if(length(multiAllele)>0){
        SNPinfo <- SNPinfo[-multiAllele,]
        geno <- geno[, -multiAllele]
      }
      # })
    }
  }else if(filetype == "vcf"){
    rawgeno <- data.table::fread(genofile, header=TRUE, sep = "\t",
                                 colClasses = "character", skip = "#CHROM")
    SNPinfo <- rawgeno[,c(1,3,2)]
    colnames(SNPinfo) <- c("chrN", "rsID", "bp")
    class(SNPinfo)<-"data.frame"
    infocol <- max(which(colnames(rawgeno) == "INFO"),which(colnames(rawgeno) == "FORMAT"))
    rawgeno <- rawgeno[, -seq_len(infocol)]
    rawgeno <- as.matrix(rawgeno)
    geno <- VCFtogeno(rawgeno)
    if(is.null(chrN)) chrN <- unique(SNPinfo[,1])
    if(startbp != -Inf | endbp != Inf){
      if((length(chrN)>1) | (length(unique(SNPinfo[,1]))>1))
        stop("If your input data include more than one chromosome or you choose more than one chromosome,
             then do not specify the startbp and endbp!")
      SNPloca <- which(SNPinfo$bp >= startbp & SNPinfo$bp <= endbp & SNPinfo$chrN == chrN)
      SNPinfo <- SNPinfo[SNPloca,]
      geno <- geno[,SNPloca]
    }
  }else {
    stop("Please check your input files.\n(supporting file format: txt, ped, map, raw, traw, vcf)")
  }
  #choose chrN
  message("data reformating done!")
  reverseind <- which(colSums(geno)>nrow(geno))
  geno[,reverseind] <- 2-geno[,reverseind]
  class(SNPinfo$bp)<-"numeric"
  class(SNPinfo$rsID) <- "character"
  #choose chrN
  index <- (SNPinfo$chrN %in% chrN)
  if(all(index == TRUE) == FALSE){
    geno <- geno[,index]
    SNPinfo <- SNPinfo[,index]
  }
  return(list(geno = geno, SNPinfo = SNPinfo))
}
cutsequence <- function(geno, SNPinfo, leng, subTaskSize, LD, clstgap, CLQcut, cutByForce)
{
  message("start to split the whole sequence into sub-task regions")
  modeNum <- 1
  lastnum <- 0
  nSNPs <- dim(SNPinfo)[1]
  geno <- as.matrix(geno)
  mincut <- max(0.5, CLQcut^2)
  # region length<=3000
  SNPbpgap <- diff(SNPinfo[,3])
  if(is.null(cutByForce)){
    precutpoints <- which(SNPbpgap>clstgap)
  }else{
    cutByForce <- vapply(cutByForce[,3], function(x) max(which(SNPinfo[,3]<=x)), 1)
    precutpoints <- sort(unique(cutByForce))
    precutpoints1 <- which(SNPbpgap>clstgap)
    precutpoints <- c(precutpoints, precutpoints1)
  }
  if (dim(geno)[2] <= subTaskSize) {
    message("there is only one sub-region!")
    return(list(dim(geno)[2], NULL))
  } else {
    calterms <- c(1, 10, leng)
    cutpoints <- NULL
    i <- leng # i :current candidate cutposition
    while (i <= (dim(geno)[2] - leng)) {
      if((i-lastnum) > 5*subTaskSize){
        modeNum <- 2
        break;
      }
      if(i %in% precutpoints){
        cutpoints <- c(cutpoints, i)
        lastnum <- i
        cat("total SNPs = ",nSNPs, " | cutpoint = ", i, "\r")
        i<-i+(leng/2)
        # message(i)
        cutnow <- FALSE
      }else{
        for(j in calterms){
          if(LD == "r2"){
            nowcm <- suppressWarnings(cor(geno[,(i-j+1):(i)], geno[,((i+1):(i+j))]
                                          ,use="pairwise.complete.obs"))
            nowr2 <- nowcm^2
            nowr2[which(nowr2 < mincut)] <- 0
          }else{
            nowr2 <- genoDp2(geno[,(i-j+1):(i), drop=FALSE], geno[,((i+1):(i+j)), drop=FALSE])
          }
          if(sum(nowr2,na.rm=TRUE)>0){
            i <- i+1
            cutnow <- FALSE
            break
          }
          if(j == leng){
            cutnow <- TRUE
            # message(i)
          }
        }
        if(cutnow == TRUE){
          cutpoints <- c(cutpoints, i)
          lastnum <- i
          cat("total SNPs = ",nSNPs, " | cutpoint = ", i, "\r")
          i<-i+(leng/2)
          cutnow <- FALSE
        }
      }
    }##end while
    if(modeNum == 1){
      cutpoints <- c(0,cutpoints, dim(geno)[2])
      # separate too big regions candi.cutpoints return(cutpoints,candi.cutpoints)
      atfcut <- NULL
      while (max(diff(cutpoints)) > subTaskSize) {
        diffseq <- diff(cutpoints)
        recutpoint <- which(diffseq > subTaskSize)
        nowmaxsize <- max(diff(cutpoints))
        tt <- cbind((cutpoints[recutpoint] + 1), cutpoints[recutpoint + 1])
        numvec <- NULL
        for (i in seq_len(dim(tt)[1])){
          st <- tt[i, 1]
          ed <- tt[i, 2]
          if (ed > (dim(geno)[2] - leng)) {
            ed <- dim(geno)[2] - leng
          }
          weakcount <- vapply(c((st + leng):(ed - leng)), function(x) {
            tick <- as.integer(leng/5)
            if(LD == "r2"){
              nowCM <- suppressWarnings(cor(geno[, (x - tick + 1):(x)], geno[, (x+ 1):(x + tick)]
                                            ,use="pairwise.complete.obs"))
              nowr2 <- nowCM^2
              length(which(nowr2>= mincut))
            }else{
              nowr2 <- genoDp2(geno[, (x - tick + 1):(x), drop=FALSE], geno[, (x+ 1):(x + tick), drop=FALSE])
              length(which(nowr2>= mincut))
            }
          }, 1)
          weakcount.s <- sort(weakcount)
          weaks <- weakcount.s[10]
          weakpoint <- which(weakcount <= weaks)
          weakpoint <- weakpoint + st + leng - 1
          nearcenter <- vapply(weakpoint, function(x) abs((ed - x) - (x - st)), 1)
          addcut <- weakpoint[which(nearcenter == min(nearcenter))][1]
          cat("total SNPs = ",nSNPs, " | additional cutpoint = ", addcut, "\r")
          numvec <- c(numvec, addcut)
          atfcut <- c(atfcut, addcut)
        }  ##end for
        cutpoints <- sort(c(cutpoints, numvec))
        newcandi <- which(diff(cutpoints) > subTaskSize)
        # remaxsize <- max(diff(cutpoints))
        # message(remaxsize) message(newcandi)
        if (length(newcandi) == 0) {
          break
        }
      }
    }
    ##end while
    if(modeNum == 2){
      message(paste("split the whole sequence into sub-task regions of length", subTaskSize))
      precutpoints <- which(SNPbpgap>clstgap)
      cutpoints <- seq(subTaskSize, dim(geno)[2], subTaskSize)
      remaincutpoints <- vapply(cutpoints, function(x) sum(abs(x-precutpoints)<(leng/2))<1, TRUE)
      atfcut <- cutpoints[remaincutpoints]
      cutpoints <- sort(c(cutpoints[remaincutpoints], precutpoints))
      if(max(atfcut) == dim(geno)[2]){
        atfcut <- atfcut[-(length(atfcut))]
      }else{
        cutpoints <- c(cutpoints, dim(geno)[2])
      }
    }
  }
  cat("\n")
  message("splitting sequence, done!")
  if(!is.null(atfcut)) atfcut <- sort(atfcut)
  cutpoints <- unique(c(cutpoints, precutpoints))
  return(list(sort(cutpoints), atfcut))
}
intervalCliqueList <- function(clstlist)
{
  clstlist <- clstlist[order(vapply(clstlist, min, 1))]
  bp.clstlist <- t(vapply(clstlist, function(x) range(x), c(1,2)))  ###
  bp.clstlist <- bp.clstlist[order(bp.clstlist[,1]),]
  IMsize <- dim(bp.clstlist)[1]  ## adjacency matrix of intervals in interval graph
  adjacencyM <- matrix(0, IMsize, IMsize)
  for (i in seq_len(IMsize-1)) {
    for (j in (i+1):IMsize) {
      if(bp.clstlist[i,2]>bp.clstlist[j,1]){
        adjacencyM[j, i]<-adjacencyM[i, j] <- 1
      }else{
        next
      }
    }
  }
  diag(adjacencyM) <- 0
  interval.graph <- igraph::graph.adjacency(adjacencyM, mode="undirected",
                                            weighted=TRUE, diag=FALSE, add.colnames=NULL)
  # message(paste("max coreness", max(coreness(interval.graph))))
  # message(paste("ecount", ecount(interval.graph), "vertex*5 ", 5*IMsize))
  if(max(igraph::coreness(interval.graph))>10){ #ecount(interval.graph)> 5*IMsize|
    interval.cliques <- igraph::maximal.cliques(interval.graph, min=1)
  }else{
    interval.cliques <- igraph::cliques(interval.graph, min=1)
  }
  interval.cliques <- interval.cliques[order(vapply(interval.cliques, min, 1))]

  intervals <- lapply(interval.cliques, function(x) unlist(clstlist[x]))
  intervals <- lapply(intervals, sort)
  intervals <- lapply(intervals, unique)
  weight.itv <- vapply(intervals, length, 1)

  intervals.range <- t(vapply(intervals, range, c(1,2)))
  unique.intervals.range <- unique(intervals.range)

  rangeinfo <- cbind(intervals.range, weight.itv)

  interval.info <- apply(unique.intervals.range, 1, function(x) {
    sameitv <- which(rangeinfo[, 1] == x[1] & rangeinfo[, 2] == x[2])
    maxweight <- max(rangeinfo[sameitv, 3])
    sameitv[which(maxweight == rangeinfo[sameitv, 3])][1]
  })
  res <- rangeinfo[interval.info,,drop=FALSE]
  return(res)
  # final.intervals <- intervals[interval.info]
  # final.intervals.w <- rangeinfo[interval.info, 3]
  # return(list(final.intervals, final.intervals.w))
}
findMaximumIndept <- function(interval.range, sample.weight)
{
  n.of.sample <- length(sample.weight)
  pre.range <- as.list(rep(NA, n.of.sample))  #pre.range : range of predecessor
  ## pre.range : n by 2, i row (x,y) : possible predecessors of i interval are from x interval to y interval
  for (i in seq_len(n.of.sample)) {
    nowstart <- interval.range[i, 1]
    if (sum(interval.range[, 2] < nowstart) > 0) {
      pre.range[[i]] <- which(interval.range[, 2] < nowstart)
    }
  }
  sources <- seq_len(n.of.sample)[(vapply(pre.range, function(x) all(is.na(x)) == TRUE, TRUE))]

  ## source of comparability graph of complement of Interval graph
  if (length(sources) < n.of.sample) {
    not.s <- setdiff(seq_len(n.of.sample), sources)
    for (i in not.s) {
      pre.pre <- sort(unique(unlist(pre.range[pre.range[[i]]])))
      pre.range[[i]] <- setdiff(pre.range[[i]], pre.pre)
    }
    names(pre.range) <- sample.weight
    n.interval <- seq_len(n.of.sample)
    route.weights <- rep(0, n.of.sample)  ##cumulative weights
    route.weights[sources] <- sample.weight[sources]
    pointers <- rep(0, n.of.sample)  ## predecessor of current interval
    pointers[sources] <- NA
    explored <- rep(0, n.of.sample)
    explored[sources] <- 1
    info <- cbind(n.interval, route.weights, pointers, explored)

    for (i in not.s) {
      maybe.pred <- pre.range[[i]]
      now.info <- info[maybe.pred, , drop=FALSE]
      max.info <- now.info[which(now.info[, 2] == max(now.info[, 2])), , drop=FALSE]
      if (dim(max.info)[1] > 1)
        max.info <- max.info[1, , drop=FALSE]
      info[i, 2] <- sample.weight[i] + max.info[2]
      info[i, 3] <- max.info[1]
      info[i, 4] <- 1
    }

    #### trace maximum independent set
    start.itv <- which(info[, 2] == max(info[, 2]))[1]
    predecessor <- info[start.itv, 3]
    indept.set <- c(predecessor, start.itv)
    while (!is.na(predecessor)) {
      predecessor <- info[predecessor, 3]
      indept.set <- c(predecessor, indept.set)
    }

    indept.set <- as.vector(indept.set)
    indept.set <- indept.set[-which(is.na(indept.set))]
    indept.set.weight <- max(info[, 2])
  } else {
    indept.set <- which(sample.weight == max(sample.weight))
    indept.set.weight <- max(sample.weight)
  }

  final.result <- list(indept.set=indept.set, indept.set.weight=indept.set.weight)
  return(final.result)
}
constructLDblock <- function(clstlist, subSNPinfo)
{
  # subfunction: intervalCliqueList, findMaximumIndept
  Totalblocks <- NULL
  while (length(clstlist) > 0) {
    if(length(clstlist)==1){
      Totalblocks <- rbind(Totalblocks, range(clstlist[[1]]))
      break
    }else{
      allsnps <- lapply(clstlist, function(x) c(min(x):max(x)))
      candi.interval <- intervalCliqueList(clstlist)
      intervals <- candi.interval[,seq_len(2),drop=FALSE]  ## list of SNPs in each cliques
      weight.itv <- candi.interval[,3]  ## weights of each cliques
      MWIS <- findMaximumIndept(intervals, weight.itv)  ##find independent set
      subLDblocks <- intervals[MWIS[[1]],, drop=FALSE]
      Totalblocks <- rbind(Totalblocks, subLDblocks)
      takenSNPs <- apply(subLDblocks, 1, function(x) c(min(x):max(x)))
      takenSNPs <- as.integer(unlist(takenSNPs))
      clstlist <- lapply(clstlist, function(x) setdiff(x, takenSNPs))
      clstlist <- clstlist[vapply(clstlist, length, 1)>0]
      if (length(clstlist) == 0) break
      # when a taken cluster located in the middle of the other cluster, we split the other cluster
      while(length(clstlist)>0){
        addinglist <- NULL
        for (n in seq_len(length(clstlist))) {
          nowbin <- clstlist[[n]]
          intersection <- intersect(c(min(nowbin):max(nowbin)), takenSNPs)
          if (length(intersection) > 0) {
            part1 <- nowbin[which(nowbin < min(intersection))]
            part2 <- setdiff(nowbin, c(min(part1):max(intersection)))
            clstlist[[n]] <- part1
            addinglist <- c(addinglist, list(part2))
          }
        }
        clstlist <- c(clstlist, addinglist)
        clstlist <- clstlist[vapply(clstlist, function(x) length(x) > 1, TRUE)]
        clstlist <- clstlist[order(vapply(clstlist, min, 1))]
        if(length(addinglist) ==0) break;
      }
    }

  }
  return(Totalblocks)
}
subBigLD <- function(subgeno, subSNPinfo,  CLQcut, clstgap, CLQmode, hrstParam, LD, hrstType)
{
  subbinvec <- CLQD(geno=subgeno, SNPinfo=subSNPinfo, CLQcut=CLQcut, clstgap=clstgap,
                    hrstType=hrstType, hrstParam=hrstParam, CLQmode=CLQmode, LD=LD)
  if(all(is.na(subbinvec))){return(NULL)}
  bins <- seq_len(max(subbinvec[which(!is.na(subbinvec))]))
  clstlist <- lapply(bins, function(x) which(subbinvec == x))
  clstlist <- lapply(clstlist, sort)  ###
  clstlist <- clstlist[order(vapply(clstlist, min, 1))]  ###
  nowLDblocks <- constructLDblock(clstlist, subSNPinfo)
  # message('constructLDblock done!')
  nowLDblocks <- nowLDblocks[order(nowLDblocks[, 1]), , drop=FALSE]
  return(nowLDblocks)
}

appendcutByForce <- function(LDblocks, Ogeno, OSNPinfo, CLQcut, clstgap,
                             CLQmode, hrstParam, LD, hrstType)
{
  expandB <- NULL
  snp1 <- which(OSNPinfo[,3]<LDblocks[1,5])
  if(length(snp1)>2){
    OSNPs <- seq_len(max(snp1))
    firstB <- LDblocks[1,]
    secondSNPs <- (as.numeric(firstB[1]):as.numeric(firstB[2]))
    if(LD == "r2"){
      cor2 <- suppressWarnings(cor(Ogeno[,OSNPs,drop=FALSE], Ogeno[,secondSNPs,drop=FALSE],
                                   use="pairwise.complete.obs")^2)
    }else if(LD == "Dprime"){
      cor2 <- genoDp2(Ogeno[,OSNPs,drop=FALSE], Ogeno[,secondSNPs,drop=FALSE])
    }
    cor2num <- apply(cor2, 1, function(x) {
      sum(x>CLQcut^2, na.rm = TRUE)
    })
    cor2ratio <- cor2num/(dim(cor2)[2])
    cor2numT <- cor2ratio>0.6
    cor2numT <- c(cor2numT, 1)
    points2 <- min(which(cor2numT>0))
    NsecondSNPs <- points2:max(secondSNPs)
    reOSNPs <- setdiff(seq_len(max(NsecondSNPs)), NsecondSNPs)
    if(length(reOSNPs)>1){
      subgeno <- Ogeno[, reOSNPs]
      subSNPinfo <- OSNPinfo[reOSNPs,]
      subBlocks <- subBigLD(subgeno, subSNPinfo,  CLQcut, clstgap, CLQmode, hrstParam, LD, hrstType)
      subBlocks <- subBlocks+min(reOSNPs)-1
      expandB <- rbind(expandB, subBlocks)
    }
    firstSNPs <- NsecondSNPs
  }else{
    firstB <- LDblocks[1,]
    firstSNPs <- (as.numeric(firstB[1]):as.numeric(firstB[2]))
  }
  if(dim(LDblocks)[1]>1){
    for(i in seq_len(dim(LDblocks)[1]-1)){
      # if(i==2192) stop("stop")
      secondB <- LDblocks[(i+1),]
      secondSNPs <- (as.numeric(secondB[1]):as.numeric(secondB[2]))
      OSNPs <- setdiff(max(firstSNPs):min(secondSNPs), c(max(firstSNPs), min(secondSNPs)))
      if(length(OSNPs)==0){
        expandB <- rbind(expandB, range(firstSNPs))
        firstSNPs <- secondSNPs
      }else{
        if(LD == "r2"){
          cor1 <- suppressWarnings(cor(Ogeno[,firstSNPs,drop=FALSE], Ogeno[,OSNPs,drop=FALSE],
                                       use="pairwise.complete.obs")^2)
        }else if(LD == "Dprime"){
          cor1 <- genoDp2(Ogeno[,firstSNPs,drop=FALSE], Ogeno[,OSNPs,drop=FALSE])
        }
        cor1num <- apply(cor1, 2, function(x) {
          sum(x>CLQcut^2)
        })
        cor1ratio <- cor1num/(dim(cor1)[1])
        # cor1num <- apply(cor1r2, 2, function(x) length(which(x>CLQcut^2))/length(firstSNPs))
        if(LD == "r2"){
          cor2 <- suppressWarnings(cor(Ogeno[,OSNPs,drop=FALSE], Ogeno[,secondSNPs,drop=FALSE],
                                       use="pairwise.complete.obs")^2)
        }else if (LD == "Dprime"){
          cor2 <- genoDp2(Ogeno[,OSNPs,drop=FALSE], Ogeno[,secondSNPs,drop=FALSE])
        }
        cor2num <- apply(cor2, 1, function(x) {
          sum(x>CLQcut^2)
        })
        cor2ratio <- cor2num/(dim(cor2)[2])
        cor1numT <- cor1ratio>0.6
        cor2numT <- cor2ratio>0.6
        cor1numT <- c(1, cor1numT, 0)
        cor2numT <- c(0, cor2numT, 1)
        points1 <- max(firstSNPs)+max(which(cor1numT>0))-1
        NfirstSNPs <- min(firstSNPs):points1
        points2 <- max(firstSNPs)+max(which(cor2numT>0))-1
        NsecondSNPs <- points2:max(secondSNPs)
        if(max(NfirstSNPs)<min(NsecondSNPs)){
          expandB <- rbind(expandB, range(NfirstSNPs))
          reOSNPs <- setdiff(c(min(NfirstSNPs):max(NsecondSNPs)), c(NfirstSNPs, NsecondSNPs))
          if(length(reOSNPs)>1){
            subgeno <- Ogeno[, reOSNPs]
            subSNPinfo <- OSNPinfo[reOSNPs,]
            subBlocks <- subBigLD(subgeno, subSNPinfo, CLQcut, clstgap, CLQmode, hrstParam, LD, hrstType)
            subBlocks <- subBlocks+min(reOSNPs)-1
            expandB <- rbind(expandB, subBlocks)
          }
          firstSNPs <- NsecondSNPs
        }else{
          #merge two blocks
          # message(paste("i", i,"now!!!!!!"))
          subgeno <- Ogeno[, c(min(firstSNPs):max(secondSNPs))]
          subSNPinfo <- OSNPinfo[c(min(firstSNPs):max(secondSNPs)),]
          if(dim(subgeno)[2]<=3000){
            subBlocks <- subBigLD(subgeno, subSNPinfo,  CLQcut, clstgap, CLQmode, hrstParam, LD, hrstType)
            subBlocks <- subBlocks+min(firstSNPs)-1
          }else{
            # Too big LD blocks
            subgeno <- Ogeno[, c((max(firstSNPs)-1500):(min(secondSNPs)+1500))]
            subSNPinfo <- Ogeno[ c((max(firstSNPs)-1500):(min(secondSNPs)+1500)),]
            subBlocks <- subBigLD(subgeno, subSNPinfo,  CLQcut, clstgap, CLQmode, hrstParam, LD, hrstType)
            for(i in seq_len(dim(subBlocks)[1]-1)){
              nowB <- subBlocks[i,]
              nextB <- subBlocks[(i+1),]
              # merging steps!!
              if(nowB[2]<nextB[1]){
                next
              }else{
                subBlocks[i,]<-c(NA, NA)
                subBlocks[(i+1),]<-range(c(nowB, nextB))
              }
            }
            subBlocks<-subBlocks[!is.na(subBlocks[,1]),]
            subBlocks <- subBlocks+min(firstSNPs)-1
          }
          if(dim(subBlocks)[1]==1) {
            firstSNPs <- subBlocks[1,1]:subBlocks[1,2]
          }else{
            expandB <- rbind(expandB, subBlocks[-(dim(subBlocks)[1]),])
            firstSNPs <- subBlocks[(dim(subBlocks)[1]),1]:subBlocks[(dim(subBlocks)[1]),2]
          }
        }
        cat(paste(i," | ", dim(LDblocks)[1], "\r", sep = ""))
        # message(tail(expandB))
        # if(i >= 30) break
      }
    }
  }
  # firstSNPs
  if(max(firstSNPs)<(dim(Ogeno)[2]-1)){
    OSNPs <- (max(firstSNPs)+1):(dim(Ogeno)[2])
    if(LD == "r2"){
      cor1 <- suppressWarnings(cor(Ogeno[,firstSNPs,drop=FALSE], Ogeno[,OSNPs,drop=FALSE],
                                   use="pairwise.complete.obs")^2)
    }else if (LD == "Dprime"){
      cor1 <- genoDp2(Ogeno[,firstSNPs,drop=FALSE], Ogeno[,OSNPs,drop=FALSE])
    }
    cor1num <- apply(cor1, 2, function(x) {
      sum(x>CLQcut^2)
    })
    cor1ratio <- cor1num/(dim(cor1)[1])
    cor1numT <- cor1ratio>0.6
    cor1numT <- c(1, cor1numT, 0)
    points1 <- max(firstSNPs)+max(which(cor1numT>0))-1
    NfirstSNPs <- min(firstSNPs):points1
    expandB <- rbind(expandB, range(NfirstSNPs))
    reOSNPs <- setdiff(c(min(NfirstSNPs):dim(Ogeno)[2]), c(NfirstSNPs))
    if(length(reOSNPs)>1){
      subgeno <- Ogeno[, reOSNPs]
      subSNPinfo <- OSNPinfo[reOSNPs,]
      subBlocks <- subBigLD(subgeno, subSNPinfo,  CLQcut, clstgap, CLQmode, hrstParam, LD, hrstType)
      subBlocks <- subBlocks+min(reOSNPs)-1
      expandB <- rbind(expandB, subBlocks)
    }
  }else{
    expandB <- rbind(expandB, range(firstSNPs))
  }
  # LDblocks <- expandB
  expandB <- expandB[(expandB[,1]!=expandB[,2]),,drop=FALSE]
  start.bp <- OSNPinfo[, 3][expandB[, 1]]
  end.bp <- OSNPinfo[, 3][expandB[, 2]]
  start.rsID <- as.character(OSNPinfo[, 2][expandB[, 1]])
  end.rsID <- as.character(OSNPinfo[, 2][expandB[, 2]])
  TexpandB <- data.frame(expandB, start.rsID, end.rsID, start.bp, end.bp)
  colnames(TexpandB) <- c("start.index", "end.index", "start.rsID", "end.rsID", "start.bp", "end.bp")
  return(TexpandB)
}

#################  GPART subFUNCTIONS ########################################################################
# sub-Functions
# BigLD,CLQD, LDblockSplit, dataPreparation,
# mergeOverlapGene, LDblockGeneMerge, splitBigLD, mergeSmallRegion, namingRegion2
#
# blocks with Big-LD result
LDblockSplit <- function(geno, LDblocks, maxsize, LD){
  LDblockSizes <- apply(LDblocks, 1, diff) +1
  LargeBlocksN <- which(LDblockSizes> maxsize)
  LargeBlocks <- LDblocks[LargeBlocksN,,drop=FALSE]
  Newblocks <- NULL
  btwr2 <- NULL
  if(length(LargeBlocksN)!=0){
    while(dim(LargeBlocks)[1]>0){
      # message(dim(LargeBlocks))
      nowblocks <- LargeBlocks[1,]
      # if(nowblocks[1]>4000) break
      if(is.null(btwr2)){
        nowgeno <- geno[,nowblocks[1]:nowblocks[2]]
        if(LD == "r2"){
          nowr2 <- suppressWarnings(cor(nowgeno, use="pairwise.complete.obs")^2)
          nowr2[is.na(nowr2)] <- 0
        }else if(LD == "Dprime"){
          nowr2 <- genoDp(nowgeno, strLD = FALSE)
          nowr2[is.na(nowr2)] <- 0
          nowr2[(nowr2==Inf)|(nowr2==-Inf)] <- 0
        }
        btwr2 <- vapply(seq_len(dim(nowr2)[1]-1), function(x) mean(nowr2[seq_len(x), (x+1):(dim(nowr2)[1]-1)]), 1)
      }
      if(length(btwr2)< ceiling(maxsize*0.8*2)){
        subbtwr2 <- btwr2[(ceiling(length(btwr2)/2)-(maxsize*0.2)):(floor(length(btwr2)/2)+(maxsize*0.2))]
        addN <- which(subbtwr2 == min(subbtwr2))[1]+ (ceiling(length(btwr2)/2)-(maxsize*0.2)) - 1
        Newblocks <- rbind(Newblocks, c(nowblocks[1], (nowblocks[1]+addN-1)))
        if(Newblocks[dim(Newblocks)[1], 2]-Newblocks[dim(Newblocks)[1], 1]>=maxsize) {
          # message(Newblocks[dim(Newblocks)[1], ])
          break}
        Newblocks <- rbind(Newblocks, c(nowblocks[1]+addN, nowblocks[2]))
        if(Newblocks[dim(Newblocks)[1], 2]-Newblocks[dim(Newblocks)[1], 1]>=maxsize) {
          # message(Newblocks[dim(Newblocks)[1], ])
          break}
        LargeBlocks <- LargeBlocks[-1,, drop=FALSE]
        btwr2<- NULL
      }else{
        # message(LargeBlocks[1,]) Error in LargeBlocks[1, 1] <- (nowblocks[1] + addN)
        subbtwr2 <- btwr2[(maxsize*0.4):maxsize]
        addN <- which(subbtwr2 == min(subbtwr2))[1]+ (maxsize*0.4) - 1
        Newblocks <- rbind(Newblocks, c(nowblocks[1], (nowblocks[1]+addN-1)))
        LargeBlocks[1,1] <- (nowblocks[1]+addN)
        btwr2<-btwr2[-seq_len(addN)]
        if(diff(LargeBlocks[1,])<maxsize) {
          Newblocks <- rbind(Newblocks, LargeBlocks[1,])
          if(Newblocks[dim(Newblocks)[1], 2]-Newblocks[dim(Newblocks)[1], 1]>=maxsize) {
            # message(Newblocks[dim(Newblocks)[1], ])
            break}
          LargeBlocks <- LargeBlocks[-1,, drop=FALSE]
          btwr2<- NULL
        }
      }
      if(Newblocks[dim(Newblocks)[1], 2]-Newblocks[dim(Newblocks)[1], 1]>=maxsize) {
        # message(Newblocks[dim(Newblocks)[1], ])
        break}
    } # end while
    FinalLDblocks <- rbind(LDblocks[-LargeBlocksN,], Newblocks)
    FinalLDblocks <- FinalLDblocks[order(FinalLDblocks[,1]),]

    return(FinalLDblocks)
  }else{
    return(LDblocks)
  }
}
# merge overlapped gene regions
mergeOverlapGene <- function(Geneblocks) {
  overlaplist <- as.list(rep(NA, dim(Geneblocks)[1]))
  Geneblocks<-Geneblocks[order(Geneblocks[,1]),]
  for (i in seq_len(dim(Geneblocks)[1] - 1)) {
    overlaplist[[i]] <- i
    for (j in (i + 1):dim(Geneblocks)[1]) {
      if(Geneblocks[i,2]>=min(Geneblocks[j,1])){
        overlaplist[[i]] <- c(overlaplist[[i]], j)
      }else{
        next
      }
    }
  }
  overlaplist[[length(overlaplist)]]<- length(overlaplist)
  if(is.null(overlaplist)){
    return(Geneblocks)
  }else{
    for(i in seq_len(length(overlaplist)-1)){
      nowbin <- overlaplist[[i]]
      if(all(is.na(nowbin))) next
      for(j in (i+1):length(overlaplist)){
        if(all(overlaplist[[j]] %in% nowbin)){
          overlaplist[[j]]<-NA
        }else if(max(nowbin)<max(overlaplist[[j]])){
          break
        }
      }
    }
    overlaplist <- overlaplist[!vapply(overlaplist, function(x) all(is.na(x)), TRUE)]
    newblocks <- matrix(NA, nrow=length(overlaplist), ncol=2)
    rownames(newblocks)<-seq_len(length(overlaplist))
    for(i in seq_len(length(overlaplist))){
      if(length(overlaplist[[i]])==1){
        newblocks[i,]<-Geneblocks[overlaplist[[i]],,drop=FALSE]
        rownames(newblocks)[i]<-rownames(Geneblocks[overlaplist[[i]],,drop=FALSE])
      }else{
        newblocks[i,]<-c(min(Geneblocks[overlaplist[[i]],,drop=FALSE]),
                         max(Geneblocks[overlaplist[[i]],,drop=FALSE]))
        rownames(newblocks)[i]<-paste(rownames(Geneblocks[overlaplist[[i]],,drop=FALSE]), collapse="/")
      }

    }
    return(newblocks)
  }
}
# LDblock construction and split Large regions
LDblockGeneMerge <- function(LDblocks.T, Geneblocks.M) {
  LDblocks.T <- LDblocks.T[order(LDblocks.T[,1]),]
  Geneblocks.M <- Geneblocks.M[order(Geneblocks.M[,1]),]
  # for each gene blocks, we find the overlapping LD blocks
  LDwithGene.list <- matrix(NA, nrow=dim(Geneblocks.M)[1], ncol=2)
  rownames(LDwithGene.list) <-rownames(Geneblocks.M)
  overlapLDBlocks <- NULL
  for(i in seq_len(dim(Geneblocks.M)[1])){
    LD1 <- max(which(LDblocks.T[,1] <= Geneblocks.M[i,1]))
    LD2 <- min(which(LDblocks.T[,2] >= Geneblocks.M[i,2]))
    LDwithGene.list[i,] <- c(LD1, LD2)
    overlapLDBlocks <- c(overlapLDBlocks, LD1:LD2)
  }
  NonOverlapLDBlocks <- LDblocks.T[-overlapLDBlocks,]
  ## LDblocks which are not overlapped with gene regions

  overlapLDGene <- NULL
  NonoverlapLDgene <- NULL
  i <- 1
  while(i <= (dim(LDwithGene.list)[1])){
    nowGeneR <- LDwithGene.list[i,, drop=FALSE]
    if((LDwithGene.list[i,2] < LDwithGene.list[i+1,1])){
      NonoverlapLDgene <- rbind(NonoverlapLDgene, nowGeneR)
      # message(i);message(nowGeneR)
      i <- i+1
    }else{
      merge.candi <- nowGeneR
      for(j in (i+1):dim(LDwithGene.list)[1]){
        if(max(LDwithGene.list[i:(j-1),])>=LDwithGene.list[(j),1]){
          merge.candi<- rbind(merge.candi, LDwithGene.list[(j),,drop=FALSE])
          if(j == dim(LDwithGene.list)[1]){
            i <- j+1
          }
        }else{
          i <- j
          break
        }
      }
      overlapLDGene <- c(overlapLDGene, list(merge.candi))
    }
    if(i == dim(LDwithGene.list)[1]){
      nowGeneR <- LDwithGene.list[i,, drop=FALSE]
      NonoverlapLDgene <- rbind(NonoverlapLDgene, nowGeneR)
      # message(i);message(nowGeneR)
      break
    }
  }
  mergeLDgene <- NULL
  if(length(overlapLDGene)>0){
    mergeLDgene <- matrix(NA, length(overlapLDGene), 2)
    rownames(mergeLDgene) <- seq_len(length(overlapLDGene))
    for(i in seq_len(length(overlapLDGene))){
      nowmerge <- overlapLDGene[[i]]
      rownames(mergeLDgene)[i] <- paste(rownames(nowmerge), collapse="+")
      mergeLDgene[i,1] <- min(nowmerge)
      mergeLDgene[i,2] <- max(nowmerge)
    }
  }

  finalGeneBlocks <- rbind(mergeLDgene, NonoverlapLDgene)
  finalGeneBlocks <- finalGeneBlocks[order(finalGeneBlocks[,1]),]

  for(i in seq_len(dim(finalGeneBlocks)[1])){
    st <- LDblocks.T[finalGeneBlocks[i,1],1]
    ed <- LDblocks.T[finalGeneBlocks[i,2],2]
    finalGeneBlocks[i,] <-c(st, ed)
  }

  rownames(NonOverlapLDBlocks) <-rep("No", dim(NonOverlapLDBlocks)[1])
  res <- rbind(finalGeneBlocks, NonOverlapLDBlocks)
  res <- res[order(res[,1]),]
  return(res)
}
# split big Gene region
splitBigLD <- function(GeneLDblocks, LDblocks.T, Geneblocks,maxsize){
  BigN <- which(GeneLDblocks[,3]>maxsize)
  BigGeneblocks <- GeneLDblocks[BigN,,drop=FALSE]
  GeneLDblocks <- GeneLDblocks[-BigN,,drop=FALSE]
  newblocks <- NULL
  for (i in seq_len(dim(BigGeneblocks)[1])){
    nowBlock <- BigGeneblocks[i,,drop=FALSE]
    # nowSNPs<- nowBlock[1,1]:nowBlock[1,2]
    intersectLD <- LDblocks.T[max(which(nowBlock[1,1]>=LDblocks.T[,1])):min(which(nowBlock[1,2]<=LDblocks.T[,2])),]
    intersectLD[1,1] <- nowBlock[1,1]
    intersectLD[dim(intersectLD)[1],2] <- nowBlock[1,2]
    NintersectLD <- NULL
    while(dim(intersectLD)[1]>0){
      st <- intersectLD[1,1]
      edposi <- max(which(intersectLD[,2]<(st+maxsize)))
      ed <- intersectLD[edposi,2]
      NintersectLD <- rbind(NintersectLD, c(st, ed))
      intersectLD <- intersectLD[-seq_len(edposi),,drop=FALSE]
    }
    intersectLD <- NintersectLD
    # if(all(apply(NintersectLD, 1, diff)<50)==FALSE) message(i)
    rownames(intersectLD) <- rep(rownames(nowBlock), dim(intersectLD)[1])
    blockL <- apply(intersectLD, 1, diff)+1
    Addblocks <- cbind(intersectLD, blockL)
    newblocks <- rbind(newblocks, Addblocks)
  }
  GeneLDblocks <- rbind(GeneLDblocks, newblocks)
  GeneLDblocks <- GeneLDblocks[order(GeneLDblocks[,1]),]
  return(GeneLDblocks)
}
# small region merging
mergeSmallRegion <- function(GeneLDblocks, maxsize, minsize) {
  gname <- rownames(GeneLDblocks)
  GeneLDblocks <- as.data.frame(GeneLDblocks, stringsAsFactors=FALSE)
  # row.names(GeneLDblocks)<- NULL
  # GeneLDblocks<- cbind(GeneLDblocks, gname, stringsAsFactors<- F)
  smallbin <- NULL
  completeBin <- NULL
  while (dim(GeneLDblocks)[1] > 0) {
    # message(dim(GeneLDblocks))
    nowbin <- GeneLDblocks[1, ,drop=FALSE]
    # if(nowbin[1,1]>75530) break
    if (is.null(smallbin) & nowbin[1,3] >= minsize) {
      completeBin <- rbind(completeBin, nowbin)
      GeneLDblocks <- GeneLDblocks[-1, ,drop=FALSE]
    } else if (is.null(smallbin) & nowbin[1,3] < minsize) {
      smallbin <- nowbin
      GeneLDblocks<- GeneLDblocks[-1, ,drop=FALSE]
    } else if (smallbin[1, 3] + nowbin[1,3] < minsize) {
      st <- min(smallbin[1, 1], smallbin[1, 2], nowbin[1, 1], nowbin[1, 2])
      ed <- max(smallbin[1, 1], smallbin[1, 2], nowbin[1, 1], nowbin[1, 2])
      blockL <- smallbin[1, 3] + nowbin[1,3]
      smallbin <- data.frame(st, ed, blockL)
      GeneLDblocks <- GeneLDblocks[-1, ,drop=FALSE]
    } else if (smallbin[1, 3] + nowbin[1,3] >= minsize & smallbin[1, 3] + nowbin[1,3] <= maxsize) {
      st <- min(smallbin[1, 1], smallbin[1, 2], nowbin[1, 1], nowbin[1, 2])
      ed <- max(smallbin[1, 1], smallbin[1, 2], nowbin[1, 1], nowbin[1, 2])
      blockL <- smallbin[1, 3] + nowbin[1,3]
      completeBin <- rbind(completeBin, data.frame(st, ed, blockL))
      smallbin <- NULL
      GeneLDblocks <- GeneLDblocks[-1,,drop=FALSE]
    } else if (smallbin[1, 3] + nowbin[1,3] > maxsize) {
      lastbin <- completeBin[dim(completeBin)[1], ]
      if(is.null(lastbin)){
        message(c("Large-small-Large"))
        message(rbind(lastbin, smallbin, nowbin))
        completeBin <- rbind(completeBin, smallbin, nowbin)
        smallbin <- NULL
      }else if (lastbin[1,3] + smallbin[1,3] <= maxsize) {
        st <- min(smallbin[1, 1], smallbin[1, 2], lastbin[1, 1], lastbin[1, 2])
        ed <- max(smallbin[1, 1], smallbin[1, 2], lastbin[1, 1], lastbin[1, 2])
        blockL <- smallbin[1, 3] + lastbin[1, 3]
        completeBin <- completeBin[-dim(completeBin)[1], ,drop=FALSE]
        completeBin <- rbind(completeBin, data.frame(st, ed, blockL))
        smallbin <- NULL
      } else {
        message(c("Large-small-Large"))
        message(rbind(lastbin, smallbin, nowbin))
        completeBin <- rbind(completeBin, smallbin, nowbin)
        smallbin <- NULL
        GeneLDblocks <- GeneLDblocks[-1, ,drop=FALSE]
      }
    }
  }
  if(!is.null(smallbin) & is.null(completeBin)){
    completeBin < -smallbin
  }
  if(!is.null(smallbin) & !is.null(completeBin)){
    lastindex <- dim(completeBin)[1]
    if(completeBin[lastindex,3,drop=FALSE]+smallbin[1,3] > maxsize){
      completeBin <- rbind(completeBin, smallbin)
    }else{
      nowbin <- c(completeBin[lastindex, 1], smallbin[1,2], completeBin[lastindex, 3] + smallbin[1,3])
      completeBin[lastindex,seq_len(3)] <- nowbin
    }
  }
  return(completeBin)
}
# naming each region
namingRegion2 <- function(GeneLDblocks, genelist, chrSNPinfo) {
  st.bp <- chrSNPinfo[GeneLDblocks[,1],3]
  ed.bp <- chrSNPinfo[GeneLDblocks[,2],3]
  bpBlocks <- cbind(GeneLDblocks[,seq_len(3), drop=FALSE], st.bp, ed.bp)
  FinalGeneLDblocks <- NULL
  PartN <- 1
  Pgene <- NULL
  gene <- NULL
  Ngene <- NULL
  gname <- rep(NA, dim(bpBlocks)[1])
  FinalGeneLDblocks<- data.frame(bpBlocks, gname)

  for (i in seq_len(dim(bpBlocks)[1])){
    # if (i%%10 == 0)
    # message(i)
    nowloca <- bpBlocks[i,4:5,drop=FALSE]
    intro.index <- (which(genelist[,4]<nowloca[1,1]))
    intro.index1 <- setdiff(seq_len(dim(genelist)[1]), intro.index)
    outro.index <- (which(genelist[,3]>nowloca[1,2]))
    outro.index1 <- setdiff(seq_len(dim(genelist)[1]), outro.index)
    common.index <- intersect(intro.index1, outro.index1)

    if(length(common.index)==0){
      if(length(intro.index)==0){
        #before gene
        gname <- paste("before-", as.character(genelist[1, 1]), sep="")
      }else if(length(outro.index)==0){
        #after gene
        gname <- paste("after-", as.character(genelist[dim(genelist)[1], 1]), sep="")
      }else{
        intronum <- which(genelist[,4]==max(genelist[intro.index,4]))
        outronum <- which(genelist[,3]==min(genelist[outro.index,3]))
        gname <- paste("inter-", paste(as.character(genelist[intronum, 1]),collapse = "/"),":",paste(as.character(genelist[outronum, 1]), collapse = "/"), sep="")
      }
    }else{
      subgenelist <- genelist[common.index,,drop=FALSE]
      gname <- as.character(subgenelist[1,1])
      nowregion <- c(subgenelist[1,3], subgenelist[1,4])
      if(length(common.index)>1){
        for(k in 2:dim(subgenelist)[1]){
          nextregion <- subgenelist[k,3:4]
          if(max(nowregion)<min(nextregion)){
            gname <- paste(c(gname, as.character(subgenelist[k,1])), collapse="+")
            nowregion <- subgenelist[k,3:4]
          }else{
            gname <- paste(c(gname, as.character(subgenelist[k,1])), collapse="/")
            candi.r <- rbind(nowregion, nextregion)
            nowregion <- c(min(candi.r)[1], max(candi.r)[1])
          }
        }
      }
    }
    FinalGeneLDblocks$gname[i]<-gname
  }
  # part numbering
  FinalGeneLDblocks <- data.frame(FinalGeneLDblocks, 0)
  FinalGeneLDblocks[1, 7] <- 1
  for (i in 2:dim(FinalGeneLDblocks)[1]) {
    ifelse(FinalGeneLDblocks[i - 1, 6] == FinalGeneLDblocks[i, 6],
           FinalGeneLDblocks[i, 7] <- FinalGeneLDblocks[(i - 1), 7] + 1, FinalGeneLDblocks[i, 7] <- 1)
    # if(i%%10==0) message(i)
  }

  Finalnames <- apply(FinalGeneLDblocks, 1, function(x) {
    # paste(as.character(x[6]), '-part', as.character(x[7]), sep='')
    paste(c(x[6], "-part", as.numeric(x[7])), collapse="")
  })
  FinalGeneLDblocks <- cbind(FinalGeneLDblocks[, seq_len(5)], Finalnames)
  return(FinalGeneLDblocks)
}

geneinfoExt <- function(geneinfofile=geneinfofile, geneDB = geneDB, assembly = assembly, geneid = geneid,
                         ensbversion = ensbversion, chrNs = chrNs){
  if(!is.null(geneinfofile)){
    message("load gene information from inputed file")
    geneinfo <- data.table::fread(geneinfofile, stringsAsFactors = FALSE)
    geneinfo <- data.frame(geneinfo)
    colnames(geneinfo) <- c("genename", "chr", "start.bp", "end.bp")
  }else if(geneDB == "ensembl"){
    if(assembly=="GRCh38"){
      message("load gene information from ensembl (assembly GRCh38)")
      grch <- 38
      ensembl <- biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=ensbversion)#, version=version
      if(geneid!="ensembl_gene_id"){
        geneinfo <- biomaRt::getBM(attributes = c(geneid, 'ensembl_gene_id', 'chromosome_name','start_position', 'end_position'), #external_gene_name
                                   filters = 'chromosome_name',
                                   value = chrNs,
                                   mart = ensembl)
        geneinfo[(geneinfo[,1]==""),1]<-geneinfo[(geneinfo[,1]==""),2]
        geneinfo<-geneinfo[,c(1,3,4,5)]
      }else{
        geneinfo <- biomaRt::getBM(attributes = c(geneid, 'chromosome_name','start_position', 'end_position'), #external_gene_name
                                   filters = 'chromosome_name',
                                   value = chrNs,
                                   mart = ensembl)
      }
      # geneinfo <- geneinfo[geneinfo[,1]!="",]
      colnames(geneinfo) <- c("genename", "chr", "start.bp", "end.bp")
    }else if(assembly=="GRCh37"){
      message("load gene information from ensembl (assembly GRCh37)")
      ensembl <- biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh = 37)#, version=version

      if(geneid!="ensembl_gene_id"){
        geneinfo <- biomaRt::getBM(attributes = c(geneid, 'ensembl_gene_id', 'chromosome_name','start_position', 'end_position'), #external_gene_name
                                   filters = 'chromosome_name',
                                   value = chrNs,
                                   mart = ensembl)
        geneinfo[(geneinfo[,1]==""),1]<-geneinfo[(geneinfo[,1]==""),2]
        geneinfo<-geneinfo[,c(1,3,4,5)]
      }else{
        geneinfo <- biomaRt::getBM(attributes = c(geneid, 'chromosome_name','start_position', 'end_position'), #external_gene_name
                                   filters = 'chromosome_name',
                                   value = chrNs,
                                   mart = ensembl)
      }
      # geneinfo <- geneinfo[geneinfo[,1]!="",]
      colnames(geneinfo) <- c("genename", "chr", "start.bp", "end.bp")
    }else{
      stop("wrong assembly")
    }
  }else if(geneDB == "ucsc"){
    if(assembly == "GRCh37"){
      message("load gene information from UCSC genome browser (assembly GRCh37)")
      Homo.sapiens <- Homo.sapiens::Homo.sapiens
      cls <- AnnotationDbi::columns(Homo.sapiens::Homo.sapiens)
      cls <- cls[match(c("TXCHROM","TXEND","TXSTART"), cls)]
      kts <- AnnotationDbi::keytypes(Homo.sapiens)
      if(geneid == 'hgnc_symbol'){
        gene_id <- "SYMBOL"
      }else{
        gene_id <- geneid
      }
      kt <- kts[match(c(gene_id), kts)]
      ks <- AnnotationDbi::keys(Homo.sapiens, keytype=kt)
      res <- AnnotationDbi::select(Homo.sapiens, keys=ks, columns=cls, keytype=kt)
      geneinfo = res[res$TXCHROM %in% paste("chr", chrNs, sep = ""),]
      geneinfo$TXCHROM <- gsub(pattern = "chr", replacement = "", geneinfo$TXCHROM)
      class(geneinfo$TXCHROM)<-"integer"
      colnames(geneinfo) <- c("genename", "chr", "start.bp", "end.bp")

    } else if(assembly == "GRCh38"){
      message("load gene information from UCSC genome browser (assembly GRCh38)")
      Homo.sapiens <- Homo.sapiens::Homo.sapiens
      OrganismDbi::TxDb(Homo.sapiens) <- r.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
      cls <- AnnotationDbi::columns(Homo.sapiens)
      cls <- cls[match(c("TXCHROM","TXEND","TXSTART"), cls)]
      kts <- AnnotationDbi::keytypes(Homo.sapiens)
      if(geneid == 'hgnc_symbol'){
        gene_id <- "SYMBOL"
      }else{
        gene_id <- geneid
      }
      kt <- kts[match(c(gene_id), kts)]
      ks <- AnnotationDbi::keys(Homo.sapiens, keytype=kt)
      res <- AnnotationDbi::select(Homo.sapiens, keys=ks, columns=cls, keytype=kt)
      geneinfo = res[res$TXCHROM %in% paste("chr", chrNs, sep = ""),]
      geneinfo$TXCHROM <- gsub(pattern = "chr", replacement = "", geneinfo$TXCHROM)
      class(geneinfo$TXCHROM)<-"integer"
      colnames(geneinfo) <- c("genename", "chr", "start.bp", "end.bp")
    }else{
      stop("wrong assembly")
    }
  }else{
    stop("wrong DB name!")
  }
  return(geneinfo)
}


#gene overlap merging : combine gene region if there are genes with same name
combineOverlapSamegene = function(genelist){
  combgene <- names(table(genelist[,1])[table(genelist[,1])>1])
  for(cgene in combgene){
    subgeneinfo <- genelist[genelist[,1]==cgene,]
    subgeneinfo <- subgeneinfo[order(subgeneinfo[,3]),]
    genelist <- genelist[-which(genelist[,1]==cgene),]
    addgene <- subgeneinfo[1,,drop=FALSE]
    for(i in 2:dim(subgeneinfo)[1]){
      if(addgene[nrow(addgene),4]<subgeneinfo[i,3]){
        addgene <- rbind(addgene, subgeneinfo[i,])
      }else{
        addgene[nrow(addgene),4] <- max(subgeneinfo[i,4])
      }
    }
    genelist <- rbind(genelist, addgene)
  }
  genelist <- genelist[order(genelist[,3]),,drop=FALSE]
  return(genelist)
}

# MAIN FUNCTIONS ------------------------------------------
# CLQD < input >
#' @title partitioning into cliques
#' @name CLQD
#' @aliases CLQD
#' @description \code{CLQD} partitioning the given data into subgroups that contain SNPs which are highly correlated.
#' @usage CLQD(geno, SNPinfo, CLQcut=0.5, clstgap=40000,
#' hrstType=c("near-nonhrst", "fast", "nonhrst"), hrstParam=200,
#' CLQmode=c("density", "maximal"), LD=c("r2", "Dprime"))
#' @param geno  Data frame or matrix of additive genotype data, each column is additive genotype of each SNP. (Use data of non-monomorphic SNPs)
#' @param SNPinfo  Data frame or matrix of SNPs information. 1st column is rsID and 2nd column is bp position.
#' @param CLQcut Numeric constant; a threshold for the LD measure value |r|, between 0 to 1. Default 0.5.
#' @param clstgap Numeric constant; a threshold of physical distance (bp) between two consecutive SNPs
#' which do not belong to the same clique, i.e., if a physical distance between two consecutive SNPs in a clique
#' greater than \code{clstgap}, then the algorithm split the cliques satisfying each
#' clique do not contain such consecutive SNPs. Default 40000.
#' @param hrstType  Character constant; heuristic methods.
#' If you want to do not use heuristic algorithm, set \code{hrstType = "nonhrst"}.
#' If you want to use heuristic algorithm suggested in Kim et al.,(2017), set \code{hrstType = "fast"}.
#' That algorithm is fastest heuristic algorithm and suitable when your memory capacity is not greater than 8GB.
#' If you want to obtain the results similar to the that of non-heuristic algorithm, set \code{hrstType = "near-nonhrst"}.
#' @param hrstParam Numeric constant; parameter for heuristic algorithm "near-nonhrst".
#'  Default is \code{200}. It is recommended that you set the parameter to greater than 150.
#' @param CLQmode Character constant; the way to give priority among detected cliques.
#' if \code{CLQmode = "density"} then the algorithm gives priority to the clique of largest value of \eqn{(Number of SNPs)/(range of clique)},
#' else if \code{CLQmode = "maximal"}, then the algorithm gives priority to the largest clique. The default is "density".
#' @param LD Character constant; LD measure to use, "r2" or "Dprime". Default "r2".
# <output>
#' @return A vector of cluster numbers of all  SNPs (\code{NA} represents singleton cluster).
#'
#' @seealso \code{\link{BigLD}}
#' @examples
#'
#' data(geno)
#' data(SNPinfo)
#' CLQD(geno=geno[,1:100],SNPinfo=SNPinfo[100,])
#' CLQD(geno=geno[,1:100],SNPinfo=SNPinfo[100,], CLQmode = 'maximal')
#' CLQD(geno=geno[,1:100],SNPinfo=SNPinfo[100,], LD='Dprime')
#'
#' @author Sun-Ah Kim <sunny03@snu.ac.kr>, Yun Joo Yoo <yyoo@snu.ac.kr>
#' @importFrom stats cor median quantile
#' @importFrom utils tail
#' @importFrom Rcpp sourceCpp
#' @export
CLQD <- function(geno, SNPinfo, CLQcut=0.5, clstgap=40000, hrstType=c("near-nonhrst", "fast", "nonhrst"), hrstParam=200,
                 CLQmode=c("density", "maximal"), LD=c("r2", "Dprime"))
{
  ########################################################################################################
  skipRatio <- 0
  CLQmode <- match.arg(CLQmode)
  LD <- match.arg(LD)
  hrstType <- match.arg(hrstType)
  geno <- as.matrix(geno)
  # filtering all NA SNPs
  allNASNPs = apply(geno, 2, function(x) all(is.na(x)))
  geno<-geno[, !allNASNPs]
  SNPinfo<-SNPinfo[!allNASNPs, ]
  # Main Function
  SNPbps <- SNPinfo[, 3]
  if(LD == "r2"){
    OCM <- suppressWarnings(cor(geno, use="pairwise.complete.obs"))
    diag(OCM) <- 0
    OCM[is.na(OCM)]<-0
    OCM[abs(OCM) < CLQcut] <- 0
    OCM[abs(OCM) >= CLQcut] <- 1
  }else if(LD == "Dprime"){
    OCM <- genoDp(geno)
    OCM[(OCM==Inf)|(OCM==-Inf)] <- 0
  }
  binvector <- rep(NA, dim(OCM)[2])
  binnum <- 1
  # re.SNPbps <- SNPbps
  if(all(OCM==0)) return(binvector)

  # test graph complexity
  OCM1 <- OCM
  # bothhigh <- (degs>=400 & cores>=400)
  if(hrstType == "nonhrst"){
    heuristic <- FALSE
  }else if(hrstType == "near-nonhrst"){
    heuristic <- TRUE
    heuristicNum <- 0
    while(heuristic == TRUE){
      g <- igraph::graph_from_adjacency_matrix(OCM1, mode="undirected", weighted=TRUE, diag=FALSE, add.colnames=NA)
      cores <- igraph::coreness(g)
      highcore <- sum(cores>=hrstParam)
      local_cores <- table(cores[cores>=(hrstParam)])
      local_cores <- local_cores[local_cores>=(hrstParam)]
      if(length(local_cores>0)){
        if(heuristicNum==0){
          message("Use near-CLQ heuristic procedure!!")
          heuristicNum <- heuristicNum+1
          message(paste("hueristic loop", heuristicNum))
        }else{
          heuristicNum <- heuristicNum+1
          message(paste("hueristic loop", heuristicNum))
        }
        # find dense region
        # local_cores <- table(cores[cores>=(hrstParam)])
        # local_cores <- local_cores[local_cores>=(hrstParam)]
        local_hrstParam <- names(local_cores[which(as.integer(names(local_cores))==max(as.integer(names(local_cores))))])
        local_hrstParam <- as.numeric(local_hrstParam)
        bothhighSNPs <- which(cores == local_hrstParam)
        SNPset1 <- which(is.na(binvector))[bothhighSNPs]
        nowOCM <- OCM[SNPset1, SNPset1]
        heuristicBins <- heuristicCLQ(nowOCM, hrstParam)
        binvector[SNPset1[heuristicBins]]<-binnum
        binnum <- binnum+1
        OCM1 <- OCM[is.na(binvector), is.na(binvector)]
      }else{
        # if(heuristicNum ==0){
        #   # message("We do not need heuristic procedure")
        # }else{
        #   message("End new heuristic procedure")
        # }
        heuristic <- FALSE
      }
    } #end while
  }else if(hrstType == "fast"){
    checkLargest <- TRUE
    r2Mat <- OCM
    re.SNPbps <- SNPinfo[,3]
    # firstTerm = T
    while(checkLargest == TRUE){
      g <- igraph::graph_from_adjacency_matrix(r2Mat, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NA)
      compo = igraph::components(g)
      componum = which(compo$csize==max(compo$csize))[1]
      compov = which(compo$membership==componum)
      compadjM = OCM1[compov, compov]
      cg = igraph::graph_from_adjacency_matrix(compadjM, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NA)
      if((median(igraph::coreness(cg))>80 & max(igraph::coreness(cg))>100)| (quantile(igraph::coreness(cg), 0.75)>100 & max(igraph::coreness(cg))>100)){
        message("use fast heuristic procedure!")
        # if(quantile(degrees, 0.7) == 1) break
        # message(head((quantile(degrees, 0.7))))
        degrees = apply(r2Mat, 1, sum)
        maxdegv = which(degrees >=max(quantile(degrees, 0.7), 80))
        # if(length(maxdegv)>=1){
        maxdegvs = maxdegv
        edgeDens = NULL
        for(maxdegv in maxdegvs){
          Bignbds = which(r2Mat[maxdegv,, drop = FALSE]>0, arr.ind = TRUE)
          Bignbds.c = unique(Bignbds[,2])
          newr2Mat = r2Mat[Bignbds.c,Bignbds.c]
          EdgeDen = sum(newr2Mat)/((dim(newr2Mat)[1])*(dim(newr2Mat)[1]-1))
          edgeDens = c(edgeDens, EdgeDen)
        }
        maxdegvs = maxdegvs[order(edgeDens, decreasing = TRUE)]
        edgeDens = edgeDens[order(edgeDens, decreasing = TRUE)]
        degv = maxdegvs[1]
        edgeD = edgeDens[1]
        Bignbds = which(r2Mat[degv,, drop = FALSE]>0, arr.ind = TRUE)
        Bignbds.c = unique(Bignbds[,2])
        # maxiC = maximal.cliques(g, min = dim(OCM1)[1]*0.9)
        # largestOneRange = range(Bignbds.c)
        # largestSNPn = diff(largestOneRange)
        # largestCsize = length(Bignbds.c)
        nowSNPsbp = re.SNPbps[Bignbds.c]
        nowSNPsbploca = match(nowSNPsbp, SNPbps)
        binvector[nowSNPsbploca] <- binnum
        binnum = binnum + 1
        r2Mat <- r2Mat[-Bignbds.c, -Bignbds.c, drop = FALSE]
        OCM1 <- OCM1[-Bignbds.c, -Bignbds.c, drop = FALSE]
        re.SNPbps <- re.SNPbps[-Bignbds.c]
        # message("case2")
        checkLargest = TRUE
        if(length(re.SNPbps)<500)  checkLargest = FALSE
      }else{
        checkLargest = FALSE
      }

    }
  }
  # binvector splitting by clstgap
  if(all(is.na(binvector))==FALSE){
    binveclist <- lapply(seq_len(max(binvector, na.rm = TRUE)),
                         function(x) SNPinfo[,3][which(binvector == x)])
    binvecbplist <- newSplitCliques(binveclist, clstgap)
    binvecbpindex <- lapply(binvecbplist, function(x) match(x, SNPinfo[,3]))
    binvecbpindex <- binvecbpindex[vapply(binvecbpindex, length, c(1))>1]
    binvector <- rep(NA, length(binvector))
    for(i in seq_len(length(binvecbpindex))){
      binvector[binvecbpindex[[i]]]<-i
    }
  }
  # take Toooo Big block First!
  r2Mat <- OCM[is.na(binvector), is.na(binvector)]
  if(sum(is.na(binvector))<=1 | (sum(r2Mat)==0)) return(binvector)
  g <- igraph::graph_from_adjacency_matrix(r2Mat, mode="undirected", weighted=TRUE, diag=FALSE, add.colnames=NA)
  max.cliques <- igraph::max_cliques(g, min=2)
  if(length(max.cliques)==0) stop("max.cliques is empty")
  re.SNPbps <- SNPbps[is.na(binvector)]
  bp.cliques <- lapply(max.cliques, function(x) re.SNPbps[x])
  split.bp.cliques <- newSplitCliques(bp.cliques, clstgap)
  # reduce the number of maximal clique? (candidate) or
  # modify density function. narrow SNP distance, so small cliques are chosen preferencely.
  repeat {
    # message(binnum)
    if (all(is.na(binvector) == FALSE)) {
      break
    }
    if(length(split.bp.cliques)==0) break
    now.split.bp.cliques <- split.bp.cliques
    if(CLQmode =="density"){
      density.v <- vapply(now.split.bp.cliques, function(x) ((length(x)))/(diff(range(x))/1000), 1)
    }else{
      density.v <- vapply(now.split.bp.cliques, length, 1)
    }

    max.d <- which(density.v == max(density.v))
    max.cluster <- now.split.bp.cliques[max.d]
    if (length(max.cluster) > 1) {
      # if there are two bins of same density, then we choose the bigger one.
      max.cluster <- max.cluster[order(vapply(max.cluster, length, 1), decreasing=TRUE)]
    }
    max.cluster <- max.cluster[[1]]
    max.cluster.od <- match(max.cluster, re.SNPbps)
    # if (codechange == TRUE) {
    #   max.cluster.od <- ChooseMaximal(max.cluster.od, CLQcut, OCM)
    #   max.cluster <- re.SNPbps[max.cluster.od]
    # }
    ## excluding all SNPs in max.cluster from re.SNPbps
    split.bp.cliques <- lapply(split.bp.cliques, function(x) setdiff(x, max.cluster))
    split.bp.cliques <- unique(split.bp.cliques)
    split.bp.cliques <- split.bp.cliques[which(vapply(split.bp.cliques, length, 1) > 1)]
    binvector[match(max.cluster, SNPbps)] <- binnum
    binnum=binnum + 1
    # r2Mat <- r2Mat[-max.cluster.od, -max.cluster.od]
    # OCM <- OCM[-max.cluster.od, -max.cluster.od]
    # re.SNPbps <- setdiff(re.SNPbps, max.cluster)
    if (length(re.SNPbps) < 2) {
      break
    }
    if(length(split.bp.cliques)==0) break
    # message(sum(is.na(binvector)))
  }  ##end repeat
  return(binvector)
}


# Big-LD <input>
#' @title Estimation of LD block regions
#' @name BigLD
#' @aliases BigLD
#' @description \code{BigLD} returns the estimation of LD block regions of given data.
#' @usage BigLD(geno=NULL, SNPinfo=NULL,genofile=NULL, SNPinfofile=NULL,
#' cutByForce=NULL, LD=c("r2", "Dprime"), CLQcut=0.5,
#' clstgap=40000, CLQmode=c("density", "maximal"),
#' leng=200, subTaskSize=1500, MAFcut=0.05, appendRare=FALSE,
#' hrstType=c("near-nonhrst", "fast", "nonhrst"),
#' hrstParam=200, chrN=NULL, startbp=-Inf, endbp=Inf)
#' @param geno Data frame or matrix of additive genotype data, each column is additive genotype of each SNP.
#' @param SNPinfo Data frame or matrix of SNPs information.  1st column is rsID and 2nd column is bp position.
#' @param genofile Character constant; Genotype data file name (supporting format: .txt, .ped, .raw, .traw, .vcf).
#' @param SNPinfofile Character constant; SNPinfo data file name (supporting format: .txt, .map).
#' @param cutByForce Data frame; information of SNPs which are forced to be the last SNP LD block boundary.
#' @param LD Character constant; LD measure to use, r2 or Dprime. Default "r2".
#' @param CLQcut Numeric constant; a threshold for the LD measure value |r|, between 0 to 1. Default 0.5.
#' @param clstgap Numeric constant; a threshold of physical distance (bp) between two consecutive SNPs
#' which do not belong to the same clique, i.e., if a physical distance between two consecutive SNPs in a clique
#' greater than \code{clstgap}, then the algorithm split the cliques satisfying each
#' clique do not contain such consecutive SNPs. Default 40000.
#' @param CLQmode Character constant; the way to give priority among detected cliques.
#' if \code{CLQmode = "density"} then the algorithm gives priority to the clique of largest value of \eqn{(Number of SNPs)/(range of clique)},
#' else if \code{CLQmode = "maximal"}, then the algorithm gives priority to the largest clique. The default is "density".
#' @param leng Numeric constant; the number of SNPs in a preceding and a following region
#' of each sub-region boundary, every SNP in a preceding and every SNP in a following region need to be in weak LD. Default 200.
#' @param subTaskSize  Numeric constant; upper bound of the number of SNPs in a one-take sub-region. Default 1500.
#' @param MAFcut Numeric constant; the MAF threshold. Default 0.05.
#' @param appendRare logical; if TRUE, the function append rare variants with MAF<MAFcut to the constructed blocks.
#' @param hrstType  Character constant; heuristic methods.
#' If you want to do not use heuristic algorithm, set \code{hrstType = "nonhrst"}.
#' If you want to use heuristic algorithm suggested in Kim et al.,(2017), set \code{hrstType = "fast"}.
#' That algorithm is fastest heuristic algorithm and suitable when your memory capacity is not greater than 8GB.
#' If you want to obtain the results similar to the that of non-heuristic algorithm, set \code{hrstType = "near-nonhrst"}.
#' @param hrstParam Numeric constant; parameter for heuristic algorithm "near-nonhrst".
#'  Default is \code{200}. It is recommended that you set the parameter to greater than 150.
#' @param chrN Numeric(or Character) constant (or vector); chromosome number to use.
#' @param startbp Numeric constant; starting bp position of the \code{chrN}. Default -Inf.
#' @param endbp Numeric constant; last bp position of the \code{chrN}. Default Inf.
# <output>
#' @return  A data frame of block estimation result.
#' Each row of data frame shows the starting SNP and end SNP of each estimated LD block.
#'
#' @author Sun-Ah Kim <sunny03@snu.ac.kr>, Yun Joo Yoo <yyoo@snu.ac.kr>
#'
#'
#'
#' @seealso \code{\link{CLQD}}, \code{\link{LDblockHeatmap}}
#'
#' @examples
#'
#' data(geno)
#' data(SNPinfo)
#' BigLD(geno[,1:100], SNPinfo[1:100,])
#'
#' \dontrun{
#' BigLD(geno, SNPinfo, LD = "Dprime")
#' BigLD(geno, SNPinfo, CLQcut = 0.5, clstgap = 40000, leng = 200, subTaskSize = 1500)
#' }
#' @importFrom stats cor median quantile
#' @importFrom utils tail
#' @importFrom Rcpp sourceCpp
#' @export

BigLD <- function(geno=NULL, SNPinfo=NULL,genofile=NULL, SNPinfofile=NULL, cutByForce=NULL,
                   LD=c("r2", "Dprime"), CLQcut=0.5, clstgap=40000, CLQmode=c("density", "maximal"),
                   leng=200, subTaskSize=1500, MAFcut=0.05, appendRare=FALSE,
                  hrstType=c("near-nonhrst", "fast", "nonhrst"), hrstParam=200,
                   chrN=NULL, startbp=-Inf, endbp=Inf)
{
  skipRatio=0.0
  CLQmode <- match.arg(CLQmode)
  hrstType <- match.arg(hrstType)
  LD <- match.arg(LD)
  # data preparation
  if(!is.null(genofile)){
    inputdata <- dataPreparation(genofile, SNPinfofile, geno, SNPinfo, chrN=chrN, startbp=startbp, endbp=endbp)
    geno <- inputdata[[1]]
    SNPinfo <- inputdata[[2]]
    chrNs <- unique(SNPinfo[,1])
  }else{
    if(is.null(geno) |is.null(SNPinfo)){
      stop("Need input data (or files)")
    }
    if(ncol(geno) != nrow(SNPinfo)){
      stop("The number of SNPs in geno and SNPinfo are not matched")
    }
    geno <- as.matrix(geno)
    # header generation
    if(is.null(chrN)) chrN <- unique(SNPinfo[,1])
    colnames(SNPinfo) <- c("chrN", "rsID", "bp")
    if(startbp != -Inf | endbp != Inf){
      if((length(chrN)>1) | (length(unique(SNPinfo[,1]))>1))
        stop("If your input data include more than one chromosome or you choose more than one chromosome,
             then do not specify the startbp and endbp!")
      SNPloca <- which(SNPinfo$bp>=startbp & SNPinfo$bp<=endbp & SNPinfo$chrN == chrN)
      SNPinfo <- SNPinfo[SNPloca,]
      geno <- geno[,SNPloca]
    }
    colnames(geno)<- SNPinfo$rsID
    chrNs <- unique(SNPinfo[,1])
  }
  # filtering all NA SNPs
  allNASNPs = apply(geno, 2, function(x) all(is.na(x)))
  geno<-geno[, !allNASNPs]
  SNPinfo<-SNPinfo[!allNASNPs, ]

  totalLDBres <- NULL
  totalCLQres <- rep(NA,dim(SNPinfo)[1])
  stindex <- 0

  if(ncol(geno)!=nrow(SNPinfo)) stop("The number of SNPs in geno data and SNPinfo data does not match")
  # for each chrN, seperately apply BigLD algorithm
  for(chrN in chrNs){
    message(paste("working chromosome name:", chrN))
    chrSNPinfo <- SNPinfo[(SNPinfo[,1]==chrN),,drop=FALSE]
    chrgeno <- geno[,(SNPinfo[,1]==chrN),drop=FALSE]
    if(dim(chrSNPinfo)[1]==1)next
    # pruning by MAF ---
    message(paste("remove SNPs with MAF <", MAFcut))
    MAF <- apply(chrgeno, 2, function(x) mean(x,na.rm=TRUE)/2)
    MAF_ok <- ifelse(MAF>=0.5,1-MAF,MAF)
    MAF <- MAF_ok
    message(paste("Number of SNPs after prunning SNPs with MAF<", MAFcut, ":", sum(MAF>=MAFcut)))
    underMAFgeno <- chrgeno[,MAF<MAFcut]
    chrgeno <- chrgeno[,MAF>=MAFcut]
    underMAFSNPinfo <- chrSNPinfo[MAF<MAFcut,]
    chrSNPinfo <- chrSNPinfo[MAF>=MAFcut,]

    chrLDBres <- matrix(NA, dim(chrSNPinfo)[1],2)
    chrCLQres <- rep(NA, dim(chrSNPinfo)[1])
    # divide whole sequence into sub-task size segments.
    chrcutbyforce = cutByForce[cutByForce[,1] == chrN, ]
    cutpoints.all<-cutsequence(geno = chrgeno, SNPinfo = chrSNPinfo, leng = leng,
                               subTaskSize = subTaskSize, LD=LD, clstgap = clstgap, CLQcut = CLQcut, cutByForce = chrcutbyforce)
    cutpoints <- cutpoints.all[[1]]
    atfcut <- (cutpoints.all[[2]])
    if (!is.null(atfcut)){
      atfcut <- sort(atfcut)
    }
    cutpoints <- setdiff(cutpoints, 0)
    cutblock <- cbind(c(1, cutpoints + 1), c(cutpoints, dim(chrgeno)[2]))
    cutblock <- cutblock[-(dim(cutblock)[1]), , drop=FALSE]
    cutblock <- cutblock[cutblock[,1]!=cutblock[,2],, drop=FALSE]
    # for each subtaskregion, adapt BigLD algorithm ------------------
    for(i in seq_len(dim(cutblock)[1])){
      nowst <- cutblock[i, 1]
      nowed <- cutblock[i, 2]
      subgeno <- chrgeno[, nowst:nowed]
      subSNPinfo <- chrSNPinfo[nowst:nowed, ]
      message(paste(Sys.time()," | ", "chr",chrN,":",min(subSNPinfo[,3]), "-" ,max(subSNPinfo[,3]), " | sub-region: ", i,"/",dim(cutblock)[1],  sep = ""))
      subbinvec <- CLQD(geno = subgeno, SNPinfo = subSNPinfo, CLQcut = CLQcut, clstgap = clstgap, hrstType=hrstType,
                        hrstParam = hrstParam, CLQmode = CLQmode, LD = LD)
      # chrCLQres <- c(chrCLQres, subbinvec)
      chrCLQres[nowst: nowed]<-subbinvec
      cat('CLQ done!\r')
      if(all(is.na(subbinvec) == TRUE)) next;
      bins <- seq_len(max(subbinvec[which(!is.na(subbinvec))]))
      clstlist <- lapply(bins, function(x) which(subbinvec  == x))
      clstlist <- lapply(clstlist, sort)  ###
      clstlist <- clstlist[order(vapply(clstlist, min, 1))]  ###
      nowLDblocks <- constructLDblock(clstlist, subSNPinfo)
      nowLDblocks <- nowLDblocks[order(nowLDblocks[, 1]), , drop = FALSE]
      cat('constructLDblock done!\r')
      nowLDblocks.bp <- cbind(subSNPinfo[nowLDblocks[,1],3], subSNPinfo[nowLDblocks[,2],3])
      nowLDblocks <- nowLDblocks + (cutblock[i, 1] - 1)
      nowLDblocks <- nowLDblocks[order(nowLDblocks[, 1]), , drop = FALSE]
      preleng1 <- length(which(!is.na(chrLDBres[, 1])))
      chrLDBres[(preleng1 + 1):(preleng1 + dim(nowLDblocks)[1]), ] <- nowLDblocks

    }

    doneLDblocks <- chrLDBres[which(!is.na(chrLDBres[, 1])), , drop = FALSE]
    largeLDblocks <- doneLDblocks[apply(doneLDblocks, 1, diff)>5,, drop = FALSE]

    # atfcut -- LD block reconst ----------
    # newLDblocks <- matrix(NA, dim(chrSNPinfo)[1], 2)
    if(length(atfcut)!=0){
      addingBlocks <- NULL
      message("handling the ambiguous sub-task regions..")
      for(i in atfcut){
        print(paste("ambiguous boundary", i, "\r"))
        if(length(which(largeLDblocks[,1]>i))==0) break
        st <- largeLDblocks[max(which(largeLDblocks[,2]<=i)),1]
        ed <- largeLDblocks[min(which(largeLDblocks[,1]>i)),2]

        # saveBlocks <- doneLDblocks[doneLDblocks[,1]<st,, drop = FALSE]
        delIndex <- which(doneLDblocks[,1]>=st & doneLDblocks[,2]<=ed)
        doneLDblocks <- doneLDblocks[-delIndex,]
        # if(dim(saveBlocks)[1] != 0){
        #   newLDblocks[which(is.na(newLDblocks)==TRUE)[1]:(which(is.na(newLDblocks)==TRUE)[1] +dim(saveBlocks)[1]),] <- saveBlocks
        # }
        subgeno <- chrgeno[, st:ed]
        subSNPinfo <- chrSNPinfo[st:ed, ]
        subbinvec <- CLQD(geno = subgeno, SNPinfo = subSNPinfo, CLQcut = CLQcut, clstgap = clstgap, LD = LD, hrstType=hrstType,
                          hrstParam = hrstParam)
        # chrCLQres <- c(chrCLQres, subbinvec)
        chrCLQres[st: ed] <- subbinvec  #### change vector procedure
        cat('CLQ done!\r')
        bins <- seq_len(max(subbinvec[which(!is.na(subbinvec))]))
        clstlist <- lapply(bins, function(x) which(subbinvec == x))
        clstlist <- lapply(clstlist, sort)  ###
        clstlist <- clstlist[order(vapply(clstlist, min, 1))]  ###
        nowLDblocks <- constructLDblock(clstlist, subSNPinfo)
        nowLDblocks <- nowLDblocks + (st - 1)
        nowLDblocks <- nowLDblocks[order(nowLDblocks[, 1]), , drop=FALSE]
        addingBlocks <- rbind(addingBlocks, nowLDblocks)
      }
      addingBlocks <- rbind(doneLDblocks, addingBlocks)
      addingBlocks <- addingBlocks[order(addingBlocks[, 2]), , drop=FALSE]
      addingBlocks <- addingBlocks[order(addingBlocks[, 1]), , drop=FALSE]
      for(i in 1:(dim(addingBlocks)[1]-1)){
        if(addingBlocks[i,2]>=addingBlocks[(i+1),2]){
          addingBlocks[(i+1),]<-c(min(addingBlocks[i:(i+1),]), max(addingBlocks[i:(i+1),]))
          addingBlocks[i,1] <- NA
          next
        }
        if(addingBlocks[i,2]>=addingBlocks[(i+1),1]){
          addingBlocks[(i+1),] <-c(min(addingBlocks[i:(i+1),]), max(addingBlocks[i:(i+1),]))
          addingBlocks[i,1] <- NA
          next;
        }
      }
      newLDblocks <- addingBlocks[!is.na(addingBlocks[,1]),]
    }else{
      newLDblocks<-doneLDblocks
    } # LDblock reconst ...atfcut

    start.index <- newLDblocks[,1]
    end.index <- newLDblocks[,2]
    start.rsID <- as.character(chrSNPinfo[newLDblocks[,1],2])
    end.rsID <- as.character(chrSNPinfo[newLDblocks[,2],2])
    start.bp <- chrSNPinfo[newLDblocks[,1],3]
    end.bp <- chrSNPinfo[newLDblocks[,2],3]
    chrLDblocks <- data.frame(start.index, end.index, start.rsID, end.rsID, start.bp, end.bp)
    chrSNPinfo <- SNPinfo[(SNPinfo[,1]==chrN),,drop=FALSE]
    chrLDblocks$start.index<- match(chrLDblocks[,5], chrSNPinfo[,3])
    chrLDblocks$end.index<- match(chrLDblocks[,6], chrSNPinfo[,3])
    #append rare variants --------------------
    if(appendRare==TRUE){
      message( "Append rare variants!")
      chrgeno <- geno[,(SNPinfo[,1]==chrN),drop=FALSE]
      chrLDblocks<-appendcutByForce(chrLDblocks, chrgeno, chrSNPinfo, CLQcut=CLQcut, clstgap=clstgap,
                                    CLQmode=CLQmode, hrstType = hrstType, hrstParam=hrstParam, LD=LD)
    }
    chrLDblocks <- data.frame(chr=chrN, chrLDblocks)
    chrLDblocks$start.index <- chrLDblocks$start.index+stindex
    chrLDblocks$end.index <- chrLDblocks$end.index+stindex
    stindex <- stindex + dim(chrSNPinfo)[1]
    totalLDBres <- rbind(totalLDBres, chrLDblocks)
  } # end for current chr
  message("\nBigLD done!")
  return(totalLDBres)
}

#-GPART
# GPART < input >
#' @title Partitioning genodata based on the result obtained by using Big-LD and gene region information.
#' @name GPART
#' @aliases GPART
#' @description
#' \code{GPART} partition the given genodata using the result obtained by Big-LD and gene region information.
#' The algorithm partition the whole sequence into sub sequences of which size do not exceed the given threshold.
#' @usage GPART(geno=NULL, SNPinfo=NULL, geneinfo=NULL, genofile=NULL,
#' SNPinfofile=NULL, geneinfofile=NULL, geneDB = c("ensembl","ucsc"),
#' assembly = c("GRCh38", "GRCh37"), geneid = "hgnc_symbol",ensbversion = NULL,
#' chrN=NULL, startbp=-Inf, endbp=Inf, BigLDresult=NULL, minsize=4, maxsize=50,
#' LD=c("r2", "Dprime"), CLQcut=0.5, CLQmode=c("density", "maximal"), MAFcut = 0.05,
#' GPARTmode=c("geneBased", "LDblockBased"),
#' Blockbasedmode=c("onlyBlocks", "useGeneRegions"))
#' @param geno Data frame or matrix of additive genotype data, each column is additive genotype of each SNP.
#' @param SNPinfo Data frame or matrix of SNPs information.  1st column is rsID and 2nd column is bp position.
#' @param geneinfo Data frame or matrix of Gene info data.
#' (1st col : Genename, 2nd col : chromosome, 3rd col : start bp,  4th col : end bp)
#' @param genofile Character constant; Genotype data file name (supporting format: .txt, .ped, .raw, .traw, .vcf).
#' @param SNPinfofile Character constant; SNPinfo data file name (supporting format: .txt, .map).
#' @param geneinfofile A Character constant; file containing the gene information
#' (1st col : Genename, 2nd col : chromosome, 3rd col : start bp,  4th col : end bp)
#' @param geneDB A Character constant; database type for gene information. Set \code{"ensembl"} to get gene info from "Ensembl",
#'  or set \code{"ucsc"} to get gene info from "UCSC genome browser" (See package "biomaRt" or
#'  package "homo.sapiens"/"TxDb.Hsapiens.UCSC.hg38.knownGene"/ "TxDb.Hsapiens.UCSC.hg19.knownGene" for details.)
#' @param assembly A character constant; set \code{"GRCh37"} for GRCh37, or set \code{"GRCh38"} for GRCh38
#' @param geneid A character constant; When you use the gene information by \code{geneDB}. specity the symbol for gene name to use.
#' default is "hgnc_symbol".
#' (eg. 'ensembl_gene_id' for \code{geneDB = "ensembl"}, "ENTREZID"/"ENSEMBL"/"REFSEQ"/"TXNAME" for \code{geneDB="ucsc"}.
#' See package 'biomaRt' or package 'Homo.sapiens' for details)
#' @param ensbversion a integer constant; you can set the release version of ensembl
#' when you use the gene information by using \code{geneDB='emsembl'} and \code{assembly='GRCh38'}
#' @param chrN Numeric(or Character) constant (or vector); chromosome number to use. If \code{NULL}(default), we use all chromosome.
#' @param startbp Numeric constant; starting bp position of the \code{chrN}. Default -Inf.
#' @param endbp Numeric constant; last bp position of the \code{chrN}. Default Inf.
#' @param BigLDresult Data frame; a result obtained by \code{BigLD} function. If \code{NULL}(default),
#' the \code{GPART} function first excute \code{BigLD} function to obtain LD blocks estimation result.
#' @param minsize Numeric constant; the lower bound of number of SNPs in a partition.
#' @param maxsize Numeric constant; the upper bound of number of SNPs in a partition.
#' @param LD Character constant; LD measure to use, r2 or Dprime.
#' @param CLQcut Numeric constant; threshold for the correlation value |r|, between 0 to 1.
#' @param CLQmode Character constant; the way to give priority among detected cliques.
#' if \code{CLQmode = "density"} then the algorithm gives priority to the clique of largest value of \eqn{(Number of SNPs)/(range of clique)},
#' else if \code{CLQmode = "maximal"}, then the algorithm gives priority to the largest clique. The default is "density".
#' @param MAFcut Numeric constant; the MAF threshold. Default 0.05.
#' @param GPARTmode Character constant; GPART algorithm methods to use, "geneBased" or "LDblockBased". Default is ??geneBased??
#' @param Blockbasedmode Character constant; When you set \code{GPARTmode = "LDblockBased"}, specify LDblock based method as
#' \code{"onlyBlocks"}("LDblock based only" algorithm) or \code{"useGeneRegions"}(LDblock based and also use gene info algorithm).
#'
# < output >
#' @return \code{GPART} returns data frame which contains 9 information of each partition
#' (chromosome, index number of the first SNP and last SNP, rsID of the first SNP and last SNP,
#' basepair position of the first SNP and last SNP, blocksize, Name of a block)
#'
#' @author Sun Ah Kim <sunny03@snu.ac.kr>, Yun Joo Yoo <yyoo@snu.ac.kr>
#' @seealso \code{\link{BigLD}}
#'
#'
#' @examples
#' data(geno)
#' data(SNPinfo)
#' data(geneinfo)
#' GPART(geno=geno[,1:100], SNPinfo=SNPinfo[1:100,], geneinfo=geneinfo)
#'
#' @importFrom stats cor median quantile
#' @importFrom utils tail
#' @importFrom Rcpp sourceCpp
#' @import Homo.sapiens
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @export
GPART <- function(geno=NULL, SNPinfo=NULL, geneinfo=NULL, genofile=NULL, SNPinfofile=NULL, geneinfofile=NULL,
                  geneDB = c("ensembl","ucsc"), assembly = c("GRCh38", "GRCh37"), geneid = "hgnc_symbol",
                  ensbversion = NULL,
                  chrN=NULL, startbp=-Inf, endbp=Inf, BigLDresult=NULL, minsize=4, maxsize=50,
                  LD=c("r2", "Dprime"), CLQcut=0.5, CLQmode=c("density", "maximal"), MAFcut = 0.05,
                  GPARTmode=c("geneBased", "LDblockBased"), Blockbasedmode=c("onlyBlocks", "useGeneRegions"))
{
  # Main part
  GPARTmode <- match.arg(GPARTmode)
  geneDB <- match.arg(geneDB)
  assembly <- match.arg(assembly)
  CLQmode <- match.arg(CLQmode)
  LD <- match.arg(LD)# if null choose "null"
  # data preparation
  # input data preparation

  if(!is.null(genofile)){
    inputdata <- dataPreparation(genofile, SNPinfofile, geno, SNPinfo, chrN=chrN, startbp=startbp, endbp=endbp)
    geno <- inputdata[[1]]
    SNPinfo <- inputdata[[2]]
  }else{
    geno <- geno
    geno <- as.matrix(geno)
    SNPinfo <- SNPinfo
    if(is.null(geno) |is.null(SNPinfo)){
      stop("Need input data (or files)")
    }
    if(ncol(geno) != nrow(SNPinfo)){
      stop("The number of SNPs in geno and SNPinfo are not matched")
    }
    if(is.null(chrN)) chrN <- unique(SNPinfo[,1])
    colnames(SNPinfo) <- c("chrN", "rsID", "bp")
    if(startbp != -Inf | endbp != Inf){
      if((length(chrN)>1) | (length(unique(SNPinfo[,1]))>1))
        stop("If your input data include more than one chromosome or you choose more than one chromosome,
             then do not specify the startbp and endbp!")
      SNPloca <- which(SNPinfo$bp>=startbp & SNPinfo$bp<=endbp & SNPinfo$chrN == chrN)
      SNPinfo <- SNPinfo[SNPloca,]
      geno <- geno[,SNPloca]
    }
    colnames(geno)<- SNPinfo$rsID
  }
  chrNs <- unique(SNPinfo[,1])

  # filtering all NA SNPs
  allNASNPs = apply(geno, 2, function(x) all(is.na(x)))
  geno<-geno[, !allNASNPs]
  SNPinfo<-SNPinfo[!allNASNPs, ]
  #SNP filtering
  MAF <- apply(geno, 2, function(x) mean(x,na.rm=TRUE)/2)
  MAF_ok <- ifelse(MAF>=0.5,1-MAF,MAF)
  MAF <- MAF_ok
  message(paste("Number of SNPs after prunning SNPs with MAF<", MAFcut, ":", sum(MAF>=MAFcut)))
  geno <- geno[,MAF>=MAFcut]
  SNPinfo <- SNPinfo[MAF>=MAFcut,]

  if(is.null(BigLDresult)){
    message("Start to execute BigLD function!")
    BigLDresult <- BigLD(geno, SNPinfo, CLQcut=CLQcut, CLQmode=CLQmode, LD=LD, MAFcut=MAFcut)#, MAFcut=MAFcut
    message("Big-LD, done!")
  }else{
    message("Use the inputted BigLD result")
    # BigLDblocks <- BigLDresult #[BigLDresult$chr==chrN,]
    newindex.st <- match(BigLDresult$start.rsID, SNPinfo[,2])
    newindex.ed <- match(BigLDresult$end.rsID, SNPinfo[,2])
    #BigLDresult check
    if(sum(is.na(newindex.st))!=0 |sum(is.na(newindex.ed))!=0){
      stop("The inputted data and the bigLD results are not paired.
           \nPlease check whether the BigLD results are obtained from the input data")
    }
    BigLDresult$start.index<-newindex.st
    BigLDresult$end.index<-newindex.ed
    }
  # gene based
  # load gene info
  if(is.null(geneinfo)){
    geneinfo <- geneinfoExt(geneinfofile=geneinfofile, geneDB = geneDB, assembly = assembly, geneid = geneid,
                             ensbversion = ensbversion, chrNs = chrNs)
  }
  if(GPARTmode == "geneBased"){
    message(paste("GPART algorithm: geneBased"))
    TotalRes <- NULL
    startindex <- 0
    for(chrN in unique(SNPinfo[,1])){
      message(paste("working chromosome name:", chrN))
      chrSNPinfo <- SNPinfo[(SNPinfo[,1]==chrN),,drop=FALSE]
      chrgeno <- geno[,(SNPinfo[,1]==chrN),drop=FALSE]
      BigLDblocks <- BigLDresult[BigLDresult$chr==chrN,,drop=FALSE]

      genelist <- geneinfo[which(geneinfo[,2] == chrN), ]
      colnames(genelist)<- c("genename", "chr", "start.bp", "end.bp")
      genelist <- genelist[order(genelist[,3]),,drop=FALSE]
      genelist <- combineOverlapSamegene(genelist)
      GeneRegionSNPs1 <- NULL
      for (i in seq_len(dim(genelist)[1])){
        test <- (which(chrSNPinfo[, 3] >= genelist[i, 3] & chrSNPinfo[, 3] <= genelist[i, 4]))
        ifelse(length(test) > 0, test <- range(test), test <- c(0, 0))
        GeneRegionSNPs1 <- rbind(GeneRegionSNPs1, test)
        # if(i%%10) message(i)
      }
      rownames(GeneRegionSNPs1) <- genelist[, 1]
      SNPexist <- apply(GeneRegionSNPs1, 1, function(x) (x[1] != 0 & x[2] != 0))
      SNPexist <- unlist(SNPexist)
      Geneblocks <- GeneRegionSNPs1[SNPexist, ,drop=FALSE]
      if(dim(Geneblocks)[1]>0){
        Geneblocks.M <- mergeOverlapGene(Geneblocks)
      }else{
        Geneblocks.M <- Geneblocks
      }
      if(dim(BigLDblocks)[1]>0){
        newindex.st <- match(BigLDblocks$start.rsID, chrSNPinfo[,2])
        newindex.ed <- match(BigLDblocks$end.rsID, chrSNPinfo[,2])
        #BigLDresult check
        # if(sum(is.na(newindex.st))!=0 |sum(is.na(newindex.ed))!=0){
        #   stop("The inputted data and the bigLD results are not paired.
        #      \nPlease check whether the BigLD results are obtained from the input data")
        # }
        BigLDblocks$start.index<-newindex.st
        BigLDblocks$end.index<-newindex.ed
        nowLDblocks <- BigLDblocks
        LDblocks <- cbind(nowLDblocks[, 2], nowLDblocks[, 3])
        # Split Large LDblocks
        FinalLDblocks <- LDblockSplit(chrgeno, LDblocks, maxsize, LD)
        LDblocks.sgt <- setdiff(seq_len(dim(chrSNPinfo)[1]), unlist(apply(FinalLDblocks, 1, function(x) min(x):max(x))))
        LDblocks.sgt <- cbind(LDblocks.sgt, LDblocks.sgt)
        LDblocks.T <- rbind(FinalLDblocks, LDblocks.sgt)
        LDblocks.T <- LDblocks.T[order(LDblocks.T[, 1]), ]
        colnames(LDblocks.T) <- c("st", "ed")
        # gene-base block partitioning
      }else{
        LDblocks.T <- cbind(seq_len(dim(chrSNPinfo)[1]), seq_len(dim(chrSNPinfo)[1]))
        colnames(LDblocks.T) <- c("st", "ed")
      }
      if(dim(Geneblocks.M)[1]>0){
        GeneLDblocks2 <- LDblockGeneMerge(LDblocks.T, Geneblocks.M)
      }else{
        GeneLDblocks2 <- LDblocks.T
        rownames(GeneLDblocks2) <- rep("No", dim(GeneLDblocks2)[1])
      }
      blockL <- GeneLDblocks2[, 2] - GeneLDblocks2[, 1] + 1
      GeneLDblocks3 <- cbind(GeneLDblocks2, blockL)
      # split big block
      if(sum(blockL>maxsize)>0){
        GeneLDblocks4 <- splitBigLD(GeneLDblocks3, LDblocks.T, Geneblocks, maxsize)
      }else{
        GeneLDblocks4 <- GeneLDblocks3
        # GeneLDblocks4 <- data.frame(GeneLDblocks4, rep("No", dim(GeneLDblocks4)[1]))
        # colnames(GeneLDblocks4)[4] <-"gname"
      }
      GeneLDblocks5 <- mergeSmallRegion(GeneLDblocks4, maxsize, minsize)
      # naming regions
      GeneLDblocks <- namingRegion2(GeneLDblocks5, genelist, chrSNPinfo)
      st.rsID <- chrSNPinfo[GeneLDblocks[, 1], 2]
      ed.rsID <- chrSNPinfo[GeneLDblocks[, 2], 2]
      Finalresult <- data.frame(chrN,GeneLDblocks$st, GeneLDblocks$ed, st.rsID, ed.rsID, GeneLDblocks$st.bp, GeneLDblocks$ed.bp,
                                GeneLDblocks$blockL, GeneLDblocks$Finalnames)
      colnames(Finalresult) <- c("chr","st", "ed", "st.rsID", "ed.rsID", "st.bp", "ed.bp", "blocksize", "Name")
      Finalresult$st <- Finalresult$st+startindex
      Finalresult$ed <- Finalresult$ed+startindex
      startindex <- startindex + dim(chrgeno)[2]
      message(paste("chr", chrN, "done!"))
      TotalRes <- rbind(TotalRes, Finalresult)
    }
    # LD block based
  }else if (GPARTmode == "LDblockBased"){
    Blockbasedmode <- match.arg(Blockbasedmode)
    message(paste("GPART algorithm: LDblockBased-", Blockbasedmode, sep = ""))
    TotalRes <- NULL
    startindex <- 0
    for(chrN in unique(SNPinfo[,1])){
      message(paste("working chromosome name:", chrN))
      chrSNPinfo <- SNPinfo[(SNPinfo[,1]==chrN),,drop=FALSE]
      chrgeno <- geno[,(SNPinfo[,1]==chrN),drop=FALSE]
      BigLDblocks <- BigLDresult[BigLDresult$chr==chrN,,drop=FALSE]

      genelist <- geneinfo[which(geneinfo[,2] == chrN), ]
      colnames(genelist)<- c("genename", "chr", "start.bp", "end.bp")
      genelist <- genelist[order(genelist[,3]),]
      genelist <- combineOverlapSamegene(genelist)
      # gene-region merging
      GeneRegionSNPs1 <- NULL
      for (i in seq_len(dim(genelist)[1])){
        test <- (which(chrSNPinfo[, 3] >= genelist[i, 3] & chrSNPinfo[, 3] <= genelist[i, 4]))
        ifelse(length(test) > 0, test <- range(test), test <- c(0, 0))
        GeneRegionSNPs1 <- rbind(GeneRegionSNPs1, test)
        # if(i%%10) message(i)
      }
      rownames(GeneRegionSNPs1) <- genelist[, 1]
      SNPexist <- apply(GeneRegionSNPs1, 1, function(x) (x[1] != 0 & x[2] != 0))
      SNPexist <- unlist(SNPexist)
      Geneblocks <- GeneRegionSNPs1[SNPexist, ,drop=FALSE]
      if(dim(Geneblocks)[1]>0){
        Geneblocks.M <- mergeOverlapGene(Geneblocks)
      }else{
        Geneblocks.M <- Geneblocks
      }

      if(dim(BigLDblocks)[1]>0){
        newindex.st <- match(BigLDblocks$start.rsID, chrSNPinfo[,2])
        newindex.ed <- match(BigLDblocks$end.rsID, chrSNPinfo[,2])
        BigLDblocks$start.index <- newindex.st
        BigLDblocks$end.index <- newindex.ed
        nowLDblocks <- BigLDblocks
        LDblocks <- cbind(nowLDblocks[, 2], nowLDblocks[, 3])
        # Split Large LDblocks
        FinalLDblocks <- LDblockSplit(chrgeno, LDblocks, maxsize, LD)
        LDblocks.sgt <- setdiff(seq_len(dim(chrSNPinfo)[1]), unlist(apply(FinalLDblocks, 1, function(x) min(x):max(x))))
        LDblocks.sgt <- cbind(LDblocks.sgt, LDblocks.sgt)
        LDblocks.T <- rbind(FinalLDblocks, LDblocks.sgt)
        LDblocks.T <- LDblocks.T[order(LDblocks.T[, 1]), ]
        colnames(LDblocks.T) <- c("st", "ed")
      }else{
        LDblocks.T <- cbind(seq_len(dim(chrSNPinfo)[1]), seq_len(dim(chrSNPinfo)[1]))
        colnames(LDblocks.T) <- c("st", "ed")
      }
      # onlyblocks? or useGeneRegions?
      if(Blockbasedmode == "onlyBlocks"){
        blockL <- apply(LDblocks.T, 1, function(x) diff(x)+1)
        LDblocks.T <- data.frame(LDblocks.T, blockL)
        LDblocks1 <- mergeSmallRegion(LDblocks.T, maxsize, minsize)
        GeneLDblocks <- namingRegion2(LDblocks1, genelist, chrSNPinfo)
        # end: if(Blockbasedmode == "onlyBlocks")
      }else if(Blockbasedmode == "useGeneRegions"){
        blockL <- apply(LDblocks.T, 1, function(x) diff(x)+1)
        LDblocks.T <- data.frame(LDblocks.T, blockL)
        gLDblocks <- data.frame(LDblocks.T, NA)
        for(k in seq_len(dim(gLDblocks)[1])){
          nowblock <- gLDblocks[k,,drop=FALSE]
          intro.index <- (which(Geneblocks[,2]<nowblock[1,1]))
          intro.index1 <- setdiff(seq_len(dim(Geneblocks)[1]), intro.index)
          outro.index <- (which(Geneblocks[,1]>nowblock[1,2]))
          outro.index1 <- setdiff(seq_len(dim(Geneblocks)[1]), outro.index)
          common.index <- intersect(intro.index1, outro.index1)
          if(length(common.index)==0) {
            gLDblocks[k,4]<-NA
          }else{
            gnames <- rownames(Geneblocks[common.index,1, drop=FALSE])
            gname <- paste(gnames, collapse="/")
            gLDblocks[k,4]<-gname
          }
        }
        ## LD block
        bigblocks <- gLDblocks[gLDblocks$blockL>=minsize,]
        smallblocks <- gLDblocks[gLDblocks$blockL<minsize,]
        small.gene <- smallblocks[!is.na(smallblocks[,4]),]
        small.nongene <- smallblocks[is.na(smallblocks[,4]),]
        #merge small blocks which are in gene regions
        addedblocks <- NULL
        generegions <- unique(small.gene[,4])
        for(g in generegions){
          # message(g)
          nowgeneblocks <- small.gene[(small.gene[,4]==g),,drop=FALSE]
          nowbin <- NULL
          while(dim(nowgeneblocks)[1]>0){
            if(is.null(nowbin)){
              nowbin <- nowgeneblocks[1,,drop=FALSE]
              nowgeneblocks <- nowgeneblocks[-1,,drop=FALSE]
            }
            if(dim(nowgeneblocks)[1]==0 & !is.null(nowbin)){
              addedblocks<- rbind(addedblocks, nowbin)
              break
            }
            newbin <- nowgeneblocks[1,,drop=FALSE]
            nowgeneblocks <- nowgeneblocks[-1,,drop=FALSE]
            if(nowbin[1,3]+newbin[1,3]>maxsize){
              addedblocks<- rbind(addedblocks, nowbin)
              nowbin <- newbin
            }else if((nowbin[1,2]+1)!=newbin[1,1]){
              addedblocks<- rbind(addedblocks, nowbin)
              nowbin <- newbin
            }else if((nowbin[1,2]+1)==newbin[1,1]){
              nowbin[1,2] <- newbin[1,2]
              nowbin[1,3] <- nowbin[1,3]+newbin[1,3]
            }
            if(dim(nowgeneblocks)[1]==0 & !is.null(nowbin)){
              addedblocks<- rbind(addedblocks, nowbin)
              break
            }
          }
        }
        GeneLDblocks <- rbind(bigblocks, addedblocks, small.nongene)
        GeneLDblocks <- GeneLDblocks[order(GeneLDblocks[,1]),]
        LDblocks1 <- mergeSmallRegion(GeneLDblocks[,seq_len(3)], maxsize, minsize)
        GeneLDblocks <- namingRegion2(LDblocks1, genelist, chrSNPinfo)
      } # end if(Blockbasedmode == "useGeneRegions")
      st.rsID <- chrSNPinfo[GeneLDblocks[, 1], 2]
      ed.rsID <- chrSNPinfo[GeneLDblocks[, 2], 2]
      Finalresult <- data.frame(chrN,GeneLDblocks$st, GeneLDblocks$ed, st.rsID, ed.rsID, GeneLDblocks$st.bp, GeneLDblocks$ed.bp,
                                GeneLDblocks$blockL, GeneLDblocks$Finalnames)
      colnames(Finalresult) <- c("chr","st", "ed", "st.rsID", "ed.rsID", "st.bp", "ed.bp", "blocksize", "Name")
      Finalresult$st <- Finalresult$st+startindex
      Finalresult$ed <- Finalresult$ed+startindex
      startindex <- startindex + dim(chrgeno)[2]
      message(paste("chr", chrN, "done!"))
      TotalRes <- rbind(TotalRes, Finalresult)
    }
  }
  colnames(TotalRes)[seq_len(7)] <- c( "chr","start.index", "end.index", "start.rsID","end.rsID", "start.bp", "end.bp")
  return(TotalRes)
}


#' @title visualization of LD block structure
#' @name LDblockHeatmap
#' @aliases LDblockHeatmap
#' @description
#' \code{LDblockHeatmap} visualize the LD structure or LD block results of the inputed data.
#' LDblockHeatmap <- function(geno=NULL, SNPinfo=NULL, genofile=NULL, SNPinfofile=NULL,geneinfo=NULL, geneinfofile = NULL,
#' geneDB = c("ensembl","ucsc","file"), assembly = c("GRCh38", "GRCh37"), geneid = "hgnc_symbol",
#' ensbversion = NULL, chrN=NULL,startbp=-Inf, endbp=Inf, blockresult=NULL, blocktype=c("bigld", "gpart"),
#' minsize=4, maxsize=50, LD=c("r2", "Dprime", "Dp-str"),  MAFcut=0.05,CLQcut=0.5,
#' CLQmode=c("density", "maximal"),  CLQshow=FALSE,
#' type=c("png", "tif"), filename="heatmap", res=300, onlyHeatmap=FALSE)
#' @param geno Data frame or matrix of additive genotype data, each column is additive genotype of each SNP.
#' @param SNPinfo Data frame or matrix of SNPs information.  1st column is rsID and 2nd column is bp position.
#' @param genofile Character constant; Genotype data file name (supporting format: .txt, .ped, .raw, .traw, .vcf).
#' @param SNPinfofile Character constant; SNPinfo data file name (supporting format: .txt, .map).
#' @param geneinfo Data frame or matrix of Gene info data.
#' (1st col : Genename, 2nd col : chromosome, 3rd col : start bp,  4th col : end bp)
#' @param geneinfofile A Character constant; file containing the gene information
#' (1st col : Genename, 2nd col : chromosome, 3rd col : start bp,  4th col : end bp)
#' @param geneDB A Character constant; database type for gene information. Set \code{"ensembl"} to get gene info from "Ensembl",
#'  or set \code{"ucsc"} to get gene info from "UCSC genome browser" (See package "biomaRt" or
#'  package "homo.sapiens"/"TxDb.Hsapiens.UCSC.hg38.knownGene"/ "TxDb.Hsapiens.UCSC.hg19.knownGene" for details.)
#' @param assembly A character constant; set \code{"GRCh37"} for GRCh37, or set \code{"GRCh38"} for GRCh38
#' @param geneid A character constant; When you use the gene information by \code{geneDB}. specity the symbol for gene name to use.
#' default is "hgnc_symbol".
#' (eg. 'ensembl_gene_id' for \code{geneDB = "ensembl"}, "ENTREZID"/"ENSEMBL"/"REFSEQ"/"TXNAME" for \code{geneDB="ucsc"}.
#' See package 'biomaRt' or package 'Homo.sapiens' for details)
#' @param ensbversion a integer constant; you can set the release version of ensembl
#' when you use the gene information by using \code{geneDB='emsembl'} and \code{assembly='GRCh38'}
#' @param geneshow logical; do not show the gene information if \code{geneshow=FALSE}.
#' Default is \code{geneshow=TRUE}
#' @param chrN Numeric(or Character) constant ; chromosome number to use.
#' If the data contains more than one chromosome, you need to specify the chromosome to show.
#' @param startbp Numeric constant; starting bp position of the \code{chrN}.
#' @param endbp Numeric constant; last bp position of the \code{chrN}.
#' @param blockresult Data frame; a result obtained by \code{BigLD} function or \code{GPART}. If \code{NULL}(default),
#' the function first excute \code{BigLD} or\code{GPART} to obtain block estimation result denpending on the \code{blocktype}.
#' @param blocktype Character constant; "bigld" for \code{Big-LD} or "gpart" for \code{GPART}. Default is \code{"gpart"}.
#' @param minsize Integer constant; when \code{blockresult=NULL, blocktype="gpart"}
#' specify the threshold for minsize of a block obtained by \code{GPART}
#' @param maxsize Integer constant; when \code{blockresult=NULL, blocktype="gpart"}
#' specify the threshold for maxsize of a block obtained by \code{GPART}
#' @param LD Character constant; LD measure to use, "r2" or "Dprime" or "Dp-str".
#' LD measures for LD heatmap (and BigLD execution when \code{BigLDresult=NULL}). When \code{LD = "Dp-str"},
#' heatmap shows only two cases, "weak LD or not-informative" and "strong LD". When \code{LD= Dprime} heatmap shows the estimated D' measures.
#' @param MAFcut Numeric constant; MAF threshold of SNPs to use.  Default 0.05
#' @param CLQcut Numeric constant; threshold for the correlation value |r|, between 0 to 1. Default 0.5.
#' @param CLQmode Character constant; the way to give priority among detected cliques.
#' if \code{CLQmode = "density"} then the algorithm gives priority to the clique of largest value of \eqn{(Number of SNPs)/(range of clique)},
#' else if \code{CLQmode = "maximal"}, then the algorithm gives priority to the largest clique. Default "density".
#' @param CLQshow logical; Show the LD bin structures, i.e. CLQ results, in each LD blocks.
#' Notice that if the region to show includes more than 200 SNPs, the function do now show LD bin structures automatically. Default false.
#' @param type Character constant; file format of output image file. set \code{png} for PNG format file, or \code{tif} for TIFF format file.
#' @param filename Character constant; filename of output image file. Default "heatmap".
#' @param res Numeric constant; resolution of image. Default 300.
#' @param onlyHeatmap logical; show the LD heatmap without the LD block boundaries. Default false.

# < output >
#' @return \code{GPART} returns data frame which contains 9 information of each partition
#' (chromosome, index number of the first SNP and last SNP, rsID of the first SNP and last SNP,
#' basepair position of the first SNP and last SNP, blocksize, Name of a block)
#'
#' @author Sun-Ah Kim <sunny03@snu.ac.kr>, Yun Joo Yoo <yyoo@snu.ac.kr>
#' @seealso \code{\link{BigLD}}
#'
#' @examples
#'
#' LDblockHeatmap(geno=geno[,1:100], SNPinfo=SNPinfo[1:100,], geneinfo=geneinfo,
#' filename="chr21Heatmap")
#'
#' @importFrom grDevices adjustcolor colorRampPalette dev.off heat.colors png tiff
#' @importFrom stats cor median quantile
#' @importFrom utils tail
#' @import grid
#' @importFrom Rcpp sourceCpp
#' @export
#-LDblockHeatmap2
LDblockHeatmap <- function(geno=NULL, SNPinfo=NULL, genofile=NULL, SNPinfofile=NULL,geneinfo=NULL, geneinfofile=NULL,
                           geneDB=c("ensembl","ucsc","file"), assembly=c("GRCh38", "GRCh37"), geneid="hgnc_symbol",
                           ensbversion=NULL, geneshow=TRUE, chrN=NULL,startbp=-Inf, endbp=Inf, blockresult=NULL, blocktype=c("bigld", "gpart"),
                           minsize=4, maxsize=50, LD=c("r2", "Dprime", "Dp-str"),  MAFcut=0.05,CLQcut=0.5,
                           CLQmode=c("density", "maximal"),  CLQshow=FALSE,
                           type=c("png", "tif"), filename="heatmap", res=300, onlyHeatmap=FALSE)
{
  # -set data and parameters ------
  tick <- "both"
  maxisize <- 20000
  longleng <- 2000
  shortleng <- ifelse(is.null(geneinfo), 1000, 100)
  CLQmode <- match.arg(CLQmode)
  blocktype <- match.arg(blocktype)
  type <- match.arg(type)
  LD <- match.arg(LD)
  LDcal <- ifelse(LD=="r2", "r2", "Dprime")
  # tick <- match.arg(tick)
  ## data preparation ----------

  # input data preparation
  if(!is.null(genofile)){
    inputdata <- dataPreparation(genofile, SNPinfofile, geno, SNPinfo, chrN=chrN, startbp=startbp, endbp=endbp)
    geno <- inputdata[[1]]
    SNPinfo <- inputdata[[2]]
  }else{
    if(is.null(geno) |is.null(SNPinfo)){
      stop("Need input data (or files)")
    }
    if(ncol(geno) != nrow(SNPinfo)){
      stop("The number of SNPs in geno and SNPinfo are not matched")
    }
    if(is.null(chrN)){
      if(length(unique(SNPinfo[,1]))==1){
        chrN <- unique(SNPinfo[,1])
      }else{
        stop("Please choose a chromosome!")
      }
    }
    colnames(SNPinfo) <- c("chrN", "rsID", "bp")
    geno <- as.matrix(geno)
  }
  if(is.null(geneinfo)){
    if(geneshow == TRUE){
      geneinfo <- geneinfoExt(geneinfofile=geneinfofile, geneDB = geneDB, assembly = assembly, geneid = geneid,
                               ensbversion = ensbversion, chrNs = chrN)
    }
  }

  # filtering all NA SNPs
  allNASNPs = apply(geno, 2, function(x) all(is.na(x)))
  geno<-geno[, !allNASNPs]
  SNPinfo<-SNPinfo[!allNASNPs, ]
  ## MAF filtering
  MAF <- apply(geno, 2, function(x) mean(x,na.rm=TRUE)/2)
  MAF_ok <- ifelse(MAF>=0.5,1-MAF,MAF)
  MAF <- MAF_ok
  message(paste("Number of SNPs after prunning SNPs with MAF<", MAFcut, ":", sum(MAF>=MAFcut)))
  # underMAFgeno <- geno[,MAF<MAFcut]
  geno <- geno[,MAF>=MAFcut]
  # underMAFSNPinfo <- SNPinfo[MAF<MAFcut,]
  SNPinfo <- SNPinfo[MAF>=MAFcut,]

  trascolor <-  adjustcolor( "red", alpha.f=0)
  col1 <- c('#7D0112','#820920','#87112B','#8C1834','#911E3D','#952445','#9A294C','#9E2F54','#A2345B','#A73962',
            '#AA3F69','#AE4470','#B24977','#B64E7D','#B95483','#BC598A','#C05E90','#C36396','#C6689C','#C86EA1',
            '#CB73A7','#CE78AC','#D07DB2','#D282B7','#D587BC','#D78CC1','#D991C5','#DB96CA','#DC9BCE','#DEA0D3',
            '#E0A5D7','#E1AADB','#E3AFDF','#E4B3E2','#E5B8E6','#E7BCE9','#E8C1EC','#E9C5EF','#EACAF1','#EBCEF4',
            '#ECD2F6','#EDD6F8','#EEDAF9','#EFDDFB','#F0E1FC','#F1E4FC','#F1E8FC','#F2EBFC','#F2EDFA','#F2F0F6')
  col2 <- c('#2D3184','#283C87','#23458A','#1E4E8E','#185692','#125E95','#0A6699','#026E9D','#0075A0','#007CA3',
            '#0083A6','#008AA9','#0790AC','#1296AE','#1C9CB0','#26A2B2','#2FA8B4','#38ADB6','#40B3B7','#49B8B9',
            '#51BCBA','#5AC1BB','#62C5BC','#6AC9BD','#73CDBE','#7BD1BE','#83D5BF','#8AD8C0','#92DBC1','#99DEC1',
            '#A0E0C2','#A7E3C3','#AEE5C4','#B5E7C5','#BBE9C6','#C1EBC7','#C7ECC9','#CCEECA','#D2EFCB','#D7F0CD',
            '#DBF0CF','#DFF1D0','#E3F2D2','#E7F2D4','#EAF2D6','#EDF2D9','#EFF2DB','#F1F2DD','#F2F2E0','#F3F1E4')
  col3 <- c('#D33F6A','#D44468','#D54866','#D74C63','#D85161','#D9545F','#DB585C','#DC5C5A','#DD6057','#DE6455',
            '#DF6752','#E06B50','#E16F4D','#E2724A','#E37647','#E47A45','#E57D42','#E5813F','#E6843C','#E78839',
            '#E78B37','#E88E34','#E89231','#E8952F','#E9992D','#E99C2B','#E99F2A','#E9A328','#EAA628','#EAAA28',
            '#EAAD29','#EAB02A','#E9B42C','#E9B72E','#E9BA31','#E9BE35','#E8C139','#E8C43D','#E8C742','#E7CB46',
            '#E7CE4C','#E6D151','#E6D457','#E5D75D','#E4DA64','#E4DD6B','#E3E073','#E2E37C','#E2E688','#E2E6BD')
  col4 <- c('#035F33','#176438','#23693C','#2D6D41','#367245','#3E764A','#457B4F','#4D7F53','#538458','#5A885C',
            '#608C61','#679165','#6D9569','#73996E','#799D72','#7EA176','#84A57A','#8AA97F','#8FAD83','#94B087',
            '#99B48B','#9EB88F','#A3BB93','#A8BE97','#ADC29B','#B1C59F','#B6C8A2','#BACBA6','#BECEAA','#C3D1AD',
            '#C7D4B1','#CAD6B4','#CED9B7','#D2DBBB','#D5DEBE','#D8E0C1','#DBE2C4','#DEE4C7','#E1E6CA','#E4E8CD',
            '#E6E9CF','#E9EBD2','#EBECD5','#EDEDD7','#EEEFD9','#F0EFDC','#F1F0DE','#F2F1E0','#F3F1E2','#F3F1E4')
  col5 <- c("#8B4497","#8D499B","#8F4E9F","#9052A4","#9257A8","#945CAC","#9560B0","#9765B4","#9869B8","#9A6DBC",
            "#9B72C0","#9D76C4","#9E7AC8","#9F7FCC","#A083D0","#A187D4","#A28BD8","#A390DB","#A494DF","#A598E2",
            "#A69CE6","#A7A0E9","#A7A4ED","#A8A8F0","#A8ACF3","#A9B0F7","#A9B4FA","#AAB8FD","#AABCFF","#AABFFF",
            "#ABC3FF","#ABC7FF","#ABCBFF","#ABCEFF","#ABD2FF","#ABD5FF","#ABD9FF","#ABDCFF","#ABDFFF","#AAE3FF",
            "#AAE6FF","#AAE9FF","#A9ECFF","#A9EFFF","#A8F1FF","#A7F4FF","#A7F6FF","#A6F9FF","#A4FBFF","#A3FCFF")


  allcol <- cbind(col5, col2, col3, col4)
  blockcol <- c(col2[5], col3[5], col4[5],col5[5])
  deepblockcol <- c(col2[1], col3[1], col4[1],col5[1])

  if(is.null(chrN)){
    chrN <- unique(SNPinfo[,1])
    if(length(chrN)>1) stop("The data contains two or more chromosomes. Please specify the chromosome.")
  }
  SNPloca <- which(SNPinfo$chrN==chrN & SNPinfo$bp>=startbp & SNPinfo$bp<=endbp)
  NSNPs <- length(SNPloca)
  if(NSNPs<5){
    stop("There are less than 5 SNPs in the chosen region!")
  }
  if(NSNPs>maxisize){
    Pmessage <- paste("There are too many SNPs in the chosen region. \n We use the first ",maxisize,"SNPs in the region!")
    message(Pmessage)
    SNPloca <- min(SNPloca):(min(SNPloca)+maxisize-1)
    NSNPs <- length(SNPloca)
  }
  if(NSNPs>200){
    if(CLQshow == TRUE){
      message("CLQ results can not be displayed because the number of SNPs in the region is too large.")
    }
    CLQshow <- FALSE
  }
  geno <- geno[,SNPloca]
  geno <- as.matrix(geno)
  SNPinfo <- SNPinfo[SNPloca,]

  SNPstbp <- min(SNPinfo[,3])
  SNPedbp <- max(SNPinfo[,3])

  SNPbp_cordi <- (SNPinfo[,3]-SNPstbp)/(SNPedbp - SNPstbp)

  ## matching BigLD index and SNP/genodata index
  if((!is.null(blockresult)) & (onlyHeatmap == FALSE)){
    if(blocktype == "gpart"){
      blockresult <- blockresult[,seq_len(7)]
      colnames(blockresult)[seq_len(7)] <- c( "chr","start.index", "end.index", "start.rsID","end.rsID", "start.bp", "end.bp")
    }
    blockresult <- blockresult[which(blockresult$chr==chrN & blockresult$end.bp>=SNPstbp & blockresult$start.bp<=SNPedbp),,drop=FALSE]
    blockresult$start.bp[blockresult$start.bp < SNPstbp] <- min(SNPinfo[,3])
    blockresult$end.bp[blockresult$end.bp > SNPedbp] <- max(SNPinfo[,3])
    stindex <- match(blockresult$start.bp, SNPinfo[,3])
    edindex <- match(blockresult$end.bp, SNPinfo[,3])
    blockresult$start.index <- stindex
    blockresult$end.index <- edindex
    blockresult$start.rsID <- SNPinfo[stindex,2]
    blockresult$end.rsID <- SNPinfo[edindex,2]
    blockresult <- blockresult[blockresult[,2]!=blockresult[,3],,drop=FALSE]
  }else if ((is.null(blockresult)) & (onlyHeatmap == FALSE)){
    if(blocktype == "bigld"){
      oblockresult <- BigLD(geno=geno, SNPinfo=SNPinfo, CLQcut=CLQcut, CLQmode=CLQmode, LD=LDcal, MAFcut=MAFcut)
      blockresult <- oblockresult
    }else if(blocktype == "gpart"){
      oblockresult <- GPART(geno=geno, SNPinfo=SNPinfo, geneinfo = geneinfo, CLQcut=CLQcut, CLQmode=CLQmode, LD=LDcal,
                            minsize = minsize, maxsize = maxsize, MAFcut = MAFcut)
      blockresult <- oblockresult[,seq_len(7)]
      colnames(blockresult)[seq_len(7)] <- c( "chr","start.index", "end.index", "start.rsID","end.rsID", "start.bp", "end.bp")
    }
  }

  if(NSNPs < 1000){
    unitleng <- 1.1
    unitlwd <- 0.4
  }else if(NSNPs >= 1000 & NSNPs < 4000){
    unitleng <- 2
    unitlwd <- 0.3
  }else if(NSNPs >= 4000 & NSNPs < 7000){
    unitleng <- 3
    unitlwd <- 0.2
  }else if(NSNPs >= 7000 & NSNPs <= maxisize){
    unitleng <- 4
    unitlwd <- 0.1
  }

  if(!is.null(geneinfo)){
    nowgeneinfo <- geneinfo[geneinfo[,2]==chrN,,drop=FALSE]
    nowgeneinfo <- nowgeneinfo[which(nowgeneinfo[,3]<SNPedbp & nowgeneinfo[,4]>SNPstbp),,drop=FALSE]
    nowgeneinfo <- nowgeneinfo[order(nowgeneinfo[,3]),,drop=FALSE]
    nowgeneinfo <- combineOverlapSamegene(nowgeneinfo)
  }else{
    message("There is no gene region information!")
    nowgeneinfo <- matrix(1,1,1)[-1,,drop=FALSE]
  }
  stratifyGene <- list()
  if(dim(nowgeneinfo)[1]>0){
    nowgeneinfo[which(nowgeneinfo[,3]<SNPstbp),3] <- SNPstbp
    nowgeneinfo[which(nowgeneinfo[,4]>SNPedbp),4] <- SNPedbp
    nowgeneinfo[,3:4] <- (nowgeneinfo[,3:4]-SNPstbp)/(SNPedbp - SNPstbp)
    stratifyGene <- list()
    # stratifyGene[[1]]<-nowgeneinfo[1,,drop=F]
    # nowgeneinfo<-nowgeneinfo[-1,]
    for(i in seq_len(dim(nowgeneinfo)[1])){
      nowgene <- nowgeneinfo[i,,drop=FALSE]
      appendlistnum <- which(vapply(stratifyGene, function(x) (max(x[,4])+0.01) < nowgene[1,3], TRUE)==TRUE)
      if(length(appendlistnum)==0){
        stratifyGene[[length(stratifyGene)+1]]<-nowgene
      }else{
        stratifyGene[[appendlistnum[1]]]<-rbind(stratifyGene[[appendlistnum[1]]], nowgene)
      }
    }
  }
  if(CLQshow == TRUE){
    geneheight <- min(0.05*length(stratifyGene), 0.35)
    bottomloca <- 0.33 + geneheight
  }else{
    geneheight <- min(0.05*length(stratifyGene), 0.35)
    bottomloca <- 0.25 + geneheight
  }
  if(NSNPs<300){
    nowcex <- 0.5
  }else if(NSNPs > 300 & NSNPs<1000) {
    nowcex <- 0.4
  }else if(NSNPs >= 1000 & NSNPs<=2000) {
    nowcex <- 0.4
  }else if(NSNPs >= 2000 & NSNPs<=3000) {
    nowcex <- 0.4
  }else if(NSNPs>=3000) {
    nowcex <- 0.3
  }
  # start drawing part ---------------------------------------------------------------------------------
  # heatmap construct---------------------------------------------------------------------------------
  ##color
  if(LD == 'Dprime'){
    RWcolor <- colorRampPalette(c("red", "white"))
    color <- RWcolor(50)
  }else if (LD == 'r2'){
    color=heat.colors(50)
  }else if (LD == "Dp-str"){
    RWcolor <- colorRampPalette(c("red", "white"))
    color <- RWcolor(50)[c(1,50)]
  }
  mybreak<-0:length(color)/length(color)
  VPheatmap<-viewport(x=0,y=0,width=unit(0.85/1.414,"snpc"),height=unit(0.85/1.414,"snpc"),just=c("left","bottom"), name="VPheatmap")
  rot.vp <- viewport(x=.11, y=bottomloca, width=unit(unitleng,"snpc"), height=unit(unitleng,"snpc"), just=c("left", "bottom"), angle=-45)
  # genodata
  if(NSNPs > longleng){
    # M<-Matrix::Matrix(0, dim(SNPinfo)[1], dim(SNPinfo)[1])
    startseq <- seq(1, dim(SNPinfo)[1],shortleng)
    # system.time({
    rectlist <- NULL
    for(i in startseq){
      if(LD == "r2"){
        imgLDmatrix <- suppressWarnings(cor(geno[,(i:min(i+shortleng-1, NSNPs)),drop=FALSE], geno[,(i:min(i+longleng+shortleng, NSNPs)),drop=FALSE],
                                            use="pairwise.complete.obs")^2)
      }else if(LD == "Dprime"){
        imgLDmatrix <- matCubeX2(geno[,(i:min(i+shortleng-1, NSNPs)),drop=FALSE], geno[,(i:min(i+longleng+shortleng, NSNPs)),drop=FALSE])
      }else if(LD == "Dp-str"){
        imgLDmatrix <- genoDp2(geno[,(i:min(i+shortleng-1, NSNPs)),drop=FALSE], geno[,(i:min(i+longleng+shortleng, NSNPs)),drop=FALSE],
                               strLD=TRUE)
      }
      imgLDmatrix[(imgLDmatrix==Inf)|(imgLDmatrix==-Inf)] <- 0
      imgLDmatrix[lower.tri(imgLDmatrix)] <- NA
      colcut <- as.character(cut(1-imgLDmatrix,mybreak,labels=as.character(color), include.lowest=TRUE, include.highest=FALSE))
      xx <- ((i:min(i+shortleng-1, NSNPs)))/NSNPs
      yy <- ((i):min(i+longleng+shortleng, NSNPs))/NSNPs
      right <- rep(xx, length((i):min(i+longleng+shortleng, NSNPs)))
      top <- rep(yy, each=length((i:min(i+shortleng-1, NSNPs))))
      nowrect <- rectGrob(x=right,y=top,width=1/NSNPs,height=1/NSNPs,just=c("right","top"),gp=gpar(col=NA,fill=colcut))
      cat(paste("Heatmap generating",i, "|", NSNPs, "\r"))
      rectlist <- gList(rectlist, nowrect)
    }
    LD.heat.map<-gTree(children=rectlist,name="LDheatmap")
    # }) #system.time
  }else{
    if(LD == "r2"){
      imgLDmatrix <- (cor(geno,use="pairwise.complete.obs"))^2
    }else if(LD == "Dprime"){
      imgLDmatrix <- matCubeX(geno)
    }else if(LD == "Dp-str"){
      imgLDmatrix <- genoDp(geno,strLD=TRUE)
    }
    imgLDmatrix[(imgLDmatrix==Inf)|(imgLDmatrix==-Inf)] <- 0
    imgLDmatrix[lower.tri(imgLDmatrix)] <- NA
    colcut <- as.character(cut(1-imgLDmatrix,mybreak,labels=as.character(color), include.lowest=TRUE, include.highest=FALSE))
    xx <- seq_len(NSNPs)/NSNPs
    yy <- seq_len(NSNPs)/NSNPs
    right <- rep(xx, NSNPs)
    top <- rep(yy, each=NSNPs)
    nowrect <- rectGrob(x=right,y=top,width=1/NSNPs,height=1/NSNPs,just=c("right","top"),gp=gpar(col=NA,fill=colcut))
    rectlist <- gList(nowrect)
    LD.heat.map <- gTree(children=rectlist,name="LDheatmap")
  }

  if(onlyHeatmap == FALSE){
    ### block boundaries
    BlockstP <- blockresult[,2]-1
    BlockedP <- blockresult[,3]
    x1s <- BlockstP*(1/NSNPs)
    x2s <- BlockstP*(1/NSNPs)
    x3s <- BlockedP*(1/NSNPs)
    y1s <- BlockstP*(1/NSNPs)
    y2s <- BlockedP*(1/NSNPs)
    y3s <- BlockedP*(1/NSNPs)
    BlockBoundaries <- polygonGrob(x=c(x1s, x2s, x3s), y=c(y1s, y2s, y3s), id=rep(seq_len(length(x1s)),3),
                                   gp=gpar(fill=trascolor, lwd=unitlwd, col="black"))

    LDblockHeatmap2 <- gTree(children=gList(LD.heat.map, BlockBoundaries),name="LDHeatmap", vp=VPheatmap)
    blocksizes <- blockresult[,3]-blockresult[,2]+1
    blocksizes <- sort(blocksizes, decreasing=TRUE)
    if(length(blocksizes) > (15*unitleng)){
      showLDsize <- blocksizes[(15*unitleng)]
      if(length(blocksizes[blocksizes>=showLDsize])>(15*unitleng)) showLDsize <- showLDsize+1
    }else{
      showLDsize <- min(blocksizes)
    }
    regionstbp <- (blockresult[1,6])
    regionedbp <- (blockresult[dim(blockresult)[1],7])
    rotLDtest <- gTree(children=gList(LDblockHeatmap2), name="rorLDheatmap", vp=rot.vp)
    # CLQ res show
    if(CLQshow == TRUE){
      CLQres <- NULL
      for(bn in seq_len(dim(blockresult)[1])){
        nowblock <- blockresult[bn,,drop=FALSE]
        nowgeno <- geno[,nowblock[1,2]:nowblock[1,3]]
        nowSNPinfo <- SNPinfo[nowblock[1,2]:nowblock[1,3],]
        nowCLQres <- CLQD(geno=nowgeno, SNPinfo=nowSNPinfo, CLQcut=0.5, clstgap=40000, hrstParam=200, hrstType="nonhrst",
                          CLQmode=CLQmode, LD=LDcal)
        CLQres <- c(CLQres, list(nowCLQres))
      }

      CLQmaxnum <- max(unlist(CLQres), na.rm=TRUE)
      CLQwidth <- 1/NSNPs
      CLQboxgL <- gList()

      CLQregionheight <- (0.015-0.002*unitleng)*CLQmaxnum
      vpCLQ <- viewport(x=0.11,y=bottomloca-CLQregionheight,width=unit(unitleng*0.85,"snpc"),height=unit(CLQregionheight,"snpc"),just=c("left","bottom"), name="vpCLQ")

      totalbox <- rectGrob(x=0, y=0, width=1, height=1, just=c("left","bottom"),
                           gp=gpar(fill="white",lwd=unitlwd, col="black"), vp=vpCLQ)
      CLQboxgL <- gList(CLQboxgL, totalbox)
      ifelse(CLQmaxnum<25, colgap <- (50%/%CLQmaxnum), colgap <- 1)
      bgcol <- c("gray60", "gray90")
      for(strt in seq_len(length(CLQres))){
        nowbgcol <- bgcol[strt%%2+1]
        nowCLQres <- CLQres[[strt]]
        stpoint <- blockresult[strt, 2]-1
        edpoint <- blockresult[strt, 3]
        leftbottomP <- (stpoint:edpoint)/(NSNPs)
        totalbox <- rectGrob(x=leftbottomP[1], y=0, width=CLQwidth*length(nowCLQres), height=1,
                             just=c("left","bottom"), gp=gpar(fill=nowbgcol, lwd=unitlwd, col="black"), vp=vpCLQ)
        CLQboxgL <- gList(CLQboxgL, totalbox)
        nowcols <- allcol[,(strt%%4+1)]
        clqnum <- max(nowCLQres, na.rm=TRUE)
        for(CLQn in seq_len(clqnum)){
          nowx <- leftbottomP[which((nowCLQres==CLQn)==TRUE)]
          nowx <- nowx[!is.na(nowx)]

          nownowcol <- nowcols[1+(colgap*(CLQn-1))]
          nowbox <- rectGrob(x=nowx, y=(1/CLQmaxnum)*(CLQn-1), width=CLQwidth, height=(1/CLQmaxnum),
                             just=c("left","bottom"), gp=gpar(fill=nownowcol ,lwd=unitlwd, col="black"), vp=vpCLQ)
          CLQboxgL <- gList(CLQboxgL, nowbox)
        }
      }
      #write.snpnumber
      binvector <- rep(0, NSNPs)
      for(strt in seq_len(length(CLQres))){
        stpoint <- blockresult[strt, 2]
        edpoint <- blockresult[strt, 3]
        binvector[stpoint:edpoint]<-CLQres[[strt]]
      }
      binvector[is.na(binvector)]<-0
      textx <- seq_len(NSNPs)/(NSNPs)-(0.5/NSNPs)
      texty <- binvector*(1/max(binvector)) - (0.5/max(binvector))
      textcol <- rep(NA, NSNPs)
      textcol[which(binvector <= round(max(binvector)/2))] <- "white"
      textcol[which(binvector > round(max(binvector)/2))] <- "black"
      clqtextsize =0.4
      if(NSNPs>50) {clqtextsize <- 0.3}
      if(NSNPs>100) {clqtextsize <- 0.2}
      CLQtext <- textGrob(seq_len(NSNPs), x=textx, y=texty, just=c("centre", "centre"), rot=90,
                          gp=gpar(col=textcol, cex=clqtextsize), vp=vpCLQ)
      CLQboxgL <- gList(CLQboxgL, CLQtext)
    }else{
      CLQboxgL <- gList()
      CLQregionheight <- 0
    }

    # oriblockresult <- blockresult
    vpSNPinfo <- viewport(x=0.11,y=0.90,width=unit(unitleng*0.85,"snpc"),height=unit(0.10,"snpc"),just=c("left","bottom"), name="vpSNPinfo")
    if(tick == "rsID"){
      tickname.st <- as.character(blockresult$start.rsID)
      tickname.ed <- as.character(blockresult$end.rsID)
      ticknames <- paste(tickname.st,"\n" ,tickname.ed, sep="")
    } else if(tick =="bp"){
      tickname.st <- (blockresult$start.bp)
      tickname.ed <- (blockresult$end.bp)
      ticknames <- paste(tickname.st,"\n" ,tickname.ed, sep="")
    }else if(tick == "both"){
      tickname.st1 <- as.character(blockresult$start.rsID)
      tickname.ed1 <- as.character(blockresult$end.rsID)
      tickname.st2 <- round(blockresult$start.bp*0.001)
      tickname.ed2 <- round(blockresult$end.bp*0.001)
      tickname.st <- paste(tickname.st1, "\n(",tickname.st2,"kb)", sep="")
      tickname.ed <- paste(tickname.ed1, "\n(",tickname.ed2,"kb)", sep="")
      ticknames <- paste(tickname.st,"\n" ,tickname.ed, sep="")
    }

    #SNP name write
    # blockresult <- blockresult[(blockresult[,3]-blockresult[,2]+1)>=showLDsize,]
    oriblockresult <- blockresult
    ticknames <- ticknames[(blockresult[,3]-blockresult[,2]+1)>=showLDsize]
    rectGrobposi <- seq_len(length(ticknames))/length(ticknames)
    Blockcol <- rep(blockcol, dim(blockresult)[1])[seq_len(nrow(blockresult))]
    Blockcol <- Blockcol[(blockresult[,3]-blockresult[,2]+1)>=showLDsize]
    LDSNPticksBox <- rectGrob(x=rectGrobposi, y=0.01, height=0.05, width=unit( 0.45, "snpc"),
                              gp=gpar(fill =Blockcol, col=Blockcol, alpha=0.7), just=c("centre","bottom"), vp=vpSNPinfo)
    # SNPnameposi=sort(c(rectGrobposi-min(2/length(rectGrobposi),0.01), rectGrobposi+min(2/length(rectGrobposi),0.01)))
    LDSNPticks <- textGrob(ticknames, x=rectGrobposi, y=0.1, just=c("left", "centre"),
                           rot=90, gp=gpar(cex=0.4, lineheight=0.8), vp= vpSNPinfo)
    #SNPinfo box mappingline draw
    vpSNPinfomapping <- viewport(x=0.11,y=bottomloca,width=unit(unitleng*0.85,"snpc"),height=unit(0.9-bottomloca,"snpc"),
                                 just=c("left","bottom"), name="vpSNPinfo")
    sblockresult <- blockresult[(blockresult[,3]-blockresult[,2]+1)>=showLDsize,]
    Blockcenter.x <- apply(cbind((sblockresult[,2]-1)*(1/NSNPs), sblockresult[,3]*(1/NSNPs)), 1, mean)
    Blockcenter.y <- apply(cbind((sblockresult[,2]-1)*(1/NSNPs), sblockresult[,3]*(1/NSNPs)), 1, function(x) (diff(x)/2))
    blocktop <- max(Blockcenter.y*unitleng*0.85)/(0.9-bottomloca)
    SNPblockmapping <- segmentsGrob(x0=rectGrobposi, x1=Blockcenter.x,
                                    y0=1, y1=blocktop, vp=vpSNPinfomapping, gp=gpar(lwd=unitlwd, col=Blockcol))
    SNPblockmapping1 <- segmentsGrob(x0=Blockcenter.x, x1=Blockcenter.x,
                                     y0=blocktop, y1=(Blockcenter.y*unitleng*0.85)/(0.9-bottomloca),
                                     vp=vpSNPinfomapping, gp=gpar(lwd=unitlwd, col=Blockcol))
    #SNP box name wright
    Blocknames <- paste("B-", which(oriblockresult[,3]-oriblockresult[,2]+1>=showLDsize), sep="")
    LDblockticks <- textGrob(Blocknames, x=rectGrobposi, y=-0.01, just=c("centre", "top"),
                             gp=gpar(cex=0.4), vp=vpSNPinfo)

    # (Blockcenter.y*unitleng*0.85)/(0.35)
    #
    # grid.draw(Locatick)
    #rotate LD heatmap
    #mapping lines
    #2
    mappingheight<- ifelse(CLQshow==FALSE, 0.13, 0.13)
    vpMappingLine <- viewport(x=0.11,y=bottomloca-(CLQregionheight+mappingheight),width=unit(unitleng*0.85,"snpc"),
                              height=unit(mappingheight,"snpc"),just=c("left","bottom"), name="vpmappingLine")
    # x0 <- c((blockresult[,2]-1)*(1/NSNPs), blockresult[,3]*(1/NSNPs))
    x0 <- apply(cbind((oriblockresult[,2]-1)*(1/NSNPs), oriblockresult[,3]*(1/NSNPs)), 1, mean)
    y0 <- rep(1, length(x0))
    # x1ori <- c(blockresult[,6], blockresult[,7])
    x1ori <- apply(cbind(oriblockresult[,6], oriblockresult[,7]), 1, mean)
    x1 <- (x1ori-SNPstbp)/(SNPedbp-SNPstbp)
    y1 <- rep(0, length(x1))
    mappingLine <- segmentsGrob(x0=x0, y0=y0, x1=x1, y1=y1,vp=vpMappingLine,gp=gpar(lwd=unitlwd, col=blockcol))
    # blockname on loca bar wright
    Blocknames <- paste("B-", which(oriblockresult[,3]-oriblockresult[,2]+1>=showLDsize), sep="")
    x1ori <- apply(cbind(sblockresult[,6], sblockresult[,7]), 1, mean)
    x1 <- (x1ori-SNPstbp)/(SNPedbp-SNPstbp)
    LDblocktickslocaBar <- textGrob(Blocknames, x=x1, y=0, just=c("left", "centre"), rot=90,
                                    gp=gpar(cex=0.4), vp=vpMappingLine)
    #location Bar
    locabarheight <- 0.02
    vpLocaBar <- viewport(x=0.11,y=(bottomloca-(CLQregionheight+mappingheight+locabarheight)),
                          width=unit(unitleng*0.85,"snpc"),height=unit(locabarheight,"snpc"),just=c("left","bottom"), name="vpLocaBar")
    locaBar <- rectGrob(x=0, y=0, width=1, height=1, just=c("left","bottom"), vp=vpLocaBar, gp=gpar(col="gray50", fill="gray50"))
    #draw LD blocks on the location bar
    x1 <- (blockresult[,6]-SNPstbp)/(SNPedbp-SNPstbp)
    blockwidth <- (blockresult[,7]-blockresult[,6])/(SNPedbp-SNPstbp)
    blockonBar <- rectGrob(x=x1, y=rep(0, dim(blockresult)[1]), height=1, width=blockwidth,
                           just=c("left","bottom"), gp=gpar(col=deepblockcol, fill=blockcol, lwd=0.1), vp=vpLocaBar)

    x0 <- seq_len(NSNPs)/NSNPs - (1/(2*NSNPs))
    y0 <- rep(1, length(x0))
    # x1ori <- c(blockresult[,6], blockresult[,7])
    x1 <- SNPbp_cordi
    y1 <- rep(0, length(x1))
    blocksizes <- blockresult[,3]-blockresult[,2]+1
    sgtblocks <- setdiff(seq_len(NSNPs), unlist(apply(blockresult[,2:3], 1, function(x) x[1]:x[2])))
    sgtblocks <- cbind(sgtblocks, sgtblocks)
    totalBigLDres <- rbind(as.matrix(blockresult[,2:3]), sgtblocks)
    totalBigLDres <- totalBigLDres[order(totalBigLDres[,1]),, drop=FALSE]
    SNPonBlockcol <- c()
    j<-0
    for(i in seq_len(dim(totalBigLDres)[1])){
      if(totalBigLDres[i, 2] == totalBigLDres[i, 1] ){
        SNPonBlockcol <-c(SNPonBlockcol, "gray50")
      }else{
        j<-j+1
        jj = ifelse(j%%4==0, 4, j%%4)
        SNPonBlockcol <-c(SNPonBlockcol, rep(blockcol[jj], (totalBigLDres[i, 2] - totalBigLDres[i, 1] + 1)))
      }
    }
    if(NSNPs>1000){
      x0 <- x0[seq(1,length(x0), 10)]
      x1 <- x1[seq(1,length(x1), 10)]
      y0 <- y0[seq(1,length(y0), 10)]
      y1 <- y1[seq(1,length(y1), 10)]
      SNPonBlockcol=SNPonBlockcol[seq(1,length(SNPonBlockcol), 10)]
    }

    deepSNPonBlockcol <- c()
    j<-0
    for(i in seq_len(dim(totalBigLDres)[1])){
      if(totalBigLDres[i, 2] == totalBigLDres[i, 1] ){
        deepSNPonBlockcol <-c(deepSNPonBlockcol, "gray50")
      }else{
        j<-j+1
        jj = ifelse(j%%4==0, 4, j%%4)
        deepSNPonBlockcol <-c(deepSNPonBlockcol, rep(deepblockcol[jj], (totalBigLDres[i, 2] - totalBigLDres[i, 1] + 1)))
      }
    }
    if(NSNPs>1000){
      deepSNPonBlockcol=deepSNPonBlockcol[seq(1,length(deepSNPonBlockcol), 10)]
    }

    SNPmappingLine <- segmentsGrob(x0=x0, y0=y0, x1=x1, y1=y1,vp=vpMappingLine,gp=gpar(lwd=unitlwd, col=SNPonBlockcol))
    # blockname on loca bar wright
    locaBar <- rectGrob(x=0, y=0, width=1, height=1, just=c("left","bottom"), vp=vpLocaBar, gp=gpar(col="gray50", fill="gray50"))
    SNPsonBar <- segmentsGrob(x0 =x1, x1=x1, y0=0.05, y1=0.95, vp=vpLocaBar, gp=gpar(lwd=unitlwd, col=deepSNPonBlockcol))
  }else{
    # no block boundaries, only heatmap
    LDblockHeatmap2 <- gTree(children=gList(LD.heat.map),name="LDHeatmap", vp=VPheatmap)
    rotLDtest <- gTree(children=gList(LDblockHeatmap2), name="rorLDheatmap", vp=rot.vp)
    CLQshow = FALSE
    CLQregionheight <- 0
    mappingheight <- 0.13
    locabarheight <- 0.02

    mappingheight<- ifelse(CLQshow==FALSE, 0.13, 0.05)
    vpMappingLine <- viewport(x=0.11,y=bottomloca-(CLQregionheight+mappingheight),width=unit(unitleng*0.85,"snpc"),
                              height=unit(mappingheight,"snpc"),just=c("left","bottom"), name="vpmappingLine")
    # x0 <- c((blockresult[,2]-1)*(1/NSNPs), blockresult[,3]*(1/NSNPs))
    x0 <- seq_len(NSNPs)/NSNPs - (1/(2*NSNPs))
    y0 <- rep(1, length(x0))
    # x1ori <- c(blockresult[,6], blockresult[,7])
    x1 <- SNPbp_cordi
    y1 <- rep(0, length(x1))
    if(NSNPs>1000){
      x0 <- x0[seq(1,length(x0), 10)]
      x1 <- x1[seq(1,length(x1), 10)]
      y0 <- y0[seq(1,length(y0), 10)]
      y1 <- y1[seq(1,length(y1), 10)]
      # SNPonBlockcol=SNPonBlockcol[seq(1,length(y1), 10)]
    }
    SNPmappingLine <- segmentsGrob(x0=x0, y0=y0, x1=x1, y1=y1,vp=vpMappingLine,gp=gpar(lwd=unitlwd, col="gray50"))
    # blockname on loca bar wright
    vpLocaBar <- viewport(x=0.11,y=(bottomloca-(CLQregionheight+mappingheight+locabarheight)),width=unit(unitleng*0.85,"snpc"),
                          height=unit(locabarheight,"snpc"),just=c("left","bottom"), name="vpLocaBar")
    locaBar <- rectGrob(x=0, y=0, width=1, height=1, just=c("left","bottom"), vp=vpLocaBar, gp=gpar(col="gray50", fill="gray50"))
    SNPsonBar <- segmentsGrob(x0 =x1, x1=x1, y0=0.05, y1=0.95, vp=vpLocaBar, gp=gpar(lwd=unitlwd, col="black"))
    CLQboxgL <- gList()
    LDSNPticksBox <- gList()
    LDSNPticks <- gList()
    SNPblockmapping <- gList()
    SNPblockmapping1 <- gList()
    LDblockticks <- gList()
    LDblocktickslocaBar <- gList()
    blockonBar <- gList()
    mappingLine <- gList()
  }

  # gene loca info ---------

  if(dim(nowgeneinfo)[1]>0){
    generegionH <- bottomloca-(CLQregionheight+mappingheight+locabarheight)-0.045
    stratNum <- length(stratifyGene)
    stratheight <- min(generegionH/stratNum, 0.05)
    geneinfoList <- gList()
    for(i in seq_len(stratNum)){
      # each info region : height 0.5
      nowvp <- viewport(x=0.11,y=(bottomloca-(CLQregionheight+mappingheight+locabarheight)-(stratheight*i)),
                        width=unit(unitleng*0.85,"snpc"),height=unit((stratheight*0.9),"snpc"),just=c("left","bottom"))
      nowgenes <- stratifyGene[[i]]
      x0s <- nowgenes[,3]
      x1s <- nowgenes[,4]
      nowwidth <- x1s-x0s
      # nowwidth[which(nowwidth<0.01)]<-0.008
      nowcol <- c("blue", "green3", "darkmagenta")
      geneinfoList <- gList(geneinfoList,
                            rectGrob(x=x0s, y=rep(0.8, length(x0s)), width=nowwidth, height=0.2,
                                     vp=nowvp, gp=gpar(fill=nowcol, col=nowcol, cex=nowcex), just=c("left","bottom")),
                            textGrob(as.character(nowgenes[,1]),x=x0s, y=c(0.75, 0.45), vp=nowvp, gp=gpar(col=nowcol, cex=0.4),
                                     just=c("left", "top")),
                            segmentsGrob(x0=x0s, x1=x0s, y0=0.8, y1=c(0.75, 0.45),vp=nowvp,
                                         gp=gpar(col=nowcol, cex=0.4, lty="dotted", alpha=0.5))
      )
    }
    geneMap <- gTree(children=geneinfoList,name="Genemap")
    # grid.draw(geneMap)
  }else{
    if(!is.null(geneinfo)){message("This region does not overlap any gene region!")}
    stratheight <- 0;stratNum <- 0
    geneMap<-NULL
    text_gene <- gList()
  }
  vpLocationTick <- viewport(x=0.11,y=bottomloca-(CLQregionheight+mappingheight),width=unit(unitleng*0.85,"snpc"),height=unit((stratheight*stratNum)+0.03,"snpc"),
                             just=c("left","top"), name="vpLocaTick")
  vpLocationTickname <- viewport(x=0.11,y=bottomloca-(CLQregionheight+mappingheight+locabarheight)-(stratheight*stratNum)-0.05,width=unit(unitleng*0.85,"snpc"),height=unit(0.05,"snpc"),
                                 just=c("left","bottom"), name="vpLocaTick")
  x0range <- SNPedbp-SNPstbp
  x0num <- length(strsplit(as.character(x0range), split="")[[1]])-1
  x0tickround <- seq(round(SNPstbp/(10^x0num)), round(SNPedbp/(10^x0num)),0.5)
  x0tickloca <- x0tickround*10^x0num
  x0tickloca <- x0tickloca[which(x0tickloca>=SNPstbp & x0tickloca<=SNPedbp)]
  x0tickloca <- c(SNPstbp, x0tickloca, SNPedbp)
  x0tickloca.adj <- (x0tickloca-SNPstbp)/x0range
  Locatick <- segmentsGrob(x0=x0tickloca.adj, y0=rep(1, length(x0tickloca.adj)),
                           x1=x0tickloca.adj, y1=rep(0, length(x0tickloca.adj)),
                           vp= vpLocationTick, gp=gpar(lty="dashed", col="gray50"))
  # locationTickname

  if(x0num >= 6){
    Locatickname <- textGrob(paste(round(x0tickloca/10^6,1), "Mb", sep=""),
                             x=x0tickloca.adj, y=0.9, rot=-45, just=c("left", "top"),vp=vpLocationTickname,gp=gpar(cex=0.5) )
  }else if((x0num >= 3) & (x0num < 6)){
    Locatickname <- textGrob(paste(round(x0tickloca/10^3,1), "Kb", sep=""),
                             x=x0tickloca.adj, y=0.9, rot=-45, just=c("left", "top"),vp=vpLocationTickname,gp=gpar(cex=0.5))
  }else if(x0num < 3){
    Locatickname <- textGrob(paste(x0tickloca, "bp", sep=""),
                             x=x0tickloca.adj, y=0.9, rot=-45, just=c("left", "top"),vp=vpLocationTickname,gp=gpar(cex=0.5))
  }

  #Locatick2
  x0tickround1 <- seq(round(SNPstbp/(10^(x0num-1))), round(SNPedbp/(10^(x0num-1))),1)
  x0tickloca1 <- x0tickround1*(10^(x0num-1))
  x0tickloca1 <- x0tickloca1[which(x0tickloca1 >= SNPstbp & x0tickloca1 <= SNPedbp)]
  x0tickloca1 <- c(SNPstbp, x0tickloca1, SNPedbp)
  x0tickloca1<-setdiff(x0tickloca1, x0tickloca)
  x0tickloca.adj1 <- (x0tickloca1-SNPstbp)/x0range
  Locatick1 <- segmentsGrob(x0=x0tickloca.adj1, y0=rep(1, length(x0tickloca.adj1)),
                            x1=x0tickloca.adj1, y1=rep(0, length(x0tickloca.adj1)),
                            vp= vpLocationTick, gp=gpar(lty="dotted", col="gray50"))

  #color key ----------------
  LDheatmapLegend<- function(color, LD) {
    colorrect <- rectGrob(x=(length(color):1/length(color)),y=0.7,width=1/length(color),height=0.3,just=c("right","top"),gp=gpar(col=NA,fill=color))
    if(LD == "r2"){
      ttt <-expression(paste("r"^{2}, " Color Key"))
    }else {
      ttt <-"D' Color Key"
    }
    title <- textGrob(ttt, x=0.5, y=0.9, name="title",
                      gp=gpar(cex=0.5))
    if(LD != "Dp-str"){
      labels <- textGrob(paste(0.2 * 0:5), x=0.2 * 0:5, y=0.2,
                         gp=gpar(cex=0.3), name="labels")
      ticks <- segmentsGrob(x0=c(0:5) * 0.2, y0=rep(0.4,6), x1=c(0:5) * 0.2, y1=rep(0.3, 6), name="ticks")
    }else{
      labels <- textGrob(c("weak LD\nor\nnot informative", "strong LD"), x=c(0.25, 0.75), y=0.35,
                         gp=gpar(cex=0.3,lineheight=0.8), name="labels",just=c("centre","top"))
      ticks=nullGrob()
    }
    box <- rectGrob(x=0, y=0.7, height=0.3, width=1,gp=gpar(col="black",fill=  trascolor), name="box", just=c("left","top"))
    key <- gTree(children=gList(colorrect, title, labels,ticks, box), name="Key")
    key
  }

  keyVP <- viewport(x=0.02, y=0.98, height=unit(0.1, "snpc"), width=unit(0.09, "snpc"), ######need modify
                    just=c("left", "top"), name="keyVP")

  ## mapname--------------
  vpmapname <- viewport(x=0.095, y=0, height=unit(1, "snpc"), width=unit(0.1, "snpc"), ######need modify
                        just=c("right", "bottom"), name="keyVP")
  text_heatmap <- textGrob("LD Heat Map", x=1, y=0.75, just=c("right", "centre"), gp=gpar(cex=0.5),  vp=vpmapname)
  if(CLQshow == TRUE){
    text_clq <- textGrob("CLQ results", x=1, y=bottomloca, just=c("right", "top"),gp=gpar(cex=0.5),  vp=vpmapname)
  }else{
    text_clq <- nullGrob(vp=vpmapname)
  }

  if(onlyHeatmap == FALSE){
    blockregionname <- ifelse(blocktype == "bigld", "block regions\n(Big-LD)", "block regions\n(GPART)")
    text_LD <- textGrob(blockregionname, x=1, y=bottomloca-(CLQregionheight+mappingheight), just=c("right", "centre"),
                        gp=gpar(cex=0.5), vp=vpmapname)
  }else{
    text_LD <- textGrob("SNP locations", x=1, y=bottomloca-(CLQregionheight+mappingheight), just=c("right", "top"),
                        gp=gpar(cex=0.5), vp=vpmapname)
  }
  if(dim(nowgeneinfo)[1]>0){
    text_gene <- textGrob("Genes", x=1, y= (bottomloca-(CLQregionheight+mappingheight+locabarheight)-0.02),
                          just=c("right", "top"), gp=gpar(cex=0.5),
                          vp=vpmapname)
  }else{
    text_gene <-gList()
  }



  mapname <- gTree(children=gList(text_heatmap,text_LD, text_clq, text_gene))

  # --------------------------------------------end ------
  filename <- filename
  if(unitleng == 1.1){
    widthleng <- unitleng *1.1
  }else  if(unitleng == 2){
    widthleng <- unitleng *1.1
  }else  if(unitleng == 3){
    widthleng <- unitleng
  }else  if(unitleng == 4){
    widthleng <- unitleng
  }

  if((NSNPs > 200) & (onlyHeatmap == FALSE)){
    SNPmappingLine <- gList()
    # SNPsonBar <- gList()
  }

  if((NSNPs <= 200) & (onlyHeatmap == FALSE)){
    mappingLine <- gList()
    # SNPsonBar <- gList()
  }

  cat("\n")
  message("Generating file...." )
  if(type == "tif"){
    filename <- paste(filename, ".tif", sep="")
    tiff(filename,res=res, units="cm", height=15, width=15*widthleng)
    grid.draw(mappingLine)
    grid.draw(SNPmappingLine)
    grid.draw(rotLDtest)
    grid.draw(CLQboxgL)#
    grid.draw(LDSNPticksBox)#
    grid.draw(Locatick)
    grid.draw(Locatick1)
    grid.draw(locaBar)
    grid.draw(Locatickname)
    grid.draw(blockonBar)
    grid.draw(SNPsonBar)
    grid.draw(geneMap)
    grid.draw(LDSNPticks)#
    grid.draw(SNPblockmapping)#
    grid.draw(SNPblockmapping1)#
    grid.draw(LDblockticks)#
    grid.draw(LDblocktickslocaBar)#
    grid.draw(mapname)

    pushViewport(keyVP)
    grid.draw(LDheatmapLegend(color, LD))
    upViewport()

    dev.off()
    message(filename)
  }else if(type == "png"){
    filename <- paste(filename, ".png", sep="")
    png(filename,res=res, units="cm", height=15, width=15*widthleng)
    grid.draw(mappingLine)
    grid.draw(SNPmappingLine)
    grid.draw(rotLDtest)
    grid.draw(CLQboxgL)#
    grid.draw(LDSNPticksBox)#
    grid.draw(Locatick)
    grid.draw(Locatick1)
    grid.draw(locaBar)
    grid.draw(Locatickname)
    grid.draw(blockonBar)
    grid.draw(SNPsonBar)
    grid.draw(geneMap)
    grid.draw(LDSNPticks)#
    grid.draw(SNPblockmapping)#
    grid.draw(SNPblockmapping1)#
    grid.draw(LDblockticks)#
    grid.draw(LDblocktickslocaBar)#
    grid.draw(mapname)

    pushViewport(keyVP)
    grid.draw(LDheatmapLegend(color, LD))
    upViewport()

    dev.off()
    message(filename)
  }
}

# data ------------
#' @title genotype data
#' @name geno
#' @docType data
#' @aliases geno
#' @description This data set gives genotype data that are subset of 1000 Genomes Project phase1 release 3 genotype data
#' with 286 individuals from JPT, CHB and CHS populations (1000 Genomes Project Consortium, 2012).
#' The dataset contains 9000 SNPs in chromosome 21
#' @usage data(geno)
#' @format A data frame with 286 rows and 9000 columns
#' @source 1000 genomes project Phase 1 dataset <http://www.internationalgenome.org/>
#' @references {
#'  1000 Genomes Project Consortium.
#'  "An integrated map of genetic variation from 1,092 human genomes." \emph{Nature} 491.7422 (2012): 56.
#'  }
#' @keywords datasets
NULL

#' @title SNP information data
#' @name SNPinfo
#' @docType data
#' @aliases SNPinfo
#' @description This data set gives information data of SNPs in \code{geno}
#' The dataset contains chromosome name, rsID and bp information of the 9000 SNPs in \code{geno}
#' @usage data(SNPinfo)
#' @format A data frame with 9000 rows and 3 columns for chromosome name, rsID and bp
#' @source 1000 genomes project Phase 1 dataset <http://www.internationalgenome.org/>
#' @references {
#'  1000 Genomes Project Consortium.
#'  "An integrated map of genetic variation from 1,092 human genomes." \emph{Nature} 491.7422 (2012): 56.
#'  }
#'  @keywords datasets
NULL

#' @title gene information data
#' @name geneinfo
#' @docType data
#' @aliases geneinfo
#' @description This data set gives gene information in chromosome 21.
#' @usage data(geneinfo)
#' @format A data frame with 736 rows and 4 columns for genename, chromosome name, start bp and end bp of each gene.
#' @source <http://grch37.ensembl.org/index.html>
#' @references {
#' Zerbino, Daniel R., et al. "Ensembl 2018." \emph{Nucleic acids research 46.D1} (2017): D754-D761.
#' }
#'
#' @keywords datasets
NULL








