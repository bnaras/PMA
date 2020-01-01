# March 10 2009 - This function does sparse multiple CCA as described in Witten & Tibshirani (2009) extensions to sparse CCA method.


UpdateW <- function(xlist, i, K, sumabsthis, ws, type="standard", ws.final){
  tots <- 0
  for(j in (1:K)[-i]){
    diagmat <- (t(ws.final[[i]])%*%t(xlist[[i]]))%*%(xlist[[j]]%*%ws.final[[j]])
    diagmat[row(diagmat)!=col(diagmat)] <- 0
    tots <- tots + t(xlist[[i]])%*%(xlist[[j]]%*%ws[[j]]) - ws.final[[i]]%*%(diagmat%*%(t(ws.final[[j]])%*%ws[[j]]))
  }
  if(type=="standard"){
    sumabsthis <- BinarySearch(tots, sumabsthis)
    w <- soft(tots, sumabsthis)/l2n(soft(tots, sumabsthis))
  } else {
    tots <- as.numeric(tots)
    tots <- tots/mean(abs(tots))
    w <- FLSA(tots,lambda1=sumabsthis,lambda2=sumabsthis)[1,1,]
#    flsa.out <- diag.fused.lasso.new(tots,lam1=sumabsthis)
#    lam2ind <- which.min(abs(flsa.out$lam2-sumabsthis))
#    w <- flsa.out$coef[,lam2ind]
    w <- w/l2n(w)
    w[is.na(w)] <- 0
  }
  return(w)
}

GetCrit <- function(xlist, ws, K){
  crit <- 0
  for(i in 2:K){
    for(j in 1:(i-1)){
      crit <- crit + t(ws[[i]])%*%t(xlist[[i]])%*%xlist[[j]]%*%ws[[j]]
    }
  }
  return(crit)
}

GetCors <- function(xlist, ws, K){
  cors <- 0
  for(i in 2:K){
    for(j in 1:(i-1)){
      thiscor  <-  cor(xlist[[i]]%*%ws[[i]], xlist[[j]]%*%ws[[j]])
      if(is.na(thiscor)) thiscor <- 0
      cors <- cors + thiscor
    }
  }
  return(cors)
}


ftrans <- function(x){ return(.5*log((1+x)/(1-x))) }





#' Select tuning parameters for sparse multiple canonical correlation analysis
#' using the penalized matrix decomposition.
#'
#' This function can be used to automatically select tuning parameters for
#' sparse multiple CCA. This is the analog of sparse CCA, when >2 data sets are
#' available. Each data set may have features of type="standard" or
#' type="ordered" (e.g. CGH data). Assume that there are K data sets, called
#' $X1,...,XK$.
#'
#' The tuning parameters are selected using a permutation scheme. For each
#' candidate tuning parameter value, the following is performed: (1) Repeat the
#' following n times, for n large: (a) The samples in $(X1,...,XK)$ are
#' randomly permuted to obtain data sets $(X1*,...,XK*)$. (b) Sparse multiple
#' CCA is run on the permuted data sets $(X1*,...,XK*)$ to get canonical
#' variates $(w1*,...,wK*)$. (c) Record $t* = sum_(i<j) Cor(Xi* wi*, Xj* wj*)$.
#' (2) Sparse CCA is run on the original data $(X1,...,XK)$ to obtain canonical
#' variates $(w1,...,wK)$. (3) Record $t = sum_(i<j) Cor(Xi wi, Xj wj)$. (4)
#' The resulting p-value is given by $mean(t* > t)$; that is, the fraction of
#' permuted totals that exceed the total on the real data. Then, choose the
#' tuning parameter value that gives the smallest value in Step 4.
#'
#' This function only selets tuning parameters for the FIRST sparse multiple
#' CCA factors.
#'
#' Note that $x1,...,xK$ must have same number of rows. This function performs
#' just a one-dimensional search in tuning parameter space.
#'
#' @aliases MultiCCA.permute print.MultiCCA.permute plot.MultiCCA.permute
#' @param xlist A list of length K, where K is the number of data sets on which
#' to perform sparse multiple CCA. Data set k should be a matrix of dimension
#' $n x p_k$ where $p_k$ is the number of features in data set k.
#' @param penalties The penalty terms to be considered in the cross-validation.
#' If the same penalty term is desired for each data set, then this should be a
#' vector of length equal to the number of penalty terms to be considered. If
#' different penalty terms are desired for each data set, then this should be a
#' matrix with rows equal to the number of data sets, and columns equal to the
#' number of penalty terms to be considered. For a given data set Xk, if type
#' is "standard" then the penalty term should be a number between 1 and
#' $sqrt(p_k)$ (the number of features in data set k); it is a L1 bound on wk.
#' If type is "ordered", on the other hand, the penalty term is of the form
#' lambda in the fused lasso penalty. Therefore, the interpretation of the
#' argument depends on whether type is "ordered" or "standard" for this data
#' set.
#' @param type A K-vector containing elements "standard" or "ordered" - or a
#' single value. If a single value, then it is assumed that all elements are
#' the same (either "standard" or "ordered").  If columns of v are ordered
#' (e.g. CGH spots ordered along the chromosome) then "ordered", otherwise use
#' "standard". "standard" will result in a lasso ($L_1$) penalty on v, which
#' will result in smoothness. "ordered" will result in a fused lasso penalty on
#' v, yielding both sparsity and smoothness.
#' @param niter How many iterations should be performed each time CCA is
#' called? Default is 3, since an approximate estimate of u and v is acceptable
#' in this case, and otherwise this function can be quite time-consuming.
#' @param ws A list of length K; the kth element contanis the first ncomponents
#' columns of the v matrix of the SVD of Xk. If NULL, then the SVD of Xk will
#' be computed inside this function. However, if you plan to run this function
#' multiple times, then save a copy of this argument so that it does not need
#' to be re-computed.
#' @param trace Print out progress?
#' @param nperms How many times should the data be permuted? Default is 25. A
#' large value of nperms is very important here, since the formula for
#' computing the z-statistics requires a standard deviation estimate for the
#' correlations obtained via permutation, which will not be accurate if nperms
#' is very small.
#' @param standardize Should the columns of X and Z be centered (to have mean
#' zero) and scaled (to have standard deviation 1)? Default is TRUE.
#' @param x not used.
#' @param \dots not used.
#' @return \item{zstat}{The vector of z-statistics, one per element of
#' penalties.} \item{pvals}{The vector of p-values, one per element of
#' penalties.} \item{bestpenalties}{The best set of penalties (the one with the
#' highest zstat).} \item{cors}{The value of $sum_(j<k) cor(Xk wk, Xj wj)$
#' obtained for each value of penalties.} \item{corperms}{The nperms values of
#' $sum_(j<k) cor(Xk* wk*, Xj* wj*)$ obtained for each value of penalties,
#' where Xk* indicates the Xk matrix with permuted rows, and wk* is the
#' canonical variate corresponding to the permuted data.}
#' \item{ws.init}{Initial values used for ws in sparse multiple CCA algorithm.}
#' @seealso \link{MultiCCA}, \link{CCA.permute}, \link{CCA}
#' @references
#' \insertRef{pmid19377034}{PMA}
#' @importFrom graphics lines par plot points segments text title
#' @importFrom stats cor lowess lsfit quantile rnorm runif sd
#' @examples
#'
#' # See examples in MultiCCA function
#'
#' @export MultiCCA.permute
MultiCCA.permute <- function(xlist, penalties=NULL, ws=NULL, type="standard", nperms=10, niter=3, trace=TRUE, standardize=TRUE){
  call <- match.call()
  K <- length(xlist)
  for(k in 1:K){
    if(ncol(xlist[[k]])<2) stop("Need at least 2 features in each data set!")
    if(standardize) xlist[[k]] <- scale(xlist[[k]], T, T)
  }
  if(length(type)==1) type <- rep(type, K) # If type is just a single element, expand to make a vector of length(xlist)
          # Or type can have standard/ordered for each elt of xlist
  if(length(type)!=K) stop("Type must be a vector of length 1, or length(xlist)")
  if(sum(type!="standard" & type!="ordered")>0) stop("Each element of type must be standard or ordered.")
  if(is.null(penalties)){
    if(sum(type=="ordered")==K) stop("Do not run MultiCCA.permute with only ordered data sets and penalties unspecified,
                                      since we only choose tuning the parameter via permutations when type='standard'.")
    penalties <- matrix(NA, nrow=K, ncol=10)
    for(k in 1:K){
      if(type[k]=="ordered"){
        lam <- ChooseLambda1Lambda2(svd(xlist[[k]])$v[,1])
        penalties[k,] <- lam
      } else {
        penalties[k,] <- pmax(seq(.1, .8, len=10)*sqrt(ncol(xlist[[k]])),1.1)
      }
    }
  }
  numnonzeros <- NULL
  if(!is.matrix(penalties)) penalties <- matrix(1,nrow=K,ncol=1)%*%matrix(penalties,nrow=1)
  permcors <- matrix(NA, nrow=nperms, ncol=ncol(penalties))
  cors <- numeric(ncol(penalties))
  for(i in 1:ncol(penalties)){
    out <- MultiCCA(xlist, penalty=penalties[,i], niter=niter, type=type, ws=ws, trace=trace)
    cors[i] <- GetCors(xlist, out$ws, K)
    numnonzeros <- c(numnonzeros, sum(out$numnonzeros))
    ws.init  <- out$ws.init
  }
  cat(fill=TRUE)
  for(j in 1:nperms){
    if(trace) cat("Permutation ", j, "of " , nperms ,fill=TRUE)
    xlistperm <- xlist
    for(k in 1:K){
      xlistperm[[k]] <- xlistperm[[k]][sample(1:nrow(xlistperm[[k]])),]
    }
    for(i in 1:ncol(penalties)){
      out <- MultiCCA(xlistperm, penalty=penalties[,i], niter=niter, type=type, ws=ws, trace=FALSE)
      permcors[j,i] <- GetCors(xlistperm, out$ws, K)
    }
  }
  pvals =zs =  NULL
  for(i in 1:ncol(penalties)){
    pvals <- c(pvals, mean(permcors[,i]>=cors[i]))
    zs <- c(zs, (cors[i]-mean(permcors[,i]))/(sd(permcors[,i])+.05))
  }
  if(trace) cat(fill=TRUE)
  out <- list(pvals=pvals, zstat=zs, bestpenalties=penalties[,which.max(zs)], cors=cors, corperms=permcors, numnonzeros=numnonzeros, ws.init=ws.init, call=call, penalties=penalties, type=type, nperms=nperms)
  class(out) <- "MultiCCA.permute"
  return(out)
}





#' Perform sparse multiple canonical correlation analysis.
#'
#' Given matrices $X1,...,XK$, which represent K sets of features on the same
#' set of samples, find sparse $w1,...,wK$ such that $sum_(i<j) (wi' Xi' Xj
#' wj)$ is large. If the columns of Xk are ordered (and type="ordered") then wk
#' will also be smooth. For $X1,...,XK$, the samples are on the rows and the
#' features are on the columns. $X1,...,XK$ must have same number of rows, but
#' may (and usually will) have different numbers of columns.
#'
#'
#' @aliases MultiCCA print.MultiCCA
#' @param xlist A list of length K, where K is the number of data sets on which
#' to perform sparse multiple CCA. Data set k should be a matrix of dimension
#' $n x p_k$ where $p_k$ is the number of features in data set k.
#' @param penalty The penalty terms to be used. Can be a single value (if the
#' same penalty term is to be applied to each data set) or a K-vector,
#' indicating a different penalty term for each data set. There are 2 possible
#' interpretations for the penalty terms: If type="standard" then this is an L1
#' bound on wk, and it must be between 1 and $sqrt(p_k)$ ($p_k$ is the number
#' of features in matrix Xk). If type="ordered" then this is the parameter for
#' the fused lasso penalty on wk.
#' @param type Are the columns of $x1,...,xK$ unordered (type="standard") or
#' ordered (type="ordered")? If "standard", then a lasso penalty is applied to
#' v, to enforce sparsity. If "ordered" (generally used for CGH data), then a
#' fused lasso penalty is applied, to enforce both sparsity and smoothness.
#' This argument can be a vector of length K (if different data sets are of
#' different types) or it can be a single value "ordered"/"standard" (if all
#' data sets are of the same type).
#' @param ncomponents How many factors do you want? Default is 1.
#' @param niter How many iterations should be performed? Default is 25.
#' @param ws A list of length K. The kth element contains the first ncomponents
#' columns of the v matrix of the SVD of Xk. If NULL, then the SVD of
#' $X1,...,XK$ will be computed inside the MultiCCA function. However, if you
#' plan to run this function multiple times, then save a copy of this argument
#' so that it does not need to be re-computed.
#' @param trace Print out progress?
#' @param standardize Should the columns of $X1,...,XK$ be centered (to have
#' mean zero) and scaled (to have standard deviation 1)? Default is TRUE.
#' @param x not used.
#' @param \dots not used.
#' @return \item{ws}{A list of length K, containg the sparse canonical variates
#' found (element k is a $p_k x ncomponents$ matrix).} \item{ws.init}{A list of
#' length K containing the initial values of ws used, by default these are the
#' v vector of the svd of matrix Xk.}
#' @seealso \link{MultiCCA.permute},\link{CCA}, \link{CCA.permute}
#' @references
#' \insertRef{pmid19377034}{PMA}
#' @examples
#'
#' # Generate 3 data sets so that first 25 features are correlated across
#' # the data sets...
#' u <- matrix(rnorm(50),ncol=1)
#' v1 <- matrix(c(rep(.5,25),rep(0,75)),ncol=1)
#' v2 <- matrix(c(rep(1,25),rep(0,25)),ncol=1)
#' v3 <- matrix(c(rep(.5,25),rep(0,175)),ncol=1)
#'
#' x1 <- u%*%t(v1) + matrix(rnorm(50*100),ncol=100)
#' x2 <- u%*%t(v2) + matrix(rnorm(50*50),ncol=50)
#' x3 <- u%*%t(v3) + matrix(rnorm(50*200),ncol=200)
#'
#' xlist <- list(x1, x2, x3)
#'
#' # Run MultiCCA.permute w/o specifying values of tuning parameters to
#' # try.
#' # The function will choose the lambda for the ordered data set.
#' # Then permutations will be used to select optimal sum(abs(w)) for
#' # standard data sets.
#' # We assume that x1 is standard, x2 is ordered, x3 is standard:
#' perm.out <- MultiCCA.permute(xlist, type=c("standard", "ordered",
#' "standard"))
#' print(perm.out)
#' plot(perm.out)
#' out <- MultiCCA(xlist, type=c("standard", "ordered", "standard"),
#' penalty=perm.out$bestpenalties, ncomponents=2, ws=perm.out$ws.init)
#' print(out)
#' # Or if you want to specify tuning parameters by hand:
#' # this time, assume all data sets are standard:
#' perm.out <- MultiCCA.permute(xlist, type="standard",
#' penalties=cbind(c(1.1,1.1,1.1),c(2,3,4),c(5,7,10)), ws=perm.out$ws.init)
#' print(perm.out)
#' plot(perm.out)
#'
#' # Making use of the fact that the features are ordered:
#' out <- MultiCCA(xlist, type="ordered", penalty=.6)
#' par(mfrow=c(3,1))
#' PlotCGH(out$ws[[1]], chrom=rep(1,ncol(x1)))
#' PlotCGH(out$ws[[2]], chrom=rep(2,ncol(x2)))
#' PlotCGH(out$ws[[3]], chrom=rep(3,ncol(x3)))
#'
#' @export MultiCCA
MultiCCA <- function(xlist, penalty=NULL, ws=NULL, niter=25, type="standard", ncomponents=1, trace=TRUE, standardize=TRUE){
  for(i in 1:length(xlist)){
    if(ncol(xlist[[i]])<2) stop("Need at least 2 features in each data set.")
  }
  call <- match.call()
  K <- length(xlist)
  if(length(type)==1) type <- rep(type, K) # If type is just a single element, expand to make a vector of length(xlist)
          # Or type can have standard/ordered for each elt of xlist
  if(length(type)!=K) stop("Type must be a vector of length 1, or length(xlist)")
  if(sum(type!="standard" & type!="ordered")>0) stop("Each element of type must be standard or ordered.")
  for(k in 1:K){
    if(standardize) xlist[[k]] <- scale(xlist[[k]], T, T)
  }
  if(!is.null(ws)){
    makenull <- FALSE
    for(i in 1:K){
      if(ncol(ws[[i]])<ncomponents) makenull <- TRUE
    }
    if(makenull) ws <- NULL
  }
  if(is.null(ws)){
    ws <- list()
    for(i in 1:K) ws[[i]] <- matrix(svd(xlist[[i]])$v[,1:ncomponents], ncol=ncomponents)
  }
  if(is.null(penalty)){
    penalty <- rep(NA, K)
    penalty[type=="standard"] <- 4 # this is the default value of sumabs
    for(k in 1:K){
      if(type[k]=="ordered"){
        v <- svd(xlist[[k]])$v[,1]
        penalty[k] <- ChooseLambda1Lambda2(v)
      }
    }
  }
  ws.init <- ws
  if(length(penalty)==1) penalty <- rep(penalty, K)
  if(sum(penalty<1 & type=="standard")) stop("Cannot constrain sum of absolute values of weights to be less than 1.")
  for(i in 1:length(xlist)){
    if(type[i]=="standard" && penalty[i]>sqrt(ncol(xlist[[i]]))) stop("L1 bound of weights should be no more than sqrt of the number of columns of the corresponding data set.", fill=TRUE)
  }
  ws.final <- list()
  for(i in 1:length(ws)) ws.final[[i]] <- matrix(0, nrow=ncol(xlist[[i]]), ncol=ncomponents)
  cors <- NULL
  for(comp in 1:ncomponents){
    ws <- list()
    for(i in 1:length(ws.init)) ws[[i]] <- ws.init[[i]][,comp]
    curiter <- 1
    crit.old <- -10
    crit <- -20
    storecrits <- NULL
    while(curiter<=niter && abs(crit.old-crit)/abs(crit.old)>.001 && crit.old!=0){
      crit.old <- crit
      crit <- GetCrit(xlist, ws, K)
      storecrits <- c(storecrits,crit)
      if(trace) cat(curiter, fill=FALSE)
      curiter <- curiter+1
      for(i in 1:K){
        ws[[i]] <- UpdateW(xlist, i, K, penalty[i], ws, type[i], ws.final)
      }
    }
    for(i in 1:length(ws)) ws.final[[i]][,comp] <- ws[[i]]
    cors <- c(cors, GetCors(xlist, ws,K))
  }
  out <- list(ws=ws.final, ws.init=ws.init, K=K, call=call, type=type, penalty=penalty, cors=cors)
  class(out) <- "MultiCCA"
  return(out)
}

print.MultiCCA <- function(x,...){
  cat("Call: ")
  dput(x$call)
  cat("\n\n")
  cat("Sum_{i<j} Cor(Xi wi, Xj wj) = ", sep="", paste(round(x$cors,4), sep="", " "),fill=TRUE)
  cat("There are ", x$K, " data sets.", fill=TRUE)
  for(i in 1:x$K){
    cat("Data set ", i, " is of type ", x$type[i], ".", fill=TRUE)
    if(x$type[i]=="ordered") cat("Tuning parameter used: Lambda was ", round(x$penalty[i],4),fill=TRUE)
    if(x$type[i]=="standard") cat("Tuning parameter used: Sum(abs(w)) was ", round(x$penalty[i],4),fill=TRUE)
    cat("Num non-zero elements of canonical variate(s) for data set ", i, ":    ")
    if(is.matrix(x$ws[[i]])) cat(apply(x$ws[[i]]!=0, 2, sum), fill=TRUE)
    if(!is.matrix(x$ws[[i]])) cat(sum(x$ws[[i]]!=0), fill=TRUE)
    cat(fill=TRUE)
  }
}




print.MultiCCA.permute <- function(x,...){
  cat("Call: ")
  dput(x$call)
  cat("\n\n")
  tab <- round(cbind(x$pvals, x$zstat, x$cors, colMeans(x$corperms)), 3)
  dimnames(tab) <- list(paste("Tuning parameter set ", sep="", 1:length(x$pvals)), c("P-Value", "Z", "Cors", "Cors Perm"))
  print(tab, quote=FALSE)
  cat("Highest z score: ", max(x$zstat), "\n")
  cat("P-value corresponding to highest z score: ", x$pvals[which.max(x$zstat)], fill=TRUE)
  cat("Tuning parameters corresponding to highest z score: ", round(x$bestpenalties,3), "\n")
  sumabslamvecs <- round(x$penalties,4)
  dimnames(sumabslamvecs) <- list(paste(paste("Data set", sep=" ", 1:nrow(sumabslamvecs)),
                                        sep="", paste(paste("; Type is ", sep="", x$type),sep="",": ")), 1:ncol(sumabslamvecs))
  cat(fill=TRUE)
  cat("Tuning parameters used: ",fill=TRUE)
  print(sumabslamvecs,quote=FALSE)
}

plot.MultiCCA.permute <- function(x,...){
  sumabss <- x$penalties
  ccs <- x$cors
  nperms <- x$nperms
  zstats <- x$zstat
  ccperms <- x$corperms
  par(mfrow=c(2,1))
  plot(1:ncol(sumabss), ccs, main="Correlations For Real/Permuted Data", xlab="Index of Tuning Parameter Set",
       ylab="Correlations", ylim=range(ccperms,ccs))
  points(1:ncol(sumabss),ccs,type="l")
  for(i in 1:nperms){
    points(1:ncol(sumabss),ccperms[i,],col="green")
  }
  plot(1:ncol(sumabss),zstats,main="Z", xlab="Index of Tuning Parameter Set", ylab="Z score")
  lines(1:ncol(sumabss),zstats)
}
