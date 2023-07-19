soft <- function(x,d){
  return(sign(x)*pmax(0, abs(x)-d))
}

mean_na <- function(vec){
  return(mean(vec[!is.na(vec)]))
}

l2n <- function(vec){
  a <- sqrt(sum(vec^2))
  if(a==0) a <- .05
  return(a)
}

safesvd <- function(x){
  i <- 1
  out <- try(svd(x), silent=TRUE)
  while(i<10 && inherits(out, "try-error")){
    out <- try(svd(x), silent=TRUE)
    i <- i+1
  }
  if(inherits(out, "try-error")) out <- svd(matrix(rnorm(nrow(x)*ncol(x)), ncol=ncol(x)))
  return(out)
}

BinarySearch <- function(argu,sumabs){
  if(l2n(argu)==0 || sum(abs(argu/l2n(argu)))<=sumabs) return(0)
  lam1 <- 0
  lam2 <- max(abs(argu))-1e-5
  iter <- 1
  while(iter < 150){
    su <- soft(argu,(lam1+lam2)/2)
    if(sum(abs(su/l2n(su)))<sumabs){
      lam2 <- (lam1+lam2)/2
    } else {
      lam1 <- (lam1+lam2)/2
    }
    if((lam2-lam1)<1e-6) return((lam1+lam2)/2)
    iter <- iter+1
  }
  warning("Didn't quite converge")
  return((lam1+lam2)/2)
}



SMD <- function(x, sumabsu, sumabsv, niter=20,trace=TRUE, v, upos, uneg, vpos, vneg){
  # This gets a single factor. Do MultiSMD to get multiple factors.
  nas <- is.na(x)
  v.init <- v
  xoo <- x
  if(sum(nas)>0) xoo[nas] <- mean(x[!nas])
  oldv <- rnorm(ncol(x))
  for(iter in 1:niter){
    if(sum(abs(oldv-v))>1e-7){
      oldv <- v
      if(trace) cat(iter,fill=F)
      # update u #
      argu <- xoo%*%v
      if(upos) argu <- pmax(argu,0)
      if(uneg) argu <- pmin(argu,0)
      lamu <- BinarySearch(argu,sumabsu)
      su <- soft(argu,lamu)
      u <- matrix(su/l2n(su),ncol=1)
      # done updating u #
      # update v #
      argv <- t(u)%*%xoo
      if(vpos) argv <- pmax(argv,0)
      if(vneg) argv <- pmin(argv,0)
      lamv <- BinarySearch(argv, sumabsv)
      sv <- soft(argv,lamv)
      v <- matrix(sv/l2n(sv),ncol=1)
      # done updating v #
    }
  }
  d <- as.numeric(t(u)%*%(xoo%*%v))
  if(trace) cat(fill=TRUE)
  return(list(d=d, u=u, v=v, v.init=v.init))
}

SMDOrth <- function(x, us, sumabsv=NULL, niter=20, trace=TRUE,v, vpos, vneg){
  # Gets the next u for sparse PCA, using Trevor's method of requiring this new u to be orthog
  # to all previous us (contained in us)
  nas <- is.na(x)
  v.init <- v
  xoo <- x
  if(sum(nas)>0) xoo[nas] <- mean(x[!nas])
  u <- rnorm(nrow(x))
  oldu <- rnorm(nrow(x))
  oldv <- rnorm(ncol(x))
  for(iter in 1:niter){
    if(sum(abs(oldu-u))>1e-6 || sum(abs(oldv-v))>1e-6){
      oldu <- u
      oldv <- v
      if(trace) cat(iter,fill=F)
      # update u #
      argu <- xoo%*%v
      numer <- lsfit(y=argu, x=us, intercept=FALSE)$res
      u <- numer/l2n(numer)
      # done updating u #
      # update v #
      argv <- t(u)%*%xoo
      if(vpos) argv <- pmax(argv,0)
      if(vneg) argv <- pmin(argv,0)
      lamv <- BinarySearch(argv, sumabsv)
      sv <- soft(argv,lamv)
      v <- matrix(sv/l2n(sv),ncol=1)
      # done updating v #
    }
  }
  d <- as.numeric(t(u)%*%(xoo%*%v))
  if(trace) cat(fill=TRUE)
  return(list(d=d,u=u,v=v, v.init=v.init))
}


CheckPMDV <-  function(v,x,K){
  if(!is.null(v) && is.matrix(v) && ncol(v)>=K){
    v <- matrix(v[,1:K], ncol=K)
  } else if(ncol(x)>nrow(x)){
    x[is.na(x)] <- mean_na(x)
    v <- matrix(t(x)%*%(safesvd(x%*%t(x))$v[,1:K]),ncol=K)
    if(sum(is.na(v))>0) v <- matrix(safesvd(x)$v[,1:K], ncol=K)
    v <- sweep(v,2,apply(v, 2, l2n), "/")
    if(sum(is.na(v))>0) stop("some are NA")
  } else if (ncol(x)<=nrow(x)){
    x[is.na(x)] <- mean_na(x)
    v <- matrix(safesvd(t(x)%*%x)$v[,1:K],ncol=K)
  }
  return(v)
}





#' Perform sparse principal component analysis
#'
#' Performs sparse principal components analysis by applying PMD to a data
#' matrix with lasso ($L_1$) penalty on the columns and no penalty on the rows.
#'
#' PMD(x,sumabsu=sqrt(nrow(x)), sumabsv=3, K=1) and SPC(x,sumabsv=3, K=1) give
#' the same result, since the SPC method is simply PMD with an L1 penalty on
#' the columns and no penalty on the rows.
#'
#' In Witten, Tibshirani, and Hastie (2008), two methods are presented for
#' obtaining multiple factors for SPC. The methods are as follows:
#'
#' (1) If one has already obtained factors $k-1$ factors then oen can compute
#' residuals by subtracting out these factors. Then $u_k$ and $v_k$ can be
#' obtained by applying the SPC/PMD algorithm to the residuals.
#'
#' (2) One can require that $u_k$ be orthogonal to $u_i$'s with $i<k$; the
#' method is slightly more complicated, and is explained in WT&H(2008).
#'
#' Method 1 is performed by running SPC with option orth=FALSE (the default)
#' and Method 2 is performed using option orth=TRUE. Note that Methods 1 and 2
#' always give identical results for the first component, and often given quite
#' similar results for later components.
#'
#' @aliases SPC print.SPC
#' @param x Data matrix of dimension $n x p$, which can contain NA for missing
#' values. We are interested in finding sparse principal components of
#' dimension $p$.
#' @param sumabsv How sparse do you want v to be? This is the sum of absolute
#' values of elements of v. It must be between 1 and square root of number of
#' columns of data. The smaller it is, the sparser v will be.
#' @param niter How many iterations should be performed. It is best to run at
#' least 20 of so. Default is 20.
#' @param K The number of factors in the PMD to be returned; default is 1.
#' @param v The first right singular vector(s) of the data. (If missing data is
#' present, then the missing values are imputed before the singular vectors are
#' calculated.) v is used as the initial value for the iterative PMD($L_1$,
#' $L_1$) algorithm. If x is large, then this step can be time-consuming;
#' therefore, if PMD is to be run multiple times, then v should be computed
#' once and saved.
#' @param trace Print out progress as iterations are performed? Default is
#' TRUE.
#' @param orth If TRUE, then use method of Section 3.2 of Witten, Tibshirani
#' and Hastie (2008) to obtain multiple sparse principal components. Default is
#' FALSE.
#' @param center Subtract out mean of x? Default is TRUE
#' @param cnames An optional vector containing a name for each column.
#' @param vpos Constrain the elements of v to be positive? TRUE or FALSE.
#' @param vneg Constrain the elements of v to be negative? TRUE or FALSE.
#' @param compute.pve Compute percent variance explained? Default TRUE. If not
#' needed, then choose FALSE to save time.
#' @return \item{u}{u is output. If you asked for multiple factors then each
#' column of u is a factor. u has dimension nxK if you asked for K factors.}
#' \item{v}{v is output. These are the sparse principal components. If you
#' asked for multiple factors then each column of v is a factor. v has
#' dimension pxK if you asked for K factors.} \item{d}{d is output; it is the
#' diagonal of the matrix $D$ in the penalized matrix decomposition. In the
#' case of the rank-1 decomposition, it is given in the formulation
#' $||X-duv'||_F^2$ subject to $||u||_1 <= sumabsu$, $||v||_1 <= sumabsv$.
#' Computationally, $d=u'Xv$ where $u$ and $v$ are the sparse factors output by
#' the PMD function and $X$ is the data matrix input to the PMD function.}
#' \item{prop.var.explained}{A vector containing the proportion of variance
#' explained by the first 1, 2, ..., K sparse principal components obtaineds.
#' Formula for proportion of variance explained is on page 20 of Shen & Huang
#' (2008), Journal of Multivariate Analysis 99: 1015-1034.} \item{v.init}{The
#' first right singular vector(s) of the data; these are returned to save on
#' computation time if PMD will be run again.} \item{meanx}{Mean of x that was
#' subtracted out before SPC was performed.}
#' @seealso \link{SPC.cv}, \link{PMD}, \link{PMD.cv}
#' @references
#' Witten D. M., Tibshirani R.,  and Hastie, T. (2009)
#' \emph{A penalized matrix decomposition, with applications to sparse principal components and canonical correlation analysis}, \emph{Biostatistics, Gol 10 (3), 515-534, Jul 2009}\cr
#' @examples
#'
#'
#' # A simple simulated example
#' #NOT RUN
#' #set.seed(1)
#' #u <- matrix(c(rnorm(50), rep(0,150)),ncol=1)
#' #v <- matrix(c(rnorm(75),rep(0,225)), ncol=1)
#' #x <- u%*%t(v)+matrix(rnorm(200*300),ncol=300)
#' ## Perform Sparse PCA - that is, decompose a matrix w/o penalty on rows
#' ## and w/ L1 penalty on columns
#' ## First, we perform sparse PCA and get 4 components, but we do not
#' ## require subsequent components to be orthogonal to previous components
#' #out <- SPC(x,sumabsv=3, K=4)
#' #print(out,verbose=TRUE)
#' ## We could have selected sumabsv by cross-validation, using function SPC.cv
#' ## Now, we do sparse PCA using method in Section 3.2 of WT&H(2008) for getting
#' ## multiple components - that is, we require components to be orthogonal
#' #out.orth <- SPC(x,sumabsv=3, K=4, orth=TRUE)
#' #print(out.orth,verbose=TRUE)
#' #par(mfrow=c(1,1))
#' #plot(out$u[,1], out.orth$u[,1], xlab="", ylab="")
#' ## Note that the first components w/ and w/o orth option are identical,
#' ## since the orth option only affects the way that subsequent components
#' ## are found
#' #print(round(t(out$u)%*%out$u,4)) # not orthogonal
#' #print(round(t(out.orth$u)%*%out.orth$u,4)) # orthogonal
#' #
#' ## Use SPC.cv to choose tuning parameters:
#' #cv.out <- SPC.cv(x)
#' #print(cv.out)
#' #plot(cv.out)
#' #out <- SPC(x, sumabsv=cv.out$bestsumabsv)
#' #print(out)
#' ## or we could do
#' #out <- SPC(x, sumabsv=cv.out$bestsumabsv1se)
#' #print(out)
#' #
#' #
#' @export SPC
SPC <- function(x, sumabsv=4, niter=20, K=1, orth=FALSE, trace=TRUE, v=NULL, center=TRUE, cnames=NULL, vpos=FALSE, vneg=FALSE, compute.pve=TRUE){
  if(vpos&&vneg) stop("Cannot constrain elements to be positive AND negative.")
  out <- PMDL1L1(x,sumabsu=sqrt(nrow(x)), sumabsv=sumabsv, niter=niter,K=K,orth=orth,trace=trace,v=v,center=center,cnames=cnames, upos=FALSE, uneg=FALSE, vpos=vpos, vneg=vneg)
  if(compute.pve){
     # Calculate percent variance explained, using definition on page 7 of Shen and Huang (2008) Journal of Multivariate Analysis vol 99: 1015-1034
    v <- matrix(out$v, ncol=K)
    ve <- NULL # variance explained
    xfill <- x
    if(center) xfill <- x-mean_na(x)
    xfill[is.na(x)] <- mean_na(xfill)
    for(k in 1:K){
      vk <- matrix(v[,1:k], ncol=k)
      xk <- xfill%*%vk%*%solve(t(vk)%*%vk)%*%t(vk)
      svdxk <- svd(xk)
      ve <- c(ve, sum(svdxk$d^2))
    }
    pve <- ve/sum(svd(xfill)$d^2) # proportion of variance explained
    out$prop.var.explained <- pve
  }
  out$vpos <- vpos
  out$vneg <- vneg
  class(out) <- "SPC"
  return(out)
}

#' @method print SPC
#' @export
print.SPC <- function(x,verbose=FALSE,...){
  cat("Call: ")
  dput(x$call)
  cat("\n\n")
  cat("Num non-zero v's: ", apply(x$v!=0, 2, sum), "\n")
  cat("Sumabsv used: ", x$sumabsv, "\n")
  cat("Proportion of variance explained by first k components: \n")
  tab <- cbind(paste(paste("First ", sep="", 1:length(x$prop.var.explained)), sep="", " Components"), round(x$prop.var.explained, 3))
  dimnames(tab) <- list(1:length(x$prop.var.explained), c("", "Prop. Var. Explained"))
  print(tab, quote=FALSE, sep="\t")
  if(!is.null(x$meanx)) cat("Subtracted out mean of x: ", x$meanx, "\n")
  if(verbose){
    for(k in 1:x$K){
      cat("\n Component ", k, ":\n")
      v <- x$v[,k]
      if(is.null(x$cnames)) x$cnames <- 1:length(v)
      cat(fill=T)
      vs <- cbind(x$cnames[v!=0], round(v[v!=0],3))
      dimnames(vs) <- list(1:sum(v!=0), c("Feature Name", "Feature Weight"))
      print(vs, quote=FALSE, sep="\t")
    }
  }
  if(x$vpos) cat("Elements of v constrained to be positive.", fill=TRUE)
  if(x$vneg) cat("Elements of v constrained to be negative.", fill=TRUE)
}





#' Get a penalized matrix decomposition for a data matrix.
#'
#' Performs a penalized matrix decomposition for a data matrix. Finds factors u
#' and v that summarize the data matrix well. u and v will both be sparse, and
#' v can optionally also be smooth.
#'
#' The criterion for the PMD is as follows: we seek vectors $u$ and $v$ such
#' that $u'Xv$ is large, subject to $||u||_2=1, ||v||_2=1$ and additional
#' penalties on $u$ and $v$. These additional penalties are as follows: If type
#' is "standard", then lasso ($L_1$) penalties (promoting sparsity) are placed
#' on u and v. If type is "ordered", then lasso penalty is placed on u and a
#' fused lasso penalty (promoting sparsity and smoothness) is placed on v.
#'
#' If type is "standard", then arguments sumabs OR sumabsu&sumabsv are used. If
#' type is "ordered", then sumabsu AND lambda are used. Sumabsu is the bound of
#' absolute value of elements of u. Sumabsv is bound of absolute value of
#' elements of v. If sumabs is given, then sumabsu is set to
#' sqrt(nrow(x))*sumabs and sumabsv is set to sqrt(ncol(x))*sumabs. $lambda$ is
#' the parameter for the fused lasso penalty on v when type is "ordered":
#' $lambda(||v||_1 + sum_j |v_j - v_(j-1))$.
#'
#' @param x Data matrix of dimension $n x p$, which can contain NA for missing
#' values.
#' @param type "standard" or "ordered": Do we want v to simply be sparse, or
#' should it also be smooth? If the columns of x are ordered (e.g. CGH spots
#' along a chromosome) then choose "ordered". Default is "standard". If
#' "standard", then the PMD function will make use of sumabs OR
#' sumabsu&sumabsv. If "ordered", then the function will make use of sumabsu
#' and lambda.
#' @param sumabs Used only if type is "standard". A measure of sparsity for u
#' and v vectors, between 0 and 1. When sumabs is specified, and sumabsu and
#' sumabsv are NULL, then sumabsu is set to $sqrt(n)*sumabs$ and sumabsv is set
#' to $sqrt(p)*sumabs$. If sumabs is specified, then sumabsu and sumabsv should
#' be NULL. Or if sumabsu and sumabsv are specified, then sumabs should be
#' NULL.
#' @param sumabsu Used for types "ordered" AND "standard". How sparse do you
#' want u to be? This is the sum of absolute values of elements of u. It must
#' be between 1 and the square root of the number of rows in data matrix. The
#' smaller it is, the sparser u will be.
#' @param sumabsv Used only if type is "standard". How sparse do you want v to
#' be? This is the sum of absolute values of elements of v. It must be between
#' 1 and square root of number of columns of data. The smaller it is, the
#' sparser v will be.
#' @param lambda Used only if type is "ordered". This is the tuning parameter
#' for the fused lasso penalty on v, which takes the form $lambda ||v||_1 +
#' lambda |v_j - v_(j-1)|$. $lambda$ must be non-negative. If NULL, then it is
#' chosen adaptively from the data.
#' @param niter How many iterations should be performed. It is best to run at
#' least 20 of so. Default is 20.
#' @param K The number of factors in the PMD to be returned; default is 1.
#' @param v The first right singular vector(s) of the data. (If missing data is
#' present, then the missing values are imputed before the singular vectors are
#' calculated.) v is used as the initial value for the iterative PMD algorithm.
#' If x is large, then this step can be time-consuming; therefore, if PMD is to
#' be run multiple times, then v should be computed once and saved.
#' @param trace Print out progress as iterations are performed? Default is
#' TRUE.
#' @param center Subtract out mean of x? Default is TRUE.
#' @param chrom If type is "ordered", then this gives the option to specify
#' that some columns of x (corresponding to CGH spots) are on different
#' chromosomes. Then v will be sparse, and smooth *within* each chromosome but
#' not *between* chromosomes. Length of chrom should equal number of columns of
#' x, and each entry in chrom should be a number corresponding to which
#' chromosome the CGH spot is on.
#' @param rnames An optional vector containing a name for each row of x.
#' @param cnames An optional vector containing a name for each column of x.
#' @param upos Constrain the elements of u to be positive? TRUE or FALSE.
#' @param uneg Constrain the elements of u to be negative? TRUE or FALSE.
#' @param vpos Constrain the elements of v to be positive? TRUE or FALSE.
#' Cannot be used if type is "ordered".
#' @param vneg Constrain the elements of v to be negative? TRUE or FALSE.
#' Cannot be used if type is "ordered."
#' @return \item{u}{u is output. If you asked for multiple factors then each
#' column of u is a factor. u has dimension nxK if you asked for K factors.}
#' \item{v}{v is output. If you asked for multiple factors then each column of
#' v is a factor. v has dimension pxK if you asked for K factors.} \item{d}{d
#' is output. Computationally, $d=u'Xv$ where $u$ and $v$ are the sparse
#' factors output by the PMD function and $X$ is the data matrix input to the
#' PMD function. When K=1, the residuals of the rank-1 PMD are given by $X -
#' duv'$.} \item{v.init}{The first right singular vector(s) of the data; these
#' are returned to save on computation time if PMD will be run again.}
#' \item{meanx}{Mean of x that was subtracted out before PMD was performed.}
#' @seealso \link{PMD.cv}, \link{SPC}
#' @references
#' Witten D. M., Tibshirani R.,  and Hastie, T. (2009)
#' \emph{A penalized matrix decomposition, with applications to sparse principal components and canonical correlation analysis}, \emph{Biostatistics, Gol 10 (3), 515-534, Jul 2009}\cr
#' @examples
#'
#'
#' # Try PMD with L1 penalty on rows and columns: type="standard"
#' # A simple simulated example
#' set.seed(1)
#' # Our data is a rank-one matrix, plus noise. The underlying components
#' # contain 50 and 75 non-zero elements, respectively.
#' u <- matrix(c(rnorm(50), rep(0,150)),
#' ncol=1)
#' v <- matrix(c(rnorm(75),rep(0,225)), ncol=1)
#' x <- u%*%t(v)+
#' matrix(rnorm(200*300),ncol=300)
#' # We can use cross-validation to try to find optimal value of sumabs
#' cv.out <- PMD.cv(x, type="standard", sumabss=seq(0.1, 0.6, len=20))
#' print(cv.out)
#' plot(cv.out)
#' # The optimal value of sumabs is 0.4157, but we can get within one
#' # standard error of that CV error using sumabs=0.337, which corresponds to
#' # an average of 45.8 and 71.8 non-zero elements in each component - pretty
#' # close to the true model.
#' # We can fit the model corresponding to the lowest cross-validation error:
#' out <- PMD(x, type="standard", sumabs=cv.out$bestsumabs, K=1, v=cv.out$v.init)
#' print(out)
#' par(mfrow=c(2,2))
#' par(mar=c(2,2,2,2))
#' plot(out$u[,1], main="Est. u")
#' plot(out$v[,1], main="Est. v")
#' plot(u, main="True u")
#' plot(v, main="True v")
#' # And if we want to control sumabsu and sumabsv separately, we can do
#' # that too. Let's get 2 components while we're at it:
#' out2 <- PMD(x, type="standard",  K=2, sumabsu=6, sumabsv=8, v=out$v.init,
#' cnames=paste("v", sep=" ", 1:ncol(x)), rnames=paste("u", sep=" ", 1:nrow(x)))
#' print(out2)
#'
#' # Now check out PMD with L1 penalty on rows and fused lasso penalty on
#' # columns: type="ordered". We'll use the Chin et al (2006) Cancer Cell
#' # data set; try "?breastdata" for more info.
#' data(breastdata)
#' attach(breastdata)
#' # dna contains CGH data and chrom contains chromosome of each CGH spot;
#' # nuc contains position of each CGH spot.
#' dna <- t(dna) # Need samples on rows and CGH spots on columns
#' # First, look for shared regions of gain/loss on chromosome 1.
#' # Use cross-validation to choose tuning parameter value
#' par(mar=c(2,2,2,2))
#' cv.out <- PMD.cv(dna[,chrom==1],type="ordered",chrom=chrom[chrom==1],
#' nuc=nuc[chrom==1],
#' sumabsus=seq(1, sqrt(nrow(dna)), len=15))
#' print(cv.out)
#' plot(cv.out)
#' out <- PMD(dna[,chrom==1],type="ordered",
#' sumabsu=cv.out$bestsumabsu,chrom=chrom[chrom==1],K=1,v=cv.out$v.init,
#' cnames=paste("Pos",sep="",
#' nuc[chrom==1]), rnames=paste("Sample", sep=" ", 1:nrow(dna)))
#' print(out, verbose=TRUE)
#' # Which samples actually have that region of gain/loss?
#' par(mfrow=c(3,1))
#' par(mar=c(2,2,2,2))
#' PlotCGH(dna[which.min(out$u[,1]),chrom==1],chrom=chrom[chrom==1],
#' main=paste(paste(paste("Sample ", sep="", which.min(out$u[,1])),
#' sep="; u=", round(min(out$u[,1]),3))),nuc=nuc[chrom==1])
#' PlotCGH(dna[88,chrom==1], chrom=chrom[chrom==1],
#' main=paste("Sample 88; u=", sep="", round(out$u[88,1],3)),
#' nuc=nuc[chrom==1])
#' PlotCGH(out$v[,1],chrom=chrom[chrom==1], main="V",nuc=nuc[chrom==1])
#'
#'
#' detach(breastdata)
#'
#'
#'
#'
#'
#' @export PMD
PMD <- function(x, type=c("standard", "ordered"), sumabs=.4, sumabsu=5, sumabsv=NULL, lambda=NULL, niter=20, K=1, v=NULL, trace=TRUE, center=TRUE, chrom=NULL, rnames=NULL, cnames=NULL, upos=FALSE, uneg=FALSE, vpos=FALSE, vneg=FALSE){
  if(upos&&uneg || vpos&&vneg) stop("Cannot contrain elements to be both positive and negative!")
  if(type=="ordered" && (vpos||vneg)) stop("Cannot constrain signs of v if type is ordered.")
  call <- match.call()
  type <- match.arg(type)
  if(type=="standard"){
    out <- PMDL1L1(x, sumabs=sumabs, sumabsu=sumabsu, sumabsv=sumabsv, niter=niter, K=K, v=v, trace=trace, center=center, rnames=rnames, cnames=cnames, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg)
    class(out) <- "PMDL1L1"
  } else if(type=="ordered"){
    out <- PMDL1FL(x,K=K,sumabsu=sumabsu,lambda=lambda,chrom=chrom,niter=niter, v=v, trace=trace, center=center, rnames=rnames, cnames=cnames, upos=upos, uneg=uneg)
    class(out) <- "PMDL1FL"
  }
  out$upos <- upos
  out$uneg <- uneg
  out$vpos <- vpos
  out$vneg <- vneg
  out$call <- call
  return(out)
}





#' Do tuning parameter selection for PMD via cross-validation
#'
#' Performs cross-validation to select tuning parameters for rank-1 PMD, the
#' penalized matrix decomposition for a data matrix.
#'
#' If type is "standard", then lasso ($L_1$) penalties (promoting sparsity) are
#' placed on u and v. If type is "ordered", then lasso penalty is placed on u
#' and a fused lasso penalty (promoting sparsity and smoothness) is placed on
#' v.
#'
#' Cross-validation of the rank-1 PMD is performed over sumabss (if type is
#' "standard") or over sumabsus (if type is "ordered"). If type is "ordered",
#' then lambda is chosen from the data without cross-validation.
#'
#' The cross-validation works as follows: Some percent of the elements of $x$
#' is removed at random from the data matrix. The PMD is performed for a range
#' of tuning parameter values on this partially-missing data matrix; then,
#' missing values are imputed using the decomposition obtained. The value of
#' the tuning parameter that results in the lowest sum of squared errors of the
#' missing values if "best".
#'
#' To do cross-validation on the rank-2 PMD, first the rank-1 PMD should be
#' computed, and then this function should be performed on the residuals, given
#' by $x-udv'$.
#'
#' @param x Data matrix of dimension $n x p$, which can contain NA for missing
#' values.
#' @param type "standard" or "ordered": Do we want v to simply be sparse, or
#' should it also be smooth? If the columns of x are ordered (e.g. CGH spots
#' along a chromosome) then choose "ordered". Default is "standard". If
#' "standard", then the PMD function will make use of sumabs OR
#' sumabsu&sumabsv. If "ordered", then the function will make use of sumabsu
#' and lambda.
#' @param sumabss Used only if type is "standard". A vector of sumabs values to
#' be used. Sumabs is a measure of sparsity for u and v vectors, between 0 and
#' 1. When sumabss is specified, and sumabsus and sumabsvs are NULL, then
#' sumabsus is set to $sqrt(n)*sumabss$ and sumabsvs is set at
#' $sqrt(p)*sumabss$. If sumabss is specified, then sumabsus and sumabsvs
#' should be NULL. Or if sumabsus and sumabsvs are specified, then sumabss
#' should be NULL.
#' @param sumabsus Used only for type "ordered". A vector of sumabsu values to
#' be used. Sumabsu measures sparseness of u - it is the sum of absolute values
#' of elements of u. Must be between 1 and sqrt(n).
#' @param lambda Used only if type is "ordered". This is the tuning parameter
#' for the fused lasso penalty on v, which takes the form $lambda ||v||_1 +
#' lambda |v_j - v_(j-1)|$. $lambda$ must be non-negative. If NULL, then it is
#' chosen adaptively from the data.
#' @param nfolds How many cross-validation folds should be performed?  Default
#' is 5.
#' @param niter How many iterations should be performed. For speed, only 5 are
#' performed by default.
#' @param v The first right singular vector(s) of the data. (If missing data is
#' present, then the missing values are imputed before the singular vectors are
#' calculated.) v is used as the initial value for the iterative PMD algorithm.
#' If x is large, then this step can be time-consuming; therefore, if PMD is to
#' be run multiple times, then v should be computed once and saved.
#' @param chrom If type is "ordered", then this gives the option to specify
#' that some columns of x (corresponding to CGH spots) are on different
#' chromosomes. Then v will be sparse, and smooth *within* each chromosome but
#' not *between* chromosomes. Length of chrom should equal number of columns of
#' x, and each entry in chrom should be a number corresponding to which
#' chromosome the CGH spot is on.
#' @param nuc If type is "ordered", can specify the nucleotide position of each
#' CGH spot (column of x), to be used in plotting. If NULL, then it is assumed
#' that CGH spots are equally spaced.
#' @param trace Print out progress as iterations are performed? Default is
#' TRUE.
#' @param center Subtract out mean of x? Default is TRUE
#' @param upos Constrain the elements of u to be positive? TRUE or FALSE.
#' @param uneg Constrain the elements of u to be negative? TRUE or FALSE.
#' @param vpos Constrain the elements of v to be positive? TRUE or FALSE.
#' Cannot be used if type is "ordered".
#' @param vneg Constrain the elements of v to be negative? TRUE or FALSE.
#' Cannot be used if type is "ordered."
#' @return \item{cv}{Average sum of squared errors obtained over
#' cross-validation folds.} \item{cv.error}{Standard error of average sum of
#' squared errors obtained over cross-validation folds.} \item{bestsumabs}{If
#' type="standard", then value of sumabss resulting in smallest CV error is
#' returned.} \item{bestsumabsu}{If type="ordered", then value of sumabsus
#' resulting in smallest CV error is returned.} \item{v.init}{The first right
#' singular vector(s) of the data; these are returned to save on computation
#' time if PMD will be run again.}
#' @seealso \link{PMD}, \link{SPC}
#' @references
#' Witten D. M., Tibshirani R.,  and Hastie, T. (2009)
#' \emph{A penalized matrix decomposition, with applications to sparse principal components and canonical correlation analysis}, \emph{Biostatistics, Gol 10 (3), 515-534, Jul 2009}\cr
#' @examples
#' # See examples in PMD help file
#' @export PMD.cv
PMD.cv <- function(x, type=c("standard", "ordered"), sumabss=seq(0.1,0.7,len=10), sumabsus=NULL, lambda=NULL, nfolds=5, niter=5, v=NULL, chrom=NULL, nuc=NULL, trace=TRUE, center=TRUE, upos=FALSE, uneg=FALSE, vpos=FALSE, vneg=FALSE){
  if(upos&&uneg || vpos&&vneg) stop("Cannot contrain elements to be both positive and negative!")
  if(type=="ordered" && (vpos||vneg)) stop("Cannot constrain signs of v if type is ordered.")
  call <- match.call()
  type <- match.arg(type)
  if(type=="standard"){
    out <- PMDL1L1.cv(x, sumabss=sumabss, nfolds=nfolds, niter=niter, v=v, trace=trace, center=center, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg)
    class(out) <- "PMDL1L1.cv"
  } else if(type=="ordered"){
    out <- PMDL1FL.cv(x, nfolds=nfolds, niter=niter, lambda=lambda, sumabsus=sumabsus, chrom=chrom, v=v, trace=trace, nuc=nuc, center=center, upos=upos, uneg=uneg)
    class(out) <- "PMDL1FL.cv"
  }
  out$upos <- upos
  out$uneg <- uneg
  out$vpos <- vpos
  out$vneg <- vneg
  out$call <- call
  return(out)
}



PMDL1L1 <- function(x,sumabs=.4,sumabsu=NULL,sumabsv=NULL,niter=20,K=1,v=NULL, trace=TRUE, orth=FALSE, center=TRUE, rnames=NULL, cnames=NULL, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg){
  if(center){
    meanx <- mean_na(x)
    x <- x-meanx
  } else {
    meanx <- NULL
  }
  if(orth){
    if(is.null(sumabsu) || sumabsu < sqrt(nrow(x))){
      orth <- FALSE
      warning("Orth option ignored because sparse PCA results only when sumabsu equals sqrt(nrow(x))")
    }
  }
#  if((is.null(sumabsu) && !is.null(sumabsv)) || (is.null(sumabsv) && !is.null(sumabsu))) warning("Sumabsu and sumabsv BOTH must be input when type=standard. Since only one was given, it was ignored and sumabs was used instead.")
  if(is.null(sumabsu) || is.null(sumabsv)){
    sumabsu <- sqrt(nrow(x))*sumabs
    sumabsv <- sqrt(ncol(x))*sumabs
  }
  call <-  match.call()
  if(trace && abs(mean_na(x)) > 1e-15) warning("PMDL1L1 was run without first subtracting out the mean of x.")
  if(!is.null(sumabsu) && (sumabsu<1 || sumabsu>sqrt(nrow(x)))) stop("sumabsu must be between 1 and sqrt(n)")
  if(!is.null(sumabsv) && (sumabsv<1 || sumabsv>sqrt(ncol(x)))) stop("sumabsv must be between 1 and sqrt(p)")
  v <- CheckPMDV(v,x,K)
  if(K>1 && !orth) out <- (MultiSMD(x,sumabsu=sumabsu,sumabsv=sumabsv,niter=niter,K=K, trace=trace, v=v, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg))
  if(K>1 && orth) out <- MultiSMDOrth(x,sumabsu=sumabsu,sumabsv=sumabsv,niter=niter,K=K, trace=trace, v=v,  vpos=vpos, vneg=vneg)
  if(K==1) out <- SMD(x,sumabsu=sumabsu,sumabsv=sumabsv,niter=niter, trace=trace, v=v, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg)
  obj <- (list(u=out$u,v=out$v, d=out$d, v.init=out$v.init, call=call, meanx=meanx,sumabsu=sumabsu, sumabsv=sumabsv, rnames=rnames, cnames=cnames, K=K, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg))
  class(obj) <- "PMDL1L1"
  return(obj)
}

#' @method print PMDL1L1
#' @export
print.PMDL1L1 <- function(x,verbose=FALSE,...){
  cat("Call: ")
  dput(x$call)
  cat("\n\n")
  cat("Num non-zero u's: ", apply(x$u!=0, 2, sum), "\n")
  cat("Num non-zero v's: ", apply(x$v!=0, 2, sum), "\n")
  cat("Sumabsu used: ", x$sumabsu, "\n")
  cat("Sumabsv used: ", x$sumabsv, "\n")
  if(!is.null(x$meanx)) cat("Subtracted out mean of x: ", x$meanx, "\n")
  if(verbose){
    for(k in 1:x$K){
      cat("\n Component ", k, ":\n")
      u <- x$u[,k]
      v <- x$v[,k]
      if(is.null(x$rnames)) x$rnames <- 1:length(u)
      if(is.null(x$cnames)) x$cnames <- 1:length(v)
      cat(fill=T)
      us <- cbind(x$rnames[u!=0], round(u[u!=0],3))
      dimnames(us) <- list(1:sum(u!=0), c("Row Feature Name", "Row Feature Weight"))
      vs <- cbind(x$cnames[v!=0], round(v[v!=0],3))
      dimnames(vs) <- list(1:sum(v!=0), c("Column Feature Name", "Column Feature Weight"))
      print(us, quote=FALSE, sep="\t")
      cat(fill=T)
      print(vs, quote=FALSE, sep="\t")
    }
  }
  if(x$upos) cat("Elements of u constrained to be positive.", fill=TRUE)
  if(x$vpos) cat("Elements of v constrained to be positive.", fill=TRUE)
  if(x$uneg) cat("Elements of u constrained to be negative.", fill=TRUE)
  if(x$vneg) cat("Elements of v constrained to be negative.", fill=TRUE)
}

#' @method print PMDL1FL
#' @export
print.PMDL1FL <- function(x,verbose=FALSE,...){
  cat("Call: ")
  dput(x$call)
  cat("\n\n")
  cat("Num non-zero u's: ", apply(x$u!=0, 2, sum), "\n")
  cat("Num non-zero v's: ", apply(x$v!=0, 2, sum), "\n")
  cat("Sumabsu used: ", x$sumabsu, "\n")
  cat("Lambda used: ", x$lambda, "\n")
  if(!is.null(x$meanx)) cat("Subtracted out mean of x: ", x$meanx, "\n")
    if(verbose){
    for(k in 1:x$K){
      cat("\n Component ", k, ":\n")
      u <- x$u[,k]
      v <- x$v[,k]
      if(is.null(x$rnames)) x$rnames <- 1:length(u)
      if(is.null(x$cnames)) x$cnames <- 1:length(v)
      cat(fill=T)
      us <- cbind(x$rnames[u!=0], round(u[u!=0],3))
      dimnames(us) <- list(1:sum(u!=0), c("Row Feature Name", "Row Feature Weight"))
      vs <- cbind(x$cnames[v!=0], round(v[v!=0],3))
      dimnames(vs) <- list(1:sum(v!=0), c("Column Feature Name", "Column Feature Weight"))
      print(us, quote=FALSE, sep="\t")
      cat(fill=T)
      print(vs, quote=FALSE, sep="\t")
    }
  }
  if(x$upos) cat("Elements of u constrained to be positive.", fill=TRUE)
  if(x$uneg) cat("Elements of u constrained to be negative.", fill=TRUE)
}

MultiSMDOrth <- function(x,sumabsu,sumabsv,niter,K, trace, v, vpos, vneg){
  if(sumabsu < sqrt(nrow(x))){
    warning("sumabsu was less than sqrt(nrow(x)), so the orthogonal option was not implemented.")
    return(MultiSMD(x,sumabsu=sumabsu,sumabsv=sumabsv,niter=niter,K=K, trace=trace, v=v, vpos=vpos, vneg=vneg))
  }
  nas <- is.na(x)
  v.init <- v
  xuse <- x
  ds <- numeric(K)
  us <- matrix(0,nrow=nrow(x),ncol=K)
  vs <- matrix(0,nrow=ncol(x),ncol=K)
  for(k in 1:K){
    if(k==1) out <- SMD(xuse,sumabsu=sumabsu,sumabsv=sumabsv,niter=niter,v=matrix(v[,k],ncol=1), trace=trace, upos=FALSE, uneg=FALSE,  vpos=vpos, vneg=vneg)
    if(k>1) out <- SMDOrth(xuse,us[,1:(k-1)],sumabsv=sumabsv,niter=niter,v=matrix(v[,k],ncol=1), trace=trace,  vpos=vpos, vneg=vneg)
    us[,k] <- out$u
    vs[,k] <- out$v
    ds[k] <- out$d
  }
  return(list(u=us,v=vs,d=ds, v.init=v.init))
}

MultiSMD <- function(x, sumabsu, sumabsv, K=3, niter=20,v, trace=TRUE, upos, uneg, vpos, vneg){
  nas <- is.na(x)
  v.init <- v
  xuse <- x
  ds <- numeric(K)
  us <- matrix(0,nrow=nrow(x),ncol=K)
  vs <- matrix(0,nrow=ncol(x),ncol=K)
  for(k in 1:K){
    out <- SMD(xuse, sumabsu=sumabsu,sumabsv=sumabsv,niter=niter,v=matrix(v[,k],ncol=1), trace=trace, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg)
    us[,k] <- out$u
    vs[,k] <- out$v
    ds[k] <- out$d
    res <- xuse - out$d*out$u%*%t(out$v)
    xuse[!nas] <- res[!nas] # [!nas] is new on July 24 2009
  }
  return(list(u=us,v=vs,d=ds, v.init=v.init))
}

PMDL1L1.cv <- function(x, sumabss=seq(0.1,0.7,len=10), nfolds=5, niter=5, v=NULL, trace=TRUE, center=TRUE, upos, uneg, vpos, vneg){
  if(nfolds<2) stop("Must run at least 2 cross-validation folds.")
  percentRemove <- min(.25, 1/nfolds)
  call <- match.call()
  if(max(sumabss)>1 || min(sumabss)<0) stop("sumabs must be between 0 and 1")
  # PMD only cross-validates the single-factor model.
  # and it requires just one sparsity  measure.
  xfill <- x
  missing <- is.na(x)
  v <- CheckPMDV(v,x,K=1)
  errs <- matrix(NA, nrow=nfolds, ncol=length(sumabss))
  nonzerous <- matrix(NA, nrow=nfolds, ncol=length(sumabss))
  nonzerovs <- matrix(NA, nrow=nfolds, ncol=length(sumabss))
  rands <- matrix(runif(nrow(x)*ncol(x)),ncol=ncol(x))
  for(i in 1:nfolds){
    if(trace) cat(" Fold ", i, " out of ", nfolds, "\n")
    rm <- ((i-1)*percentRemove < rands)&(rands < i*percentRemove)
    xrm <- x
    xrm[rm] <- NA
    for(j in 1:length(sumabss)){
      out <- PMDL1L1(xrm, sumabs=sumabss[j], niter=niter, v=v, trace=FALSE, center=center,K=1, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg)
      xhat <- as.numeric(out$d)*out$u%*%t(out$v)
      errs[i,j] <- sum(((xhat-x)[rm & !missing])^2)
      nonzerous[i,j] <- sum(out$u!=0)
      nonzerovs[i,j] <- sum(out$v!=0)
    }
  }
  if(trace) cat(fill=TRUE)
  err.means <- apply(errs, 2, mean)
  err.sds <- apply(errs, 2, sd)/sqrt(nfolds)
  nonzerous.mean <- apply(nonzerous,2,mean)
  nonzerovs.mean <- apply(nonzerovs,2,mean)
  bestsumabs <- sumabss[which.min(err.means)]
  object <- (list(cv=err.means, cv.error=err.sds,bestsumabs=bestsumabs, nonzerous=nonzerous.mean, nonzerovs=nonzerovs.mean, v.init=v, call=call, sumabss=sumabss, nfolds=nfolds))
  class(object) <- "PMDL1L1.cv"
  return(object)
}





#' Perform cross-validation on sparse principal component analysis
#'
#' Selects tuning parameter for the sparse principal component analysis method
#' of Witten, Tibshirani, and Hastie (2008), which involves applying PMD to a
#' data matrix with lasso ($L_1$) penalty on the columns and no penalty on the
#' rows. The tuning parameter controls the sum of absolute values - or $L_1$
#' norm - of the elements of the sparse principal component.
#'
#' This method only performs cross-validation for the first sparse principal
#' component. It does so by performing the following steps nfolds times: (1)
#' replace a fraction of the data with missing values, (2) perform SPC on this
#' new data matrix using a range of tuning parameter values, each time getting
#' a rank-1 approximationg $udv'$ where $v$ is sparse, (3) measure the mean
#' squared error of the rank-1 estimate of the missing values created in step
#' 1.
#'
#' Then, the selected tuning parameter value is that which resulted in the
#' lowest average mean squared error in step 3.
#'
#' In order to perform cross-validation for the second sparse principal
#' component, apply this function to $X-udv'$ where $udv'$ are the output of
#' running SPC on the raw data $X$.
#'
#' @aliases SPC.cv plot.SPC.cv print.SPC.cv
#' @param x Data matrix of dimension $n x p$, which can contain NA for missing
#' values. We are interested in finding sparse principal components of
#' dimension $p$.
#' @param sumabsvs Range of sumabsv values to be considered in
#' cross-validation. Sumabsv is the sum of absolute values of elements of v. It
#' must be between 1 and square root of number of columns of data. The smaller
#' it is, the sparser v will be.
#' @param nfolds Number of cross-validation folds performed.
#' @param niter How many iterations should be performed. By default, perform
#' only 5 for speed reasons.
#' @param v The first right singular vector(s) of the data. (If missing data is
#' present, then the missing values are imputed before the singular vectors are
#' calculated.) v is used as the initial value for the iterative PMD($L_1$,
#' $L_1$) algorithm. If x is large, then this step can be time-consuming;
#' therefore, if PMD is to be run multiple times, then v should be computed
#' once and saved.
#' @param trace Print out progress as iterations are performed? Default is
#' TRUE.
#' @param orth If TRUE, then use method of Section 3.2 of Witten, Tibshirani
#' and Hastie (2008) to obtain multiple sparse principal components. Default is
#' FALSE.
#' @param center Subtract out mean of x? Default is TRUE
#' @param vpos Constrain elements of v to be positive? Default is FALSE.
#' @param vneg Constrain elements of v to be negative? Default is FALSE.
#' @return \item{cv}{Average sum of squared errors that results for each tuning
#' parameter value.} \item{cv.error}{Standard error of the average sum of
#' squared error that results for each tuning parameter value.}
#' \item{bestsumabsv}{Value of sumabsv that resulted in lowest CV error.}
#' \item{nonzerovs}{Average number of non-zero elements of v for each candidate
#' value of sumabsvs.} \item{v.init}{Initial value of v that was passed in. Or,
#' if that was NULL, then first right singular vector of X.}
#' \item{bestsumabsv1se}{The smallest value of sumabsv that is within 1
#' standard error of smallest CV error.}
#' @seealso \link{SPC}, \link{PMD}, \link{PMD.cv}
#' @references
#' Witten D. M., Tibshirani R.,  and Hastie, T. (2009)
#' \emph{A penalized matrix decomposition, with applications to sparse principal components and canonical correlation analysis}, \emph{Biostatistics, Gol 10 (3), 515-534, Jul 2009}\cr
#' @examples
#'
#' #NOT RUN
#' ## A simple simulated example
#' #set.seed(1)
#' #u <- matrix(c(rnorm(50), rep(0,150)),ncol=1)
#' #v <- matrix(c(rnorm(75),rep(0,225)), ncol=1)
#' #x <- u%*%t(v)+matrix(rnorm(200*300),ncol=300)
#' ## Perform Sparse PCA - that is, decompose a matrix w/o penalty on rows
#' ## and w/ L1 penalty on columns
#' ## First, we perform sparse PCA and get 4 components, but we do not
#' ## require subsequent components to be orthogonal to previous components
#' #cv.out <- SPC.cv(x, sumabsvs=seq(1.2, sqrt(ncol(x)), len=6))
#' #print(cv.out)
#' #plot(cv.out)
#' #out <- SPC(x,sumabsv=cv.out$bestsumabs, K=4) # could use
#' ## cv.out$bestsumabvsv1se instead
#' #print(out,verbose=TRUE)
#' ## Now, we do sparse PCA using method in Section 3.2 of WT&H(2008) for getting
#' ## multiple components - that is, we require components to be orthogonal
#' #cv.out <- SPC.cv(x, sumabsvs=seq(1.2, sqrt(ncol(x)), len=6), orth=TRUE)
#' #print(cv.out)
#' #plot(cv.out)
#' #out.orth <- SPC(x,sumabsv=cv.out$bestsumabsv, K=4, orth=TRUE)
#' #print(out.orth,verbose=TRUE)
#' #par(mfrow=c(1,1))
#' #plot(out$u[,1], out.orth$u[,1], xlab="", ylab="")
#' #
#' #
#' @export SPC.cv
SPC.cv <- function(x, sumabsvs=seq(1.2,5,len=10), nfolds=5, niter=5, v=NULL, trace=TRUE, orth=FALSE,  center=TRUE, vpos=FALSE, vneg=FALSE){
  if(vpos&&vneg) stop("Cannot constrain elements of v to be both positive and negative.")
  if(nfolds<2) stop("Must run at least 2 cross-validation folds.")
  percentRemove <- min(.25, 1/nfolds)
  call <- match.call()
  if(max(sumabsvs)>sqrt(ncol(x)) || min(sumabsvs)<1) stop("sumabs must be between 1 and sqrt(ncol(x))")
  xfill <- x
  missing <- is.na(x)
  v <- CheckPMDV(v,x,K=1)
  errs <- matrix(NA, nrow=nfolds, ncol=length(sumabsvs))
  nonzerovs <- matrix(NA, nrow=nfolds, ncol=length(sumabsvs))
  rands <- matrix(runif(nrow(x)*ncol(x)),ncol=ncol(x))
  for(i in 1:nfolds){
    if(trace) cat(" Fold ", i, " out of ", nfolds, "\n")
    rm <- ((i-1)*percentRemove < rands)&(rands < i*percentRemove)
    xrm <- x
    xrm[rm] <- NA
    for(j in 1:length(sumabsvs)){
      out <- SPC(xrm, sumabsv=sumabsvs[j], orth=orth, niter=niter, v=v, trace=FALSE, center=center,K=1, vpos=vpos, vneg=vneg)
      xhat <- as.numeric(out$d)*out$u%*%t(out$v)
      errs[i,j] <- sum(((xhat-x)[rm & !missing])^2)
      nonzerovs[i,j] <- sum(out$v!=0)
    }
  }
  if(trace) cat(fill=TRUE)
  err.means <- apply(errs, 2, mean)
  err.sds <- apply(errs, 2, sd)/sqrt(nfolds)
  nonzerovs.mean <- apply(nonzerovs,2,mean)
  bestsumabsv <- sumabsvs[which.min(err.means)]
  bestsumabsv1se <- sumabsvs[min(which(err.means < min(err.means) + err.sds[which.min(err.means)]))]
  object <- (list(cv=err.means, cv.error=err.sds,bestsumabsv=bestsumabsv,nonzerovs=nonzerovs.mean, v.init=v, call=call, sumabsvs=sumabsvs, nfolds=nfolds, bestsumabsv1se=bestsumabsv1se, vpos=vpos, vneg=vneg))
  class(object) <- "SPC.cv"
  return(object)
}

#' @method plot SPC.cv
#' @export
plot.SPC.cv <- function(x,...){
  sumabsvs <- x$sumabsvs
  err.means <- x$cv
  err.sds <- x$cv.error
  plot(sumabsvs, err.means, ylim=range(c(err.means+err.sds, err.means-err.sds)), main="CV Error for SPC", xlab="Value of sumabsv", ylab="Average SSE", type="l")
  points(sumabsvs, err.means)
  lines(sumabsvs, err.means+err.sds, lty="dashed")
  lines(sumabsvs, err.means-err.sds, lty="dashed")
}

#' @method print SPC.cv
#' @export
print.SPC.cv <- function(x,...){
  cat("Call:\n")
  dput(x$call)
  cat("\n \n")
  cat('Cross-validation errors: \n')
  mat <- round(cbind(round(x$sumabsvs, 3), x$cv, x$cv.error, x$nonzerovs),3)
  dimnames(mat) <- list(1:length(x$sumabsvs), c("Sumabsvs", "CV Error", "CV S.E.", "# non-zero v's"))
  print(mat, quote=F)
  cat('\n Best sumabsv value (lowest CV error): ', x$bestsumabsv, "\n")
  cat('\n Smallest sumabsv value that has CV error within 1 SE of best CV error: ', x$bestsumabsv1se, "\n")
  if(x$vpos) cat("Elements of v constrained to be positive.", fill=TRUE)
  if(x$vneg) cat("Elements of v constrained to be negative.", fill=TRUE)
}


#' @method print PMDL1L1.cv
#' @export
print.PMDL1L1.cv <- function(x,...){
  cat("Call:\n")
  dput(x$call)
  cat("\n \n")
  cat('Cross-validation errors: \n')
  mat <- round(cbind(round(x$sumabss, 3), x$cv, x$cv.error, x$nonzerous, x$nonzerovs),3)
  dimnames(mat) <- list(1:length(x$sumabss), c("Sumabss", "CV Error", "CV S.E.", "# non-zero u's", "# non-zero v's"))
  print(mat, quote=F)
  cat('\n Best sumabs value (lowest CV error): ', x$bestsumabs, "\n")
  if(x$upos) cat("Elements of u constrained to be positive.", fill=TRUE)
  if(x$uneg) cat("Elements of u constrained to be negative.", fill=TRUE)
  if(x$vpos) cat("Elements of v constrained to be positive.", fill=TRUE)
  if(x$vneg) cat("Elements of v constrained to be negative.", fill=TRUE)
}

#' @method plot PMDL1L1.cv
#' @export
plot.PMDL1L1.cv <- function(x,...){
  sumabss <- x$sumabss
  err.means <- x$cv
  err.sds <- x$cv.error
  plot(sumabss, err.means, ylim=range(c(err.means+err.sds, err.means-err.sds)), main="CV Error for PMD(L1,L1)", xlab="Value of sumabs", ylab="Average SSE", type="l")
  points(sumabss, err.means)
  lines(sumabss, err.means+err.sds, lty="dashed")
  lines(sumabss, err.means-err.sds, lty="dashed")
}


PMDL1FL.cv <- function(x, nfolds=5, niter=5, lambda=NULL, sumabsus=NULL, chrom=NULL, v=NULL, trace=TRUE, nuc=NULL, center=TRUE, upos, uneg){
  if(is.null(sumabsus)) sumabsus <- seq(1, .7*sqrt(nrow(x)),len=10)
  percentRemove <- min(.25, 1/nfolds)
  call <- match.call()
  if(is.null(chrom)) chrom <- rep(1, ncol(x))
  v <- CheckPMDV(v,x,K=1)
  if(is.null(lambda)) lambda <- ChooseLambda1Lambda2(as.numeric(v[,1]))
  errs <- matrix(NA, nrow=nfolds, ncol=length(sumabsus))
  nonzerous <- matrix(NA, nrow=nfolds, ncol=length(sumabsus))
  nonzerovs <- matrix(NA, nrow=nfolds, ncol=length(sumabsus))
  rands <- matrix(runif(nrow(x)*ncol(x)), nrow=nrow(x))
  missing <- is.na(x)
  for(i in 1:nfolds){
    if(trace) cat(" Fold ", i, " out of ", nfolds, "\n")
    rm <- ((i-1)*percentRemove < rands)&(rands < i*percentRemove)
    xrm <- x
    xrm[rm] <- NA
    for(j in 1:length(sumabsus)){
      if(trace) cat(j,fill=F)
      out <- PMDL1FL(xrm, K=1, lambda=lambda,sumabsu=sumabsus[j], niter=niter, chrom=chrom, v=v, trace=FALSE, center=center, upos=upos, uneg=uneg)
      xhat <- as.numeric(out$d)*out$u%*%t(out$v)
      errs[i,j] <- sum(((xhat-x)[rm & !missing])^2)
      nonzerous[i,j] <- sum(out$u!=0)
      nonzerovs[i,j] <- sum(out$v!=0)
    }
  }
  if(trace) cat(fill=TRUE)
  err.means <- apply(errs, 2, mean)
  err.sds <- apply(errs, 2, sd)/sqrt(nfolds)
  nonzerous.mean <- apply(nonzerous,2,mean)
  nonzerovs.mean <- apply(nonzerovs,2,mean)
  bestind <- which.min(err.means)
  bestsumabsu <- sumabsus[bestind]
  object <- (list(cv=err.means, lambda=lambda,cv.error=err.sds,bestsumabsu=bestsumabsu, nonzerous=nonzerous.mean, nonzerovs=nonzerovs.mean, v.init=v, call=call, sumabsus=sumabsus, nfolds=nfolds,x=x,chrom=chrom,nuc=nuc, niter=niter))
  class(object) <- "PMDL1FL.cv"
  return(object)
}

#' @method plot PMDL1FL.cv
#' @export
plot.PMDL1FL.cv <- function(x,...){
  mat <- x$x
  sumabsus <- x$sumabsus
  err.means <- x$cv
  err.sds <- x$cv.error
  chrom <- x$chrom
  niter <- x$niter
  nuc <- x$nuc
  lambda <- x$lambda
  v <- x$v.init
  bestsumabsu <- x$bestsumabsu
  par(mfrow=c(2,1))
  plot(sumabsus, err.means, ylim=range(c(err.means+err.sds, err.means-err.sds)), main="CV Error for PMD(L1, FL)", xlab="Value of sumabsu", ylab="Average SSE", type="l")
  points(sumabsus, err.means)
  lines(sumabsus, err.means+err.sds, lty="dashed")
  lines(sumabsus, err.means-err.sds, lty="dashed")
  out <- PMDL1FL(mat, sumabsu=bestsumabsu,  K=1, chrom=chrom, niter=niter, v=v, trace=FALSE, lambda=lambda, upos=x$upos, uneg=x$uneg)
  PlotCGH(out$v, main="Best V", chrom=chrom, nuc=nuc)
}

#' @method print PMDL1FL.cv
#' @export
print.PMDL1FL.cv <- function(x,...){
  cat("Call:\n")
  dput(x$call)
  cat("\n \n")
  cat('Cross-validation errors: \n')
  mat <- round(cbind(round(x$sumabsus, 3), x$cv, x$cv.error, x$nonzerous, x$nonzerovs),3)
  dimnames(mat) <- list(1:length(x$sumabsus), c("Sumabsus", "CV Error", "CV S.E.", "# non-zero u's", "# non-zero v's"))
  print(mat, quote=F)
  cat('\n Best sumabs value : ', x$bestsumabsu, "\n")
  cat('\n Lambda Used : ', x$lambda, "\n")
  if(x$upos) cat("Elements of u constrained to be positive.", fill=TRUE)
  if(x$uneg) cat("Elements of u constrained to be negative.", fill=TRUE)
}
