# This contains what used to be in CGH.SparseCCA.R and MultiSparseCCA.R

#' Perform sparse canonical correlation analysis using the penalized matrix
#' decomposition.
#'
#' Given matrices X and Z, which represent two sets of features on the same set
#' of samples, find sparse u and v such that u'X'Zv is large.  For X and Z, the
#' samples are on the rows and the features are on the columns. X and Z must
#' have same number of rows, but may (and usually will) have different numbers
#' of columns. The columns of X and/or Z can be unordered or ordered. If
#' unordered, then a lasso penalty will be used to obtain the corresponding
#' canonical vector. If ordered, then a fused lasso penalty will be used; this
#' will result in smoothness.
#'
#' This function is useful for performing an integrative analysis of two sets
#' of measurements taken on the same set of samples: for instance, gene
#' expression and CGH measurements on the same set of patients. It takes in two
#' data sets, called x and z, each of which have (the same set of) samples on
#' the rows. If z is a matrix of CGH data with *ordered* CGH spots on the
#' columns, then use typez="ordered". If z consists of unordered columns, then
#' use typez="standard". Similarly for typex.
#'
#' This function performs the penalized matrix decomposition on the data matrix
#' $X'Z$. Therefore, the results should be the same as running the PMD function
#' on t(x)\%*\% z. However, when ncol(x)>>nrow(x) and ncol(z)>>nrow(z) then
#' using the CCA function is much faster because it avoids computation of
#' $X'Z$.
#'
#' The CCA criterion is as follows: find unit vectors $u$ and $v$ such that
#' $u'X'Zv$ is maximized subject to constraints on $u$ and $v$.  If
#' typex="standard" and typez="standard" then the constraints on $u$ and $v$
#' are lasso ($L_1$). If typex="ordered" then the constraint on $u$ is a fused
#' lasso penalty (promoting sparsity and smoothness). Similarly if
#' typez="ordered".
#'
#' When type x is "standard": the L1 bound of u is penaltyx*sqrt(ncol(x)).
#'
#' When typex is "ordered": penaltyx controls the amount of sparsity and
#' smoothness in u, via the fused lasso penalty: $lambda sum_j |u_j| + lambda
#' sum_j |u_j - u_(j-1)|$. If NULL, then it will be chosen adaptively from the
#' data.
#'
#' @aliases CCA print.CCA
#' @param x Data matrix; samples are rows and columns are features. Cannot
#' contain missing values.
#' @param z Data matrix; samples are rows and columns are features.  Cannot
#' contain missing values.
#' @param typex Are the columns of x unordered (type="standard") or ordered
#' (type="ordered")? If "standard", then a lasso penalty is applied to u, to
#' enforce sparsity. If "ordered" (generally used for CGH data), then a fused
#' lasso penalty is applied, to enforce both sparsity and smoothness.
#' @param typez Are the columns of z unordered (type="standard") or ordered
#' (type="ordered")? If "standard", then a lasso penalty is applied to v, to
#' enforce sparsity. If "ordered" (generally used for CGH data), then a fused
#' lasso penalty is applied, to enforce both sparsity and smoothness.
#' @param penaltyx The penalty to be applied to the matrix x, i.e. the penalty
#' that results in the canonical vector u. If typex is "standard" then the L1
#' bound on u is penaltyx*sqrt(ncol(x)). In this case penaltyx must be between
#' 0 and 1 (larger L1 bound corresponds to less penalization). If "ordered"
#' then it's the fused lasso penalty lambda, which must be non-negative (larger
#' lambda corresponds to more penalization).
#' @param penaltyz The penalty to be applied to the matrix z, i.e. the penalty
#' that results in the canonical vector v. If typez is "standard" then the L1
#' bound on v is penaltyz*sqrt(ncol(z)). In this case penaltyz must be between
#' 0 and 1 (larger L1 bound corresponds to less penalization). If "ordered"
#' then it's the fused lasso penalty lambda, which must be non-negative (larger
#' lambda corresponds to more penalization).
#' @param K The number of u's and v's desired; that is, the number of canonical
#' vectors to be obtained.
#' @param niter How many iterations should be performed? Default is 15.
#' @param v The first K columns of the v matrix of the SVD of X'Z. If NULL,
#' then the SVD of X'Z will be computed inside the CCA function. However, if
#' you plan to run this function multiple times, then save a copy of this
#' argument so that it does not need to be re-computed (since that process can
#' be time-consuming if X and Z both have high dimension).
#' @param trace Print out progress?
#' @param standardize Should the columns of x and z be centered (to have mean
#' zero) and scaled (to have standard deviation 1)? Default is TRUE.
#' @param xnames An optional vector of column names for x, defaults to `colnames(x)`
#' @param znames An optional vector of column names for z, defaults to `colnames(z)`
#' @param chromx Used only if typex is "ordered"; allows user to specify a
#' vector of length ncol(x) giving the chromosomal location of each CGH spot.
#' This is so that smoothness will be enforced within each chromosome, but not
#' between chromosomes.
#' @param chromz Used only if typez is "ordered"; allows user to specify a
#' vector of length ncol(z) giving the chromosomal location of each CGH spot.
#' This is so that smoothness will be enforced within each chromosome, but not
#' between chromosomes.
#' @param upos If TRUE, then require elements of u to be positive. FALSE by
#' default. Can only be used if type is "standard".
#' @param uneg If TRUE, then require elements of u to be negative. FALSE by
#' default.  Can only be used if type is "standard".
#' @param vpos If TRUE, require elements of v to be positive. FALSE by default.
#' Can only be used if type is "standard".
#' @param vneg If TRUE, require elements of v to be negative. FALSE by default.
#' Can only be used if type is "standard".
#' @param outcome If you would like to incorporate a phenotype into CCA
#' analysis - that is, you wish to find features that are correlated across the
#' two data sets and also correlated with a phenotype - then use one of
#' "survival", "multiclass", or "quantitative" to indicate outcome type.
#' Default is NULL.
#' @param y If outcome is not NULL, then this is a vector of phenotypes - one
#' for each row of x and z. If outcome is "survival" then these are survival
#' times; must be non-negative. If outcome is "multiclass" then these are class
#' labels (1,2,3,...). Default NULL.
#' @param cens If outcome is "survival" then these are censoring statuses for
#' each observation. 1 is complete, 0 is censored. Default NULL.
#' @return \item{u}{u is output. If you asked for multiple factors then each
#' column of u is a factor. u has dimension nxK if you asked for K factors.}
#' \item{v}{v is output. If you asked for multiple factors then each column of
#' v is a factor. v has dimension pxK if you asked for K factors.} \item{d}{A
#' vector of length K, which can alternatively be computed as the diagonal of
#' the matrix $u'X'Zv$.} \item{v.init}{The first K factors of the v matrix of
#' the SVD of x'z. This is saved in case this function will be re-run later.}
#' @seealso \link{PMD},\link{CCA.permute}
#' @references
#' Witten D. M., Tibshirani R.,  and Hastie, T. (2009)
#' \emph{A penalized matrix decomposition, with applications to sparse principal components and canonical correlation analysis}, \emph{Biostatistics, Gol 10 (3), 515-534, Jul 2009}\cr
#' @examples
#'
#' # first, do CCA with type="standard"
#' # A simple simulated example
#' set.seed(3189)
#' u <- matrix(c(rep(1,25),rep(0,75)),ncol=1)
#' v1 <- matrix(c(rep(1,50),rep(0,450)),ncol=1)
#' v2 <- matrix(c(rep(0,50),rep(1,50),rep(0,900)),ncol=1)
#' x <- u%*%t(v1) + matrix(rnorm(100*500),ncol=500)
#' z <- u%*%t(v2) + matrix(rnorm(100*1000),ncol=1000)
#' # Can run CCA with default settings, and can get e.g. 3 components
#' out <- CCA(x,z,typex="standard",typez="standard",K=3)
#' print(out,verbose=TRUE) # To get less output, just print(out)
#' # Or can use CCA.permute to choose optimal parameter values
#' perm.out <- CCA.permute(x,z,typex="standard",typez="standard",nperms=7)
#' print(perm.out)
#' plot(perm.out)
#' out <- CCA(x,z,typex="standard",typez="standard",K=1,
#' 	   penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz,
#' 	   v=perm.out$v.init)
#' print(out)
#'
#'
#' ##### The remaining examples are commented out, but uncomment to run: ######
#'
#' # Not run, to save time:
#' \dontrun{
#' ## Now try CCA with a constraint that elements of u must be negative and
#' ## elements of v must be positive:
#' perm.out <- CCA.permute(x,z,typex="standard",typez="standard",nperms=7,
#' penaltyxs=seq(.1,.7,len=10), penaltyzs=seq(.1,.7,len=10), uneg=TRUE, vpos=TRUE)
#' print(perm.out)
#' plot(perm.out)
#' out <- CCA(x,z,typex="standard",typez="standard",K=1,
#' 	   penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz,
#'            v=perm.out$v.init, uneg=TRUE, vpos=TRUE)
#' print(out)
#'
#'
#' ## Suppose we also have a quantitative outcome, y, and we want to find
#' ## features in x and z that are correlated with each other and with the
#' ## outcome:
#' y <- rnorm(nrow(x))
#' perm.out <- CCA.permute(x,z,typex="standard",typez="standard",
#' 			outcome="quantitative",y=y, nperms=6)
#' print(perm.out)
#' out<-CCA(x,z,typex="standard",typez="standard",outcome="quantitative",
#' 	 y=y,penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz)
#' print(out)
#'
#' ## now, do CCA with type="ordered"
#' ## Example involving the breast cancer data: gene expression + CGH
#' set.seed(22)
#' data(breastdata)
#' attach(breastdata)
#' dna <- t(dna)
#' rna <- t(rna)
#' perm.out <- CCA.permute(x=rna,z=dna[,chrom==1],typex="standard",
#' 		       	typez="ordered",nperms=5,penaltyxs=seq(.02,.7,len=10))
#' ## We run CCA using all gene exp. data, but CGH data on chrom 1 only.
#' print(perm.out)
#' plot(perm.out)
#' out <- CCA(x=rna,z=dna[,chrom==1], typex="standard", typez="ordered",
#' 	   penaltyx=perm.out$bestpenaltyx,
#'            v=perm.out$v.init, penaltyz=perm.out$bestpenaltyz,
#'            xnames=substr(genedesc,1,20),
#'            znames=paste("Pos", sep="", nuc[chrom==1]))
#' # Save time by inputting  lambda and v
#' print(out) # could do print(out,verbose=TRUE)
#' print(genechr[out$u!=0]) # Cool! The genes associated w/ gain or loss
#' ## on chrom 1 are located on chrom 1!!
#' par(mfrow=c(1,1))
#' PlotCGH(out$v, nuc=nuc[chrom==1], chrom=chrom[chrom==1],
#' main="Regions of gain/loss on Chrom 1 assoc'd with gene expression")
#' detach(breastdata)
#' }
#'
#' @export CCA
CCA <- function(x, z, typex=c("standard", "ordered"), typez=c("standard","ordered"), penaltyx=NULL, penaltyz=NULL, K=1, niter=15, v=NULL, trace=TRUE, standardize=TRUE, xnames=colnames(x), znames=colnames(z), chromx=NULL, chromz=NULL, upos=FALSE, uneg=FALSE, vpos=FALSE, vneg=FALSE, outcome=NULL, y=NULL, cens=NULL){
  if(ncol(x)<2) stop("Need at least two features in dataset x.")
  if(ncol(z)<2) stop("Need at least two features in dataset z.")
  if(upos && uneg) stop("At most one of upos and uneg should be TRUE!")
  if(vpos && vneg)  stop("At most one of vpos and vneg should be TRUE!")
  if(typez=="ordered" && (vpos||vneg)) stop("Cannot require elements of v to be positive or negative if typez is ordered")
  if(typex=="ordered" && (upos||uneg)) stop("Cannot require elements of u to be positive or negative if typex is ordered")
  typex <- match.arg(typex)
  typez <- match.arg(typez)
  call <- match.call()
  if(sum(is.na(x))+sum(is.na(z)) > 0) stop("Cannot have NAs in x or z")
  if(nrow(x)!=nrow(z)) stop("x and z must have same number of rows")
  if(standardize){
    sdx <- apply(x,2,sd)
    sdz <- apply(z,2,sd)
    if(min(sdx)==0) stop("Cannot standardize because some of the columns of x have std. dev. 0")
    if(min(sdz)==0) stop("Cannot standardize because some of the columns of z have std. dev. 0")
    x <- scale(x,TRUE,sdx)
    z <- scale(z,TRUE,sdz)
  }
  if(!is.null(outcome)){
    pheno.out <- CCAPhenotypeZeroSome(x,z,y,qt=.8, cens=cens, outcome=outcome, typex=typex, typez=typez)
    x <- pheno.out$x
    z <- pheno.out$z
  }
  if(typex=="standard" && !is.null(chromx)) warning("Chromx has no effect for type standard")
  if(typez=="standard" && !is.null(chromz)) warning("Chromz has no effect for type standard")
  v <- CheckVs(v,x,z,K)
  if(is.null(penaltyx)){
    if(typex=="standard") penaltyx <- .3#pmax(1.001,.3*sqrt(ncol(x)))/sqrt(ncol(x))
    if(typex=="ordered")  penaltyx <- ChooseLambda1Lambda2(as.numeric(CheckVs(NULL,z,x,1))) # v[,1] used to be NULL
  }
  if(is.null(penaltyz)){
    if(typez=="standard") penaltyz <- .3#pmax(1.001,.3*sqrt(ncol(z)))/sqrt(ncol(z))
    if(typez=="ordered")  penaltyz <-  ChooseLambda1Lambda2(as.numeric(CheckVs(NULL,x,z,1))) # ChooseLambda1Lambda2(as.numeric(v[,1]))
  }
  if(!is.null(penaltyx)){
    if(typex=="standard" && (penaltyx<0 || penaltyx>1)) stop("Penaltyx must be between 0 and 1 when typex is standard.")
    if(typex=="ordered" && penaltyx<0) stop("Penaltyx must be non-negative when typex is standard.")
  }
  if(!is.null(penaltyz)){
    if(typez=="standard" && (penaltyz<0 || penaltyz>1)) stop("Penaltyz must be between 0 and 1 when typez is standard.")
    if(typez=="ordered" && penaltyz<0) stop("Penaltyz must be non-negative when typez is standard.")
  }
  out <- CCAAlgorithm(x=x,z=z,v=v,typex=typex,typez=typez,penaltyx=penaltyx,penaltyz=penaltyz,K=K,niter=niter,trace=trace,chromx=chromx,chromz=chromz,upos=upos,uneg=uneg,vpos=vpos,vneg=vneg)
  out$outcome <- outcome
  out$call <- call
  out$xnames <- xnames
  out$znames <- znames
  out$typex<-typex
  out$typez<-typez
  out$penaltyx<-penaltyx
  out$penaltyz<-penaltyz
  out$K <- K
  out$niter <- niter
  out$upos <- upos
  out$uneg <- uneg
  out$vpos <- vpos
  out$vneg <- vneg
  out$xnames <- xnames
  out$znames <- znames
  out$v.init <- v
  out$cors <- numeric(K)
  for(k in 1:K){
    if(sum(out$u[,k]!=0)>0 && sum(out$v[,k]!=0)>0) out$cors[k] <- cor(x%*%out$u[,k],z%*%out$v[,k])
  }
  class(out) <- "CCA"
  return(out)
}

#' @method print CCA
#' @export
print.CCA <- function(x,verbose=FALSE,...){
  cat("Call: ")
  dput(x$call)
  cat("\n\n")
  cat("Num non-zeros u's: ", apply(x$u!=0,2,sum), "\n")
  cat("Num non-zeros v's: ", apply(x$v!=0,2,sum), "\n")
  cat("Type of x: ", x$typex,"\n")
  cat("Type of z: ", x$typez,"\n")
  if(x$typex=="standard") cat("Penalty for x: L1 bound is ", x$penaltyx, "\n")
  if(x$typez=="standard") cat("Penalty for z: L1 bound is ", x$penaltyz, "\n")
  if(x$typex=="ordered") cat("Penalty for x: Lambda is ", x$penaltyx, "\n")
  if(x$typez=="ordered") cat("Penalty for z: Lambda is ", x$penaltyz, "\n")
  if(x$upos) cat("U's constrained to be positive", fill=TRUE)
  if(x$uneg) cat("U's constrained to be negative", fill=TRUE)
  if(x$vpos) cat("V's constrained to be positive", fill=TRUE)
  if(x$vneg) cat("V's constrained to be negative", fill=TRUE)
  if(!is.null(x$outcome)) cat("Outcome used: ", x$outcome, fill=TRUE)
  cat("Cor(Xu,Zv): ", x$cors, fill=TRUE)
  if(verbose){
    for(k in 1:x$K){
      cat("\n Component ", k, ":\n")
      u <- x$u[,k]
      v <- x$v[,k]
      if(is.null(x$xnames)) x$xnames <- 1:length(u)
      if(is.null(x$znames)) x$znames <- 1:length(v)
      cat(fill=T)
      us <- cbind(x$xnames[u!=0], round(u[u!=0],3))
      dimnames(us) <- list(1:sum(u!=0), c("Row Feature Name", "Row Feature Weight"))
      vs <- cbind(x$znames[v!=0], round(v[v!=0],3))
      dimnames(vs) <- list(1:sum(v!=0), c("Column Feature Name", "Column Feature Weight"))
      print(us, quote=FALSE, sep="\t")
      cat(fill=T)
      print(vs, quote=FALSE, sep="\t")
    }
  }
}


CCAAlgorithm <- function(x,z,v,typex,typez,penaltyx,penaltyz,K,niter,trace,chromx,chromz,upos,uneg,vpos,vneg){
  if(typez!="ordered"){
    if(K>1) v.init <- v[apply(z^2,2,sum)!=0,]
    if(K==1) v.init <- v[apply(z^2,2,sum)!=0]
  } else {
    v.init <- v
  }
  v.init <- matrix(v.init,ncol=K)
  u=v=d=NULL
  xres <- x; zres <- z
  if(typex!="ordered") xres <- x[,apply(x^2,2,sum)!=0]
  if(typez!="ordered") zres <- z[,apply(z^2,2,sum)!=0]
  for(k in 1:K){
    if(vpos && sum(abs(v.init[v.init[,k]>0,k]))<sum(abs(v.init[v.init[,k]<0,k]))) v.init[,k] <- -v.init[,k]
    if(vneg && sum(abs(v.init[v.init[,k]<0,k]))<sum(abs(v.init[v.init[,k]>0,k]))) v.init[,k] <- -v.init[,k]
    out <- SparseCCA(xres,zres,v.init[,k],typex,typez,penaltyx, penaltyz,niter,trace, upos, uneg, vpos, vneg,chromx,chromz)
    coef <- out$d
    d <- c(d, coef)
    xres <- rbind(xres, sqrt(coef)*t(out$u))
    zres <- rbind(zres, -sqrt(coef)*t(out$v))
    u <- cbind(u, out$u)
    v <- cbind(v, out$v)
  }
  ubig <- u
  vbig <- v
  if(typex!="ordered"){
    ubig <- matrix(0,nrow=ncol(x),ncol=K)
    ubig[apply(x^2,2,sum)!=0,] <- u
  }
  if(typez!="ordered"){
    vbig <- matrix(0,nrow=ncol(z),ncol=K)
    vbig[apply(z^2,2,sum)!=0,] <- v
  }
  return(list(u=ubig,v=vbig,d=d))
}

fastsvd <- function(x,z){
  # fast svd of t(x)%*%z, where ncol(x)>>nrow(x) and same for z
  xx=x%*%t(x)
  xx2=msqrt(xx)
  y=t(z)%*%xx2
  a=try(svd(y), silent=TRUE)
  iter <- 1
  if(inherits(a, "try-error") && iter<10){
    a=try(svd(y), silent=TRUE)
    iter <- iter+1
  }
  if(iter==10) stop("too many tries.")
  v=a$u
  d=a$d
  zz=z%*%t(z)
  zz2=msqrt(zz)
  y=t(x)%*%zz2
  a=try(svd(y), silent=TRUE)
  iter <- 1
  if(inherits(a, "try-error") && iter<10){
    a=try(svd(y), silent=TRUE)
    iter <- iter+1
  }
  if(iter==10) stop("too many tries.")
  u=a$u
  return(list(u=u,v=v,d=d))
}

msqrt <- function(x){
  eigenx <- eigen(x)
  return(eigenx$vectors%*%diag(sqrt(pmax(0,eigenx$values)))%*%t(eigenx$vectors))
}




SparseCCA <- function(x,y,v,typex,typez,penaltyx, penaltyz,niter,trace, upos, uneg, vpos, vneg,chromx,chromz){
  vold <- rnorm(length(v))
  u <- rnorm(ncol(x))
  for(i in 1:niter){
    if(sum(is.na(u))>0 || sum(is.na(v))>0){
      v <- rep(0, length(v))
      vold <- v
    }
    if(sum(abs(vold-v))>1e-6){
      if(trace) cat(i,fill=F)
      # Update u #
      unew <- rep(NA, ncol(x))
      if(typex=="standard"){
        #argu <- t(x)%*%(y%*%v)
        argu <- matrix(y%*%v,nrow=1)%*%x
        if(upos) argu <- pmax(argu,0)
        if(uneg) argu <- pmin(argu,0)
        lamu <- BinarySearch(argu,penaltyx*sqrt(ncol(x)))
        su <- soft(argu,lamu)
        u <-  matrix(su/l2n(su), ncol=1)
      }else if(typex=="ordered"){
        yv <- y%*%v
        if(is.null(chromx)) chromx <- rep(1, ncol(x))
        for(j in unique(chromx)){
          xyv <- as.numeric(t(yv)%*%x[,chromx==j])#as.numeric(t(x[,chromx==j])%*%yv)
          if(penaltyx!=0){
            coefs <- FLSA(xyv/l2n(xyv),lambda1=penaltyx,lambda2=penaltyx)
#            diagfl.out <- diag.fused.lasso.new(xyv/l2n(xyv), lam1=penaltyx)
#            lam2ind <- which.min(abs(diagfl.out$lam2-penaltyx))
#            coefs <- diagfl.out$coef[,lam2ind]
          }
          if(penaltyx==0){
            coefs <- xyv/l2n(xyv)
          }
          unew[chromx==j] <- coefs
        }
        u <- unew
        if(sum(is.na(u))==0 && sum(abs(u))>0) u <- u/l2n(u)
        u <- matrix(u,ncol=1)
      }
      # Done updating u #
      # Update v #
      vnew <- rep(NA, ncol(y))
      if(typez=="standard"){
        vold <- v
        #argv <- (t(u)%*%t(x))%*%y
        argv <- matrix(x%*%u,nrow=1)%*%y
        if(vpos) argv <- pmax(argv,0)
        if(vneg) argv <- pmin(argv,0)
        lamv <- BinarySearch(argv,penaltyz*sqrt(ncol(y)))
        sv <- soft(argv, lamv)
        v <-  matrix(sv/l2n(sv),ncol=1)
      } else if (typez=="ordered"){
        xu <- x%*%u
        if(is.null(chromz)) chromz <- rep(1, ncol(y))
        for(j in unique(chromz)){
          yxu <- as.numeric(t(xu)%*%y[,chromz==j])#as.numeric(t(y[,chromz==j])%*%xu)
          if(penaltyz!=0){
            coefs <- FLSA(yxu/l2n(yxu),lambda1=penaltyz,lambda2=penaltyz)
#            diagfl.out <- diag.fused.lasso.new(yxu/l2n(yxu), lam1=penaltyz)
#            lam2ind <- which.min(abs(diagfl.out$lam2-penaltyz))
#            coefs <- diagfl.out$coef[,lam2ind]
          }
          if(penaltyz==0){
            coefs <- yxu/l2n(yxu)
          }
          vnew[chromz==j] <- coefs
        }
        v <- vnew
        if(sum(is.na(v))==0 && sum(abs(v))>0) v <- v/l2n(v)
        v <- matrix(v,ncol=1)
      }
      # Done updating v #
    }
  }
  if(trace) cat(fill=T)
  # Update d #
  d <-  sum((x%*%u)*(y%*%v))
  # Done updating d #
  if(sum(is.na(u))>0 || sum(is.na(v))>0){
    u <- matrix(rep(0,ncol(x)),ncol=1)
    v <- matrix(rep(0,ncol(y)),ncol=1)
    d <- 0
  }
  return(list(u=u,v=v,d=d))
}

CheckVs <- function(v,x,z,K){ # If v is NULL, then get v as appropriate.
    ##print(list(v=v, x = x, z = z, K = K))
    if(!is.null(v) && !is.matrix(v)) v <- matrix(v,nrow=ncol(z))
  if(!is.null(v) && ncol(v)<K) v <- NULL
  if(!is.null(v) && ncol(v)>K) v <- matrix(v[,1:K],ncol=K)
  if(is.null(v) && ncol(z)>nrow(z) && ncol(x)>nrow(x)){
    v <- try(matrix(fastsvd(x,z)$v[,1:K],ncol=K), silent=TRUE)
    attempt <- 1
    while(("try-error" %in% class(v))  && attempt < 10){
      v <- try(matrix(fastsvd(x,z)$v[,1:K],ncol=K), silent=TRUE)
      attempt <- attempt+1
    }
    if(attempt==10) stop("Problem computing SVD.")
  } else if (is.null(v) && (ncol(z)<=nrow(z) || ncol(x)<=nrow(x))){
    attempt <- 1
    v <- try(matrix(svd(t(x)%*%z)$v[,1:K],ncol=K), silent=TRUE)
    while(("try-error" %in% class(v)) && attempt<10){
      v <- try(matrix(svd(t(x)%*%z)$v[,1:K],ncol=K), silent=TRUE)
      attempt <- attempt+1
    }
    if(attempt==10) stop("Problem computing SVD.")
  }
  return(v)
}


ftrans <- function(a){
  return(log((1+a)/(1-a)))
}



CCA.permute.both <- function(x,z,typex,typez,penaltyxs,penaltyzs,niter,v,trace,nperms,standardize,chromx,chromz,upos,uneg,vpos,vneg,outcome,y,cens){
  call <- match.call()
  if(standardize){
    x <- scale(x,TRUE,TRUE)
    z <- scale(z,TRUE,TRUE)
  }
  v <- CheckVs(v,x,z,1)
  ccperms=nnonzerous.perms=nnonzerovs.perms=matrix(NA, length(penaltyxs), nperms)
  ccs=nnonzerous=nnonzerovs=numeric(length(penaltyxs))
  for(i in 1:nperms){
    if(trace && .Platform$OS.type!="windows") cat("\n Permutation ",i," out of ", nperms, " ")
#   #  if(trace && .Platform$OS.type=="windows" && i==1) pb <- winProgressBar(title="Doing Permutations", min=0, max=1, initial=(i/nperms))
#     # if(trace && .Platform$OS.type=="windows" && i>1) setWinProgressBar(pb, value=(i/nperms))
    sampz <- sample(1:nrow(z))
    sampx <- sample(1:nrow(x))
    for(j in 1:length(penaltyxs)){
      if(trace && .Platform$OS.type!="windows") cat(j,fill=FALSE)
      if(i==1){
        out <- CCA(x,z,typex=typex,typez=typez,penaltyx=penaltyxs[j], penaltyz=penaltyzs[j],y=y,outcome=outcome,cens=cens,niter=niter,v=v,trace=FALSE, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg, standardize=FALSE,chromz=chromz,chromx=chromx)
        nnonzerous[j] <- sum(out$u!=0)
        nnonzerovs[j] <- sum(out$v!=0)
        if(mean(out$u==0)!=1 && mean(out$v==0)!=1){
          ccs[j] <- cor(x%*%out$u,z%*%out$v)
        } else {
          ccs[j] <- 0
        }
      }
      out <- CCA(x[sampx,],z[sampz,],typex=typex,typez=typez,penaltyx=penaltyxs[j], penaltyz=penaltyzs[j],y=y,outcome=outcome,cens=cens,niter=niter,v=v,trace=FALSE, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg, standardize=FALSE,chromz=chromz,chromx=chromx)
      nnonzerous.perms[j,i] <- sum(out$u!=0)
      nnonzerovs.perms[j,i] <- sum(out$v!=0)
      if(mean(out$u==0)!=1 && mean(out$v==0)!=1){
        ccperms[j,i] <- cor(x[sampx,]%*%out$u,z[sampz,]%*%out$v)
      } else {
        ccperms[j,i] <- 0
      }
    }
  }
#   if(trace && .Platform$OS.type=="windows") close(pb)
  cc.norm <- ftrans(ccs)
  ccperm.norm <- ftrans(ccperms)
  zstats <- (cc.norm - rowMeans(ccperm.norm))/(apply(ccperm.norm,1,sd) + .05)
    # 0.05 added to the denominator to avoid getting zstat of INFINITY
  if(trace) cat(fill=T)
  pvals <- apply(sweep(ccperms,1,ccs,"-")>=0,1,mean)
  results <- list(zstats=zstats,penaltyxs=penaltyxs, penaltyzs=penaltyzs,bestpenaltyx=penaltyxs[which.max(zstats)], bestpenaltyz=penaltyzs[which.max(zstats)], cors=ccs, corperms=ccperms, ft.cors=cc.norm,ft.corperms=rowMeans(ccperm.norm),nnonzerous=nnonzerous,nnonzerovs=nnonzerovs, nnonzerous.perm=rowMeans(nnonzerous.perms),nnonzerovs.perm=rowMeans(nnonzerovs.perms),call=call,v.init=v,pvals=pvals,nperms=nperms,chromz=chromz,chromx=chromx,typex=typex,typez=typez, pvalbestz=pvals[which.max(zstats)])
  return(results)
}


#' @method plot CCA.permute
#' @export
plot.CCA.permute <- function(x,...){
  penaltyxs <- x$penaltyxs
  penaltyzs <- x$penaltyzs
  if(length(penaltyxs)==1 && length(penaltyzs)==1) stop("Cannot plot output of CCA.permute if only 1 tuning parameter was considered.")
  ccs <- x$cors
  nperms <- x$nperms
  zstats <- x$zstats
  ccperms <- x$corperms
  par(mfrow=c(2,1))
  if(length(unique(penaltyxs))==1 && length(unique(penaltyzs))>1){
    plot(penaltyzs, ccs, main="Correlations For Real/Permuted Data", xlab="Penalty on data set 2", ylab="Correlations", ylim=range(ccperms,ccs))
    points(penaltyzs, ccs, type="l")
    for(i in 1:nperms) points(penaltyzs, ccperms[,i], col="green")
    plot(penaltyzs,zstats,main="Z-Statistics", xlab="Penalty on data set 2", ylab="Z-statistic")
    lines(penaltyzs,zstats)
  }
  if(length(unique(penaltyzs))==1 && length(unique(penaltyxs))>1){
    plot(penaltyxs, ccs, main="Correlations For Real/Permuted Data", xlab="Penalty on data set 1", ylab="Correlations", ylim=range(ccperms,ccs))
    points(penaltyxs, ccs, type="l")
    for(i in 1:nperms) points(penaltyxs, ccperms[,i], col="green")
    plot(penaltyxs,zstats,main="Z-Statistics", xlab="Penalty on data set 1", ylab="Z-statistic")
    lines(penaltyxs,zstats)
  }
  if(length(unique(penaltyzs))>1 && length(unique(penaltyxs))>1 && sum(penaltyxs!=penaltyzs)>0){
    plot(1:length(penaltyxs), ccs, main="Correlations For Real/Permuted Data", xlab="Index of Tuning Parameters Considered", ylab="Correlations", ylim=range(ccperms,ccs))
    points(1:length(penaltyxs), ccs, type="l")
    for(i in 1:nperms) points(1:length(penaltyxs), ccperms[,i], col="green")
    plot(1:length(penaltyxs),zstats,main="Z-Statistics", xlab="Index of Tuning Parameters Considered", ylab="Z-statistic")
    lines(1:length(penaltyxs),zstats)
  }
  if(length(unique(penaltyzs))>1 && length(unique(penaltyxs))>1 && sum(penaltyxs!=penaltyzs)==0){
    plot(penaltyxs, ccs, main="Correlations For Real/Permuted Data", xlab="Penalty on data sets 1 and 2", ylab="Correlations", ylim=range(ccperms,ccs))
    points(penaltyxs, ccs, type="l")
    for(i in 1:nperms) points(penaltyxs, ccperms[,i], col="green")
    plot(penaltyxs,zstats,main="Z-Statistics", xlab="Penalty on data sets 1 and 2", ylab="Z-statistic")
    lines(penaltyxs,zstats)
  }
}

CCA.permute.zonly<- function(x,z,typex,typez,penaltyx,penaltyzs,niter,v,trace,nperms,standardize,chromx,chromz,upos,uneg,vpos,vneg,outcome,y,cens){
  call <- match.call()
  if(standardize){
    x <- scale(x,TRUE,TRUE)
    z <- scale(z, TRUE, TRUE)
  }
  v <- CheckVs(v,x,z,1)
  ccperms=nnonzerous.perms=nnonzerovs.perms=matrix(NA, length(penaltyzs), nperms)
  ccs=nnonzerous=nnonzerovs=numeric(length(penaltyzs))
  storevs <- NULL
  for(i in 1:nperms){
        if(trace && .Platform$OS.type!="windows") cat("\n Permutation ",i," out of ", nperms, " ")
#     #if(trace && .Platform$OS.type=="windows" && i==1) pb <- winProgressBar(title="Doing Permutations", min=0, max=1, initial=(i/nperms))
#     #if(trace && .Platform$OS.type=="windows" && i>1) setWinProgressBar(pb, value=(i/nperms))
    sampz <- sample(1:nrow(z))
    sampx <- sample(1:nrow(x))
    for(j in 1:length(penaltyzs)){
      if(trace && .Platform$OS.type!="windows") cat(j,fill=FALSE)
      if(i==1){
        out <- CCA(x,z,typex=typex,typez=typez,penaltyx=penaltyx, penaltyz=penaltyzs[j],y=y,outcome=outcome,cens=cens,niter=niter,v=v,trace=FALSE, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg, standardize=FALSE,chromz=chromz,chromx=chromx)
        nnonzerous[j] <- sum(out$u!=0)
        nnonzerovs[j] <- sum(out$v!=0)
        if(mean(out$u==0)!=1 && mean(out$v==0)!=1){
          ccs[j] <- cor(x%*%out$u,z%*%out$v)
        } else {
          ccs[j] <- 0
        }
        storevs <- cbind(storevs, out$v)
      }
      out <- CCA(x[sampx,],z[sampz,],typex=typex,typez=typez,penaltyx=penaltyx, penaltyz=penaltyzs[j],y=y,outcome=outcome,cens=cens,niter=niter,v=v,trace=FALSE, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg, standardize=FALSE,chromx=chromx,chromz=chromz)
      nnonzerous.perms[j,i] <- sum(out$u!=0)
      nnonzerovs.perms[j,i] <- sum(out$v!=0)
      if(mean(out$u==0)!=1 && mean(out$v==0)!=1){
        ccperms[j,i] <- cor(x[sampx,]%*%out$u,z[sampz,]%*%out$v)
      } else {
        ccperms[j,i] <- 0
      }
    }
  }
#   if(trace && .Platform$OS.type=="windows") close(pb)
  cc.norm <- ftrans(ccs)
  ccperm.norm <- ftrans(ccperms)
  zstats <- (cc.norm - rowMeans(ccperm.norm))/(apply(ccperm.norm,1,sd) + .05)
  if(trace) cat(fill=T)
  pvals <- apply(sweep(ccperms,1,ccs,"-")>=0,1,mean)
  results <- list(zstats=zstats,typex=typex,typez=typez,penaltyxs=rep(penaltyx,length(penaltyzs)),penaltyzs=penaltyzs,bestpenaltyx=penaltyx,bestpenaltyz=penaltyzs[which.max(zstats)], cors=ccs, corperms=ccperms, ft.cors=cc.norm,ft.corperms=rowMeans(ccperm.norm),nnonzerous=nnonzerous,nnonzerovs=nnonzerovs, nnonzerous.perm=rowMeans(nnonzerous.perms),nnonzerovs.perm=rowMeans(nnonzerovs.perms),call=call,v.init=v, pvals=pvals,nperms=nperms,chromx=chromx,chromz=chromz,storevs=storevs, outcome=outcome, pvalbestz=pvals[which.max(zstats)])
  return(results)
}

CCA.permute.justone <- function(x,z,typex,typez,penaltyx,penaltyz,niter,v,trace,nperms,standardize,chromx,chromz,upos,uneg,vpos,vneg,outcome,y,cens){
  call <- match.call()
  if(standardize){
    x <- scale(x,TRUE,TRUE)
    z <- scale(z, TRUE, TRUE)
  }
  v <- CheckVs(v,x,z,1)
  storevs <- NULL
  for(i in 1:nperms){
        if(trace && .Platform$OS.type!="windows") cat("\n Permutation ",i," out of ", nperms, " ")
#     #if(trace && .Platform$OS.type=="windows" && i==1) pb <- winProgressBar(title="Doing Permutations", min=0, max=1, initial=(i/nperms))
#     #if(trace && .Platform$OS.type=="windows" && i>1) setWinProgressBar(pb, value=(i/nperms))
    sampz <- sample(1:nrow(z))
    sampx <- sample(1:nrow(x))
      if(i==1){
        out <- CCA(x,z,typex=typex,typez=typez,penaltyx=penaltyx, penaltyz=penaltyz,y=y,outcome=outcome,cens=cens,niter=niter,v=v,trace=FALSE, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg, standardize=FALSE,chromz=chromz,chromx=chromx)
  ccperms=nnonzerous.perms=nnonzerovs.perms=rep(NA, nperms)
        nnonzerou <- sum(out$u!=0)
        nnonzerov <- sum(out$v!=0)
        if(mean(out$u==0)!=1 && mean(out$v==0)!=1){
          cc <- cor(x%*%out$u,z%*%out$v)
        } else {
          cc <- 0
        }
        storevs <- cbind(storevs, out$v)
      }
      out <- CCA(x[sampx,],z[sampz,],typex=typex,typez=typez,penaltyx=penaltyx, penaltyz=penaltyz,y=y,outcome=outcome,cens=cens,niter=niter,v=v,trace=FALSE, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg, standardize=FALSE,chromx=chromx,chromz=chromz)
      nnonzerous.perms[i] <- sum(out$u!=0)
      nnonzerovs.perms[i] <- sum(out$v!=0)
      if(mean(out$u==0)!=1 && mean(out$v==0)!=1){
        ccperms[i] <- cor(x[sampx,]%*%out$u,z[sampz,]%*%out$v)
      } else {
        ccperms[i] <- 0
      }
    }
#   if(trace && .Platform$OS.type=="windows") close(pb)
  cc.norm <- ftrans(cc)
  ccperm.norms <- ftrans(ccperms)
  zstat <- (cc.norm - mean(ccperm.norms))/(sd(ccperm.norms) + .05)
  if(trace) cat(fill=T)
  cc <- as.numeric(cc)
  ccperms <- as.numeric(ccperms)
  pval <- mean(ccperms>=cc)
  results <- list(zstats=zstat,typex=typex,typez=typez,penaltyxs=penaltyx, penaltyzs=penaltyz,bestpenaltyx=penaltyx,bestpenaltyz=penaltyz, cors=cc, corperms=ccperms, ft.cors=cc.norm,ft.corperms=mean(ccperm.norms),nnonzerous=nnonzerou,nnonzerovs=nnonzerov, nnonzerous.perm=mean(nnonzerous.perms),nnonzerovs.perm=mean(nnonzerovs.perms),call=call,v.init=v, pvals=pval,nperms=nperms,chromx=chromx,chromz=chromz,storevs=storevs, outcome=outcome, pvalbestz=pval)
  return(results)
}


CCA.permute.xonly<- function(x,z,typex,typez,penaltyxs,penaltyz,niter,v,trace,nperms=25,standardize,chromx,chromz,upos,uneg,vpos,vneg,outcome,y,cens){
  call <- match.call()
  if(standardize){
    x <- scale(x,TRUE,TRUE)
    z <- scale(z, TRUE, TRUE)
  }
  v <- CheckVs(v,x,z,1)
  ccperms=nnonzerous.perms=nnonzerovs.perms=matrix(NA, length(penaltyxs), nperms)
  ccs=nnonzerous=nnonzerovs=numeric(length(penaltyxs))
  storevs <- NULL
  for(i in 1:nperms){
    if(trace && .Platform$OS.type!="windows") cat("\n Permutation ",i," out of ", nperms, " ")
#     #if(trace && .Platform$OS.type=="windows" && i==1) pb <- winProgressBar(title="Doing Permutations", min=0, max=1, initial=(i/nperms))
#     #if(trace && .Platform$OS.type=="windows" && i>1) setWinProgressBar(pb, value=(i/nperms))
    sampz <- sample(1:nrow(z))
    sampx <- sample(1:nrow(x))
    for(j in 1:length(penaltyxs)){
      if(trace && .Platform$OS.type!="windows") cat(j,fill=FALSE)
      if(i==1){
        out <- CCA(x,z,typex=typex,typez=typez,penaltyx=penaltyxs[j], penaltyz=penaltyz,y=y,outcome=outcome,cens=cens,niter=niter,v=v,trace=FALSE, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg, standardize=FALSE,chromz=chromz,chromx=chromx)
        nnonzerous[j] <- sum(out$u!=0)
        nnonzerovs[j] <- sum(out$v!=0)
        if(mean(out$u==0)!=1 && mean(out$v==0)!=1){
          ccs[j] <- cor(x%*%out$u,z%*%out$v)
        } else {
          ccs[j] <- 0
        }
        storevs <- cbind(storevs, out$v)
      }
      out <- CCA(x[sampx,],z[sampz,],typex=typex,typez=typez,penaltyx=penaltyxs[j], penaltyz=penaltyz,y=y,outcome=outcome,cens=cens,niter=niter,v=v,trace=FALSE, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg, standardize=FALSE,chromx=chromx,chromz=chromz)
      nnonzerous.perms[j,i] <- sum(out$u!=0)
      nnonzerovs.perms[j,i] <- sum(out$v!=0)
      if(mean(out$u==0)!=1 && mean(out$v==0)!=1){
        ccperms[j,i] <- cor(x[sampx,]%*%out$u,z[sampz,]%*%out$v)
      } else {
        ccperms[j,i] <- 0
      }
    }
  }
#   if(trace && .Platform$OS.type=="windows") close(pb)
  cc.norm <- ftrans(ccs)
  ccperm.norm <- ftrans(ccperms)
  zstats <- (cc.norm - rowMeans(ccperm.norm))/(apply(ccperm.norm,1,sd) + .05)
  if(trace) cat(fill=T)
  pvals <- apply(sweep(ccperms,1,ccs,"-")>=0,1,mean)
  results <- list(zstats=zstats,typex=typex,typez=typez,penaltyxs=penaltyxs,penaltyzs=rep(penaltyz, length(penaltyxs)),bestpenaltyx=penaltyxs[which.max(zstats)],bestpenaltyz=penaltyz, cors=ccs, corperms=ccperms, ft.cors=cc.norm,ft.corperms=rowMeans(ccperm.norm),nnonzerous=nnonzerous,nnonzerovs=nnonzerovs, nnonzerous.perm=rowMeans(nnonzerous.perms),nnonzerovs.perm=rowMeans(nnonzerovs.perms),call=call,v.init=v, pvals=pvals,nperms=nperms,chromx=chromx,chromz=chromz,storevs=storevs, outcome=outcome, pvalbestz=pvals[which.max(zstats)])
  return(results)
}







#' Select tuning parameters for sparse canonical correlation analysis using the
#' penalized matrix decomposition.
#'
#' This function can be used to automatically select tuning parameters for
#' sparse CCA using the penalized matrix decompostion. For each data set x and
#' z, two types are possible: (1) type "standard", which does not assume any
#' ordering of the columns of the data set, and (2) type "ordered", which
#' assumes that columns of the data set are ordered and thus that corresponding
#' canonical vector should be both sparse and smooth (e.g. CGH data).
#'
#' For X and Z, the samples are on the rows and the features are on the
#' columns.
#'
#' The tuning parameters are selected using a permutation scheme. For each
#' candidate tuning parameter value, the following is performed: (1) The
#' samples in X are randomly permuted nperms times, to obtain matrices
#' $X*_1,X*_2,...$. (2) Sparse CCA is run on each permuted data set $(X*_i,Z)$
#' to obtain factors $(u*_i, v*_i)$. (3) Sparse CCA is run on the original data
#' (X,Z) to obtain factors u and v. (4) Compute $c*_i=cor(X*_i u*_i,Z v*_i)$
#' and $c=cor(Xu,Zv)$. (5) Use Fisher's transformation to convert these
#' correlations into random variables that are approximately normally
#' distributed. Let Fisher(c) denote the Fisher transformation of c. (6)
#' Compute a z-statistic for Fisher(c), using
#' $(Fisher(c)-mean(Fisher(c*)))/sd(Fisher(c*))$. The larger the z-statistic,
#' the "better" the corresponding tuning parameter value.
#'
#' This function also gives the p-value for each pair of canonical variates
#' (u,v) resulting from a given tuning parameter value. This p-value is
#' computed as the fraction of $c*_i$'s that exceed c (using the notation of
#' the previous paragraph).
#'
#' Using this function, only the first left and right canonical variates are
#' considered in selection of the tuning parameter.
#'
#' Note that x and z must have same number of rows. This function
#' performs just a one-dimensional search in tuning parameter space,
#' even if penaltyxs and penaltyzs both are vectors: the pairs
#' `(penaltyxs[1],penaltyzs[1])`,
#' `(penaltyxs[2],penaltyzs[2])`,.... are considered.
#' @param x Data matrix; samples are rows and columns are features.
#' @param z Data matrix; samples are rows and columns are features. Note that x
#' and z must have the same number of rows, but may (and generally will) have
#' different numbers of columns.
#' @param typex Are the columns of x unordered (type="standard") or ordered
#' (type="ordered")? If "standard", then a lasso penalty is applied to v, to
#' enforce sparsity. If "ordered" (generally used for CGH data), then a fused
#' lasso penalty is applied, to enforce both sparsity and smoothness.
#' @param typez Are the columns of z unordered (type="standard") or ordered
#' (type="ordered")? If "standard", then a lasso penalty is applied to v, to
#' enforce sparsity. If "ordered" (generally used for CGH data), then a fused
#' lasso penalty is applied, to enforce both sparsity and smoothness.
#' @param penaltyxs The set of x penalties to be considered. If
#' typex="standard", then the L1 bound on u is penaltyxs*sqrt(ncol(x)). If
#' "ordered", then it's the lambda for the fused lasso penalty. The user can
#' specify a single value or a vector of values. If penaltyxs is a vector and
#' penaltyzs is a vector, then the vectors must have the same length. If NULL,
#' then the software will automatically choose a single lambda value if type is
#' "ordered", or a grid of (L1 bounds)/sqrt(ncol(x)) if type is "standard".
#' @param penaltyzs The set of z penalties to be considered. If
#' typez="standard", then the L1 bound on v is penaltyzs*sqrt(ncol(z)). If
#' "ordered", then it's the lambda for the fused lasso penalty. The user can
#' specify a single value or a vector of values. If penaltyzs is a vector and
#' penaltyzs is a vector, then the vectors must have the same length. If NULL,
#' then the software will automatically choose a single lambda value if type is
#' "ordered", or a grid of (L1 bounds)/sqrt(ncol(z)) if type is "standard".
#' @param niter How many iterations should be performed each time CCA is
#' called? Default is 3, since an approximate estimate of u and v is acceptable
#' in this case, and otherwise this function can be quite time-consuming.
#' @param v The first K columns of the v matrix of the SVD of X'Z. If NULL,
#' then the SVD of X'Z will be computed inside this function. However, if you
#' plan to run this function multiple times, then save a copy of this argument
#' so that it does not need to be re-computed (since that process can be
#' time-consuming if X and Z both have high dimension).
#' @param trace Print out progress?
#' @param nperms How many times should the data be permuted? Default is 25. A
#' large value of nperms is very important here, since the formula for
#' computing the z-statistics requires a standard deviation estimate for the
#' correlations obtained via permutation, which will not be accurate if nperms
#' is very small.
#' @param standardize Should the columns of X and Z be centered (to have mean
#' zero) and scaled (to have standard deviation 1)? Default is TRUE.
#' @param chromx Used only if typex="ordered"; a vector of length ncol(x) that
#' allows you to specify which chromosome each CGH spot is on. If NULL, then it
#' is assumed that all CGH spots are on same chromosome.
#' @param chromz Used only if typex="ordered"; a vector of length ncol(z) that
#' allows you to specify which chromosome each CGH spot is on. If NULL, then it
#' is assumed that all CGH spots are on same chromosome.
#' @param upos If TRUE, then require all elements of u to be positive in sign.
#' Default is FALSE. Can only be used if type is standard.
#' @param uneg If TRUE, then require all elements of u to be negative in sign.
#' Default is FALSE. Can only be used if type is standard.
#' @param vpos If TRUE, then require all elements of v to be positive in sign.
#' Default is FALSE.  Can only be used if type is standard.
#' @param vneg If TRUE, then require all elements of v to be negative in sign.
#' Default is FALSE. Can only be used if type is standard.
#' @param outcome If you would like to incorporate a phenotype into CCA
#' analysis - that is, you wish to find features that are correlated across the
#' two data sets and also correlated with a phenotype - then use one of
#' "survival", "multiclass", or "quantitative" to indicate outcome type.
#' Default is NULL.
#' @param y If outcome is not NULL, then this is a vector of phenotypes - one
#' for each row of x and z. If outcome is "survival" then these are survival
#' times; must be non-negative. If outcome is "multiclass" then these are class
#' labels. Default NULL.
#' @param cens If outcome is "survival" then these are censoring statuses for
#' each observation. 1 is complete, 0 is censored. Default NULL.
#' @return \item{zstat}{The vector of z-statistics, one per element of
#' sumabss.} \item{pvals}{The vector of p-values, one per element of sumabss.}
#' \item{bestpenaltyx}{The x penalty that resulted in the highest z-statistic.}
#' \item{bestpenaltyz}{The z penalty that resulted in the highest z-statistic.}
#' \item{cors}{The value of cor(Xu,Zv) obtained for each value of sumabss.}
#' \item{corperms}{The nperms values of cor(X*u*,Zv*) obtained for each value
#' of sumabss, where X* indicates the X matrix with permuted rows, and u* and
#' v* are the output of CCA using data (X*,Z).} \item{ft.cors}{The result of
#' applying Fisher transformation to cors.} \item{ft.corperms}{The result of
#' applying Fisher transformation to corperms.} \item{nnonzerous}{Number of
#' non-zero u's resulting from applying CCA to data (X,Z) for each value of
#' sumabss.} \item{nnonzerouv}{Number of non-zero v's resulting from applying
#' CCA to data (X,Z) for each value of sumabss.} \item{v.init}{The first factor
#' of the v matrix of the SVD of x'z. This is saved in case this function (or
#' the CCA function) will be re-run later.}
#' @seealso \link{PMD},\link{CCA}
#'
#' @references
#' Witten D. M., Tibshirani R.,  and Hastie, T. (2009)
#' \emph{A penalized matrix decomposition, with applications to sparse principal components and canonical correlation analysis}, \emph{Biostatistics, Gol 10 (3), 515-534, Jul 2009}\cr
#' @examples
#'
#' # See examples in CCA function
#'
#' @export CCA.permute
CCA.permute <- function(x,z,typex=c("standard", "ordered"), typez=c("standard","ordered"), penaltyxs=NULL, penaltyzs=NULL, niter=3,v=NULL,trace=TRUE,nperms=25, standardize=TRUE, chromx=NULL, chromz=NULL,upos=FALSE, uneg=FALSE, vpos=FALSE, vneg=FALSE, outcome=NULL, y=NULL, cens=NULL){
  if(ncol(x)<2) stop("Need at least 2 features in data set x.")
  if(ncol(z)<2) stop("Need at least 2 features in data set z.")
  u <- NULL
  typex <- match.arg(typex)
  typez <- match.arg(typez)
  call <- match.call()
  if(!is.null(penaltyxs) && !is.null(penaltyzs) && length(penaltyxs)>1 && length(penaltyzs)>1 && length(penaltyxs)!=length(penaltyzs)) stop("Penaltyxs and Penaltyzs must be same length, or one must have length 1. This is because tuning parameters are considered in pairs.")
  if(is.null(penaltyxs) && typex=="ordered"){
    u <- CheckVs(NULL,z,x,1)
    penaltyxs <- c(ChooseLambda1Lambda2(as.numeric(u)))
    warning("Since type of x is ordered, the penalty for x was chosen w/o permutations.")
  }
  if(is.null(penaltyzs) && typez=="ordered"){
    v <- CheckVs(v,x,z,1)
    penaltyzs <- c(ChooseLambda1Lambda2(as.numeric(v)))
    warning("Since type of z is ordered, the penalty for z was chosen w/o permutations.")
  }
  if(is.null(penaltyxs)) penaltyxs <- seq(.1,.7,len=10)
  if(is.null(penaltyzs)) penaltyzs <- seq(.1,.7,len=10)
  if(typex=="ordered" && (upos||uneg)) stop("If type=ordered then you cannot require elements of u to be positive or negative!")
  if(typez=="ordered" && (vpos||vneg)) stop("If type=ordered then you cannot require elements of v to be positive or negative!")
  if(length(unique(penaltyxs))==1 && length(unique(penaltyzs))==1){
    out <- CCA.permute.justone(x=x,z=z,typex=typex,typez=typez,penaltyx=penaltyxs[1],penaltyz=penaltyzs[1],niter=niter,v=v,trace=trace,nperms=nperms,standardize=standardize,chromx=chromx,chromz=chromz,upos=upos,uneg=uneg,vpos=vpos,vneg=vneg,outcome=outcome,y=y,cens=cens)
  }
  if(length(penaltyxs)==1 && length(penaltyzs)>1) out <- CCA.permute.zonly(x=x,z=z,typex=typex,typez=typez,penaltyx=penaltyxs,penaltyzs=penaltyzs,niter=niter,v=v,trace=trace,nperms=nperms,standardize=standardize,chromx=chromx,chromz=chromz,upos=upos,uneg=uneg,vpos=vpos,vneg=vneg,outcome=outcome,y=y,cens=cens)
  if(length(penaltyxs)>1 && length(penaltyzs)==1) out <- CCA.permute.xonly(x=x,z=z,typex=typex,typez=typez,penaltyxs=penaltyxs,penaltyz=penaltyzs,niter=niter,v=v,trace=trace,nperms=nperms,standardize=standardize,chromx=chromx,chromz=chromz,upos=upos,uneg=uneg,vpos=vpos,vneg=vneg,outcome=outcome,y=y,cens=cens)
  if(length(penaltyzs)>1 && length(penaltyxs)>1) out <- CCA.permute.both(x=x,z=z,typex=typex,typez=typez,penaltyxs=penaltyxs,penaltyzs=penaltyzs,niter=niter,v=v,trace=trace,nperms=nperms,standardize=standardize,chromx=chromx,chromz=chromz,upos=upos,uneg=uneg,vpos=vpos,vneg=vneg,outcome=outcome,y=y,cens=cens)
  out$call <- call
  out$upos <- upos
  out$uneg <- uneg
  out$vpos <- vpos
  out$vneg <- vneg
  class(out) <- "CCA.permute"
  return(out)
}


#' @method print CCA.permute
#' @export
print.CCA.permute <- function(x,...){
  cat("Call: ")
  dput(x$call)
  cat("\n")
  cat("Type of x: ", x$typex,"\n")
  cat("Type of z: ", x$typez,"\n")
  if(x$upos) cat("U's constrained to be positive", fill=TRUE)
  if(x$uneg) cat("U's constrained to be negative", fill=TRUE)
  if(x$vpos) cat("V's constrained to be positive", fill=TRUE)
  if(x$vneg) cat("V's constrained to be negative", fill=TRUE)
  if(!is.null(x$outcome)) cat("Outcome used: ", x$outcome, fill=TRUE)
  if(length(x$penaltyxs)>1 && length(x$penaltyzs)>1){
    tab <- round(cbind(x$penaltyxs,x$penaltyzs, x$zstats,x$pvals,x$cors,rowMeans(x$corperms),x$ft.cors,x$ft.corperms,x$nnonzerous,x$nnonzerovs),3)
#    if(x$typex=="ordered" && x$typez=="ordered") dimnames(tab) <- list(1:length(x$penaltyxs), c("X Lambda", "Z Lambda","Z-Stat","P-Value","Cors","Cors Perm", "FT(Cors)", "FT(Cors Perm)", "# U's non-zero", "# V's non-zero"))
#    if(x$typex=="standard" && x$typez=="ordered") dimnames(tab) <- list(1:length(x$penaltyxs), c("X L1 Bound", "Z Lambda","Z-Stat","P-Value","Cors","Cors Perm", "FT(Cors)", "FT(Cors Perm)", "# U's non-zero", "# V's non-zero"))
#    if(x$typex=="ordered" && x$typez=="standard") dimnames(tab) <- list(1:length(x$penaltyxs), c("X Lambda", "Z L1 Bound","Z-Stat","P-Value","Cors","Cors Perm", "FT(Cors)", "FT(Cors Perm)", "# U's non-zero", "# V's non-zero"))
#    if(x$typex=="standard" && x$typez=="standard") dimnames(tab) <- list(1:length(x$penaltyxs), c("X L1 Bound", "Z L1 Bound","Z-Stat","P-Value","Cors","Cors Perm", "FT(Cors)", "FT(Cors Perm)", "# U's non-zero", "# V's non-zero"))
    dimnames(tab) <- list(1:length(x$penaltyxs), c("X Penalty", "Z Penalty", "Z-Stat", "P-Value", "Cors", "Cors Perm", "FT(Cors)", "FT(Cors Perm)", "# U's Non-Zero", "# Vs Non-Zero"))
    print(tab)
    if(x$typex=="standard") cat("Best L1 bound for x: ", x$bestpenaltyx,fill=TRUE)
    if(x$typex=="ordered") cat("Best lambda for x: ", x$bestpenaltyx,fill=TRUE)
    if(x$typez=="standard") cat("Best L1 bound for z: ", x$bestpenaltyz,fill=TRUE)
    if(x$typez=="ordered") cat("Best lambda for z: ", x$bestpenaltyz,fill=TRUE)
 } else {
   cat("P-value is ", x$pvals, fill=TRUE)
   cat("Z-stat is ", x$zstats, fill=TRUE)
   cat("Correlation is ", x$cors, fill=TRUE)
   cat("Average correlation of permuted data is ", mean(x$corperms),fill=TRUE)
 }
}

CCAPhenotypeZeroSome <- function(x,z,y,qt=.8,cens=NULL,outcome=c("quantitative", "survival", "multiclass"), typex,typez){
  outcome <- match.arg(outcome)
  if(outcome=="quantitative"){
    score.x <- quantitative.func(t(x)[,!is.na(y)],y[!is.na(y)])$tt
    score.z <- quantitative.func(t(z)[,!is.na(y)],y[!is.na(y)])$tt
  } else if (outcome=="survival"){
    score.x <- cox.func(t(x)[,!is.na(y)],y[!is.na(y)],cens[!is.na(y)])$tt
    score.z <- cox.func(t(z)[,!is.na(y)],y[!is.na(y)],cens[!is.na(y)])$tt
  } else if (outcome=="multiclass"){
    score.x <- multiclass.func(t(x)[,!is.na(y)],y[!is.na(y)])$tt
    score.z <- multiclass.func(t(z)[,!is.na(y)],y[!is.na(y)])$tt
  }
  if(typex=="standard"){
    keep.x <- abs(score.x)>=quantile(abs(score.x),qt)
  } else if(typex=="ordered"){
    lam <- ChooseLambda1Lambda2(as.numeric(score.x))
    flsa.out <- FLSA(as.numeric(score.x),lambda1=lam, lambda2=lam)
#    diagfl.out <- diag.fused.lasso.new(as.numeric(score.x), lam1=lam)
#    lam2ind <- which.min(abs(diagfl.out$lam2-lam))
#    flsa.out <- diagfl.out$coef[,lam2ind]
    par(mfrow=c(2,1))
    keep.x <- abs(flsa.out)>=quantile(abs(flsa.out), qt)
    if(mean(keep.x)==1 | mean(keep.x)==0) keep.x <- (abs(score.x) >= quantile(abs(score.x), qt))
  }
  if(typez=="standard"){
    keep.z <- abs(score.z)>=quantile(abs(score.z),qt)
  } else if(typez=="ordered"){
    lam <- ChooseLambda1Lambda2(as.numeric(score.z))
    flsa.out <- FLSA(as.numeric(score.z),lambda1=lam, lambda2=lam)
#    diagfl.out <- diag.fused.lasso.new(as.numeric(score.z), lam1=lam)
#    lam2ind <- which.min(abs(diagfl.out$lam2-lam))
#    flsa.out <- diagfl.out$coef[,lam2ind]
    par(mfrow=c(2,1))
    keep.z <- abs(flsa.out)>=quantile(abs(flsa.out), qt)
    if(mean(keep.z)==1 | mean(keep.z)==0) keep.z <- (abs(score.z) >= quantile(abs(score.z), qt))
  }
  xnew <- x
  xnew[,!keep.x] <- 0
  znew <- z
  znew[,!keep.z] <- 0
  return(list(x=xnew,z=znew))
}

