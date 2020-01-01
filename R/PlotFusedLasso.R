#' Plot CGH data
#'
#' Given a vector of gains/losses at CGH spots, this makes a plot of gain/loss
#' on each chromosome.
#'
#' This function makes a plot of regions of genomic gain/loss.
#'
#' @param array A vector containing the chromosomal location of each CGH spot.
#' @param chrom A numeric vector of the same length as "array"; its values
#' should indicate the chromosome that each CGH spot is on (for instance, for
#' human genomic data, values of chrom should range from 1 to 24). If NULL,
#' then it is assumed that all elements of 'array' are on the same chromosome.
#' @param nuc A numeric vector of same length as "array", indicating the
#' nucleotide position of each CGH spot. If NULL, then the function assumes
#' that each CGH spot corresponds to a consecutive position. E.g. if there are
#' 200 CGH spots on chromosome 1, then they are located at positions
#' 1,2,...,199,200.
#' @param main Give your plot a title.
#' @param scaleEachChrom Default is TRUE. This means that each chromosomes CGH
#' spots are divided by 1.1 times the max of the CGH spots on that chromosome.
#' This way, the CGH spots on each chromosome of the plot are as big as
#' possible (i.e. easy to see). If FALSE, then all of the CGH spots are divided
#' by 1.1 times the max of ALL the CGH spots. This means that on some
#' chromosomes CGH spots might be hard to see, but has the advantage that now
#' relative magnitudes of CGH spots on different chromosomes can be seen from
#' figure.
#' @seealso \link{PMD}, \link{PMD.cv}, \link{CCA}, \link{CCA.permute}
#' @references
#' \insertRef{pmid19377034}{PMA}
#' @examples
#'
#' # Use breast data
#' data(breastdata)
#' attach(breastdata)
#'
#' # dna contains CGH data and chrom contains chromosome of each CGH spot;
#' # nuc contains position of each CGH spot.
#' dna <- t(dna)
#' PlotCGH(dna[1,],chrom=chrom,nuc=nuc,main="Sample 1: All Chromosomes")
#' PlotCGH(dna[1,chrom==1], chrom=chrom[chrom==1], nuc=nuc[chrom==1],
#' main= "Sample 1: Chrom 1")
#' PlotCGH(dna[1,chrom<=3], chrom=chrom[chrom<=3], nuc=nuc[chrom<=3],
#'  main="Sample 1: Chroms 1, 2, and 3")
#' detach(breastdata)
#'
#' @export PlotCGH
PlotCGH <- function(array,chrom=NULL,nuc=NULL,main="", scaleEachChrom=TRUE){
  if(is.null(chrom)){
    chrom <- rep(1,length(array))
    warning("Since chrom was not entered, PlotCGH assumed that all CGH spots in array are on the same chromosome.")
  }
#  if(!is.numeric(chrom)) stop("Chrom must be numeric.")
  if(is.null(nuc)){
    nuc <- rep(NA, length(chrom))
    for(i in unique(chrom)){
      nuc[chrom==i] <- 1:sum(chrom==i)
    }
  }
  scaledarray <- numeric(length(array))
  if(scaleEachChrom){
    for(i in (unique(chrom))) scaledarray[chrom==i] <- array[chrom==i]/(1.1*max(abs(array[chrom==i])))
  } else {
    scaledarray <- array/(.9*max(abs(array)))
  }
  plot.CGH.FL.Single(scaledarray,chrom,nuc,main)
}

plot.CGH.FL.Single<-function(array, chr, nucposi, main=""){
  if(length(array)!=length(chr) || length(array)!=length(nucposi)) stop("Array, chrom, and nuc must all have the same length (or chrom & nuc can be NULL).")
  plot(0,0,type="n",axes=F,ylim=c(0,length(unique(chr))+1),xlim=c(-.05*max(nucposi), max(nucposi)),ylab="",xlab="",main=main,cex.main=1)
  for(j in 1:length((unique(chr)))){
    chrj <- (unique(chr))[j]
    jp=length((unique(chr)))-j+1
    nuc=nucposi[chr==chrj]
    y=array[chr==chrj]
    y[is.na(y)]<-0
    yposi=ynega=y
    yposi[y<0]<-0
    ynega[y>0]<-0
    pick<-(1:length(y))[y!=0]
    if(length(pick)>0){
      segments(nuc[pick],jp,nuc[pick],jp+yposi[pick],col=2)
      segments(nuc[pick],jp,nuc[pick],jp+ynega[pick],col=3)
    }
    segments(0,jp,max(nuc),jp)
    text(-.05*max(nucposi),jp,labels=chrj,cex=.7)
  }
}
