#'Plot Count Matrix
#'@description Creates a plot of the count matrix (see \link{towercountsummary}) in which daily counts are color-coded by water clarity.  Any hourly interpolations that differ from the raw counts are highlighted in red.
#'
#'Also included in the plot are empirical and cumulative run density by date (shown at right) and hour (shown at bottom).  The minimum continuous period containing more than 80 percent of the run is highlighted in green on the hourly density plot at the bottom, and on the hour labels at the top.
#'
#'At bottom right, the abundance estimate and SE are given, as are an approximate 95 percent Wald confidence interval.
#'@param countsummary A countsummary object from \link{towercountsummary}.
#'@param cex The character expansion factor to use for the counts.  Defaults to 0.6.
#'@param cexsummary The character expansion factor to use for the summary at bottom right.  Defaults to 0.7.
#'@param showsummary Whether to show the summary values for the abundance estimate.  Defaults to TRUE.
#'@param ... Additional plotting arguments.
#'@author Matt Tyers
#'@examples
#' data(Gulk_2015, Gulk_2014)
#'
#' chin2015 <- towercountsummary(date=Gulk_2015$date, hour=Gulk_2015$hour,
#'      count=Gulk_2015$chin, clarity=Gulk_2015$clarity)
#' plotcountmatrix(chin2015)
#'
#' sock2015 <- towercountsummary(date=Gulk_2015$date, hour=Gulk_2015$hour,
#'      count=Gulk_2015$sock, clarity=Gulk_2015$clarity)
#' plotcountmatrix(sock2015)
#'
#' chin2014 <- towercountsummary(date=Gulk_2014$date, hour=Gulk_2014$hour,
#'      count=Gulk_2014$chin, clarity=Gulk_2014$clarity)
#' plotcountmatrix(chin2014)
#'@export
plotcountmatrix <- function(countsummary,cex=.6,cexsummary=.7,showsummary=TRUE,showinterp=TRUE,...) {
  if(class(countsummary)!="countsummary") stop("Input must be a complete countsummary object from towercountsummary().")

  rownames <- dimnames(countsummary$countmatrix)[[1]]
  colnames <- dimnames(countsummary$countmatrix)[[2]]

  if(!any(countsummary$calcmatrix$interprow)) showinterp <- FALSE

  # plotting the clarity matrix
  colors <- c("#BBBBFF22","#BBBBFF55","#BBBBFF88","#BBBBFFBB","#FF444477","#FF333399")
  clarities <- c(1:4,4.5,5)
  countsummary$claritymatrix[countsummary$claritymatrix==5] <- 6
  countsummary$claritymatrix[countsummary$claritymatrix==4.5] <- 5
  plot(NA,xlim=c(-4,30+3*showinterp),ylim=c(min((-dim(countsummary$claritymatrix)[1]-6),-12),2),xaxt='n',yaxt='n',xlab="",ylab="",...=...)
  for(i in 1:(dim(countsummary$claritymatrix)[1])) {
    for(j in 1:24) {
      rect(j-1,-i,j,-i+1,border=NA,col=colors[countsummary$claritymatrix[i,j]])
      if(!is.na(countsummary$rawcountmatrix[i,j]) & !is.na(countsummary$countmatrix[i,j])) {
        if(countsummary$rawcountmatrix[i,j]==countsummary$countmatrix[i,j]&!countsummary$calcmatrix$interprow[i]) text((j-.5),(-i+.5),labels=countsummary$rawcountmatrix[i,j],cex=cex)
        if(countsummary$rawcountmatrix[i,j]!=countsummary$countmatrix[i,j]&!countsummary$calcmatrix$interprow[i]) {
          text((j-.5),(-i+.5),labels=paste0(countsummary$rawcountmatrix[i,j],"->",countsummary$countmatrix[i,j]),cex=cex,col="red")
        }
        if(countsummary$calcmatrix$interprow[i]) text((j-.5),(-i+.5),labels=countsummary$rawcountmatrix[i,j],cex=cex,col="gray45")
      }
    }
  }

  # highlighting max passage hours
  rect(countsummary$maxpass[1]-1,1,countsummary$maxpass[2],2,col="#BBFFBBBB",border=NA)
  rect(countsummary$maxpass[1]-1,(-7-dim(countsummary$claritymatrix)[1]),countsummary$maxpass[2],(-1-dim(countsummary$claritymatrix)[1]),col="#BBFFBBBB",border=NA)
  text(mean(countsummary$maxpass)-.5,(-8-dim(countsummary$claritymatrix)[1]),
       labels=paste0(round(100*sum(countsummary$propbyhour[(countsummary$maxpass[1]):(countsummary$maxpass[2])]),1),"% of total run"),cex=cex)

  # color legend
  for(i in 1:6) {
    rect(-4,-i,-3,-i+1,border=NA,col=colors[i])
    text(-3.5,-i+.5,labels=clarities[i],cex=cex)
  }
  text(-3.5,1.5,"Clarity",cex=cex)

  # row names (hours) and column names (dates)
  text(seq(.5,23.5,by=1),1.5,labels=colnames,cex=cex)
  text(-1.5,seq(-.5,(-dim(countsummary$claritymatrix)[1]),by=-1),labels=rownames,cex=cex)

  # expanded total
  text(25,1.5,"Expanded",cex=cex)
  text(25,seq(-.5,(-dim(countsummary$claritymatrix)[1]),by=-1),round(countsummary$calcmatrix$Yd),cex=cex,
       font=1+countsummary$calcmatrix$interprow,col=1+countsummary$calcmatrix$interprow)

  # daily run plot
  rect(26+3*showinterp,(-.5),30+3*showinterp,(-dim(countsummary$claritymatrix)[1]))
  lines(26+3*showinterp+(4*countsummary$calcmatrix$Yd/max(countsummary$calcmatrix$Yd,na.rm=T)),seq(-.5,(-dim(countsummary$claritymatrix)[1]),by=-1),col="blue",lty=2)
  YdCDF <- NA
  for(i in 1:dim(countsummary$claritymatrix)[1]) {
    YdCDF[i] <- sum(countsummary$calcmatrix$Yd[1:i])/sum(countsummary$calcmatrix$Yd)
  }
  lines(26+3*showinterp+(4*YdCDF),seq(-.5,(-dim(countsummary$claritymatrix)[1]),by=-1),col="blue",lwd=2)

  # daily hour plot
  rect(0,(-7-dim(countsummary$claritymatrix)[1]),24,(-1-dim(countsummary$claritymatrix)[1]))
  lines(seq(.5,23.5,by=1),((-7-dim(countsummary$claritymatrix)[1])+(6*countsummary$propbyhour/max(countsummary$propbyhour))),col="blue",lty=2)
  propCDF <- NA
  for(i in 1:24) {
    propCDF[i] <- sum(countsummary$propbyhour[1:i])
  }
  lines(seq(.5,23.5,by=1),((-7-dim(countsummary$claritymatrix)[1])+(6*propCDF)),col="blue",lwd=2)

  # summary numbers
  if(showsummary) {
    text(26+3*showinterp,-2-length(rownames),labels="N est:",cex=cexsummary,pos=4)
    text(26+3*showinterp,-4-length(rownames),labels="SE:",cex=cexsummary,pos=4)
    text(26+3*showinterp,-6-length(rownames),labels="95% CI:",cex=cexsummary,pos=4)
    text(27.4+3*showinterp,-2-length(rownames),labels=round(countsummary$Nhat),cex=cexsummary,pos=4)
    text(27.4+3*showinterp,-4-length(rownames),labels=round(countsummary$SE),cex=cexsummary,pos=4)
    text(27.4+3*showinterp,-6-length(rownames),labels=paste0("(",round(countsummary$CI[1]),", ",round(countsummary$CI[2]),")"),cex=cexsummary,pos=4)
  }

  # interpolation lines
  if(showinterp) {
    whichrows <- ((1:length(rownames))[countsummary$calcmatrix$interprow])
    delt <- seq(0,2,by=2/length(whichrows))
    for(i in 1:length(whichrows)) {
      lines(rep(26+delt[i],2),c(-whichrows[i],1-whichrows[i]),lwd=3,lend=1,col=2)
      whichtointerp1 <- (1:length(rownames))[(1:length(rownames))>=(whichrows[i]-countsummary$calcmatrix$ninterp[whichrows[i]]) & (1:length(rownames))<(whichrows[i]) & !countsummary$calcmatrix$interprow]
      whichtointerp2 <- (1:length(rownames))[(1:length(rownames))<=(whichrows[i]+countsummary$calcmatrix$ninterp[whichrows[i]]) & (1:length(rownames))>(whichrows[i]) & !countsummary$calcmatrix$interprow]
      for(j in 1:length(whichtointerp1)) lines(rep(26+delt[i],2),c(-whichtointerp1[j],1-whichtointerp1[j]),lwd=2,lend=1)
      for(j in 1:length(whichtointerp2)) lines(rep(26+delt[i],2),c(-whichtointerp2[j],1-whichtointerp2[j]),lwd=2,lend=1)
      lines(rep(26+delt[i],2),c(1-whichrows[i],-whichtointerp1[length(whichtointerp1)]),lty=3)
      lines(rep(26+delt[i],2),c(-whichrows[i],1-whichtointerp2[1]),lty=3)
    }
  }
}
