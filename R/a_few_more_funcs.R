#'Save Summary
#'@description Saves output from a countsummary object as a set of .csv files,
#'  in the current working directory.
#'
#'  \itemize{ \item{\code{NAME_counts_calcs.csv: }}{All counts (including
#'  interpolations) and calculations leading to expanded estimate and variance
#'  calculations.}
#'  \item{\code{NAME_rawcounts.csv: }}{All recorded visual counts.}
#'  \item{\code{NAME_clarity.csv: }}{A matrix of clarity values
#'  corresponding to the count matrices.}
#'  \item{\code{NAME_hourprops.csv:
#'  }}{The run proportion and cumulative run proportion by hour, as well as the
#'  interval considered to be the peak daily run.}
#'  \item{\code{NAME_summary.csv:
#'  }}{Summary values: the escapement estimate and its associated variance,
#'  standard error, and 95 percent Wald confidence interval endpoints; also the
#'  total net number of fish visually counted, the number of periods
#'  successfully counted, the number of periods missed due to poor visibility,
#'  the total number of count days, and the number of complete count days with no
#'  compromised counts.} }
#'@param countsummary A countsummary object from \link{towercountsummary}.
#'@param summaryname The name to use in place of \code{NAME} in the file names
#'  given above.  Accepting the default (\code{NULL}) results in the name of the
#'  countsummary object being used.
#'@author Matt Tyers
#'@examples
#' data(Gulk_2015)
#'
#' chin2015 <- towercountsummary(date=Gulk_2015$date, hour=Gulk_2015$hour,
#'      count=Gulk_2015$chin, clarity=Gulk_2015$clarity)
#' savesummary(chin2015,summaryname="Gulk_chin_2015")
#'@export
savesummary <- function(countsummary,summaryname=NULL) {
  if(is.null(summaryname)) summaryname <- deparse(substitute(countsummary))

  counts_calcs <- cbind(countsummary$countmatrix,countsummary$calcmatrix)
  write.csv(counts_calcs,file=paste0(summaryname,"_counts_calcs.csv"))

  write.csv(countsummary$rawcountmatrix,file=paste0(summaryname,"_rawcounts.csv"))
  write.csv(countsummary$claritymatrix,file=paste0(summaryname,"_clarity.csv"))

  cumul <- NA
  for(i in 1:24) cumul[i] <- sum(countsummary$propbyhour[1:i])
  peak <- rep(FALSE,24)
  peak[countsummary$maxpass[1]:countsummary$maxpass[2]] <- TRUE
  hourprop <- data.frame(c(16:23,0:15),countsummary$propbyhour,cumul,peak)
  names(hourprop) <- c("hour","proportion","cumulative","peak_run")
  write.csv(hourprop,file=paste0(summaryname,"_hourprops.csv"))

  allclear <- NA
  for(i in 1:dim(countsummary$countmatrix)[1]) {
    allclear[i] <- all(countsummary$claritymatrix[i,]<=countsummary$maxclarity,na.rm=T) &
      all(!is.na(countsummary$rawcountmatrix[i,]))
  }
  summarynums <- c(countsummary$Nhat,countsummary$Varhat,countsummary$SE,
                   countsummary$CI[1],countsummary$CI[2],countsummary$totcount,
                   sum(countsummary$claritymatrix<=countsummary$maxclarity,na.rm=T),
                   sum(countsummary$claritymatrix>countsummary$maxclarity,na.rm=T),
                   dim(countsummary$countmatrix)[1],sum(allclear))
  summarynames <- c("Nhat","Vhat","SE","95% CI low","95% CI high","Net escapement visually counted",
                    "Number of periods successfully counted",
                    "Number of periods missed due to poor visibility",
                    "Total number of count days",
                    "Complete count days with no compromised counts")
  summarystuff <- data.frame(summarynums,row.names=summarynames)
  names(summarystuff)<-"value"
  write.csv(summarystuff,file=paste0(summaryname,"_summary.csv"))
  
  counts <- data.frame(vis_count=rowSums(countsummary$rawcountmatrix,na.rm=T), 
                       expanded_interpolated=round(countsummary$calcmatrix$Yd), 
                       cumulative=round(cumsum(countsummary$calcmatrix$Yd)))
  write.csv(counts,file=paste0(summaryname,"_counts.csv"))

  cat("files written to",getwd())
}


#'Converting From a Decimal Datetime
#'@description Separates the date and hour components of a single column of decimal datetime.
#'@param datetime A vector of decimal datetime
#'@param origin A character string specifying the origin of the date numbering system.  This must be structured as "yyyy-mm-dd".  Defaults to \code{"1899-12-30"}
#'@return A two-element list:
#'\itemize{
#'\item{\code{$date: }}{Dates in calendar date format, see \link{as.Date}.}
#'\item{\code{$hour: }}{Hour block recorded as an integer between 0 and 23.}}
#'@author Matt Tyers
#'@examples
#' data(datetime2015)
#' head(datetime2015)
#'
#' date_hour <- decdatetime(datetime2015)
#' head(date_hour$date)
#' head(date_hour$hour)
#'@export
decdatetime <- function(datetime,origin="1899-12-30") {
  time <- round((datetime-floor(datetime))*24)
  date <- as.Date(floor(datetime),origin=as.Date(origin))
  OUT <- list(date=date,hour=time)
  return(OUT)
}
