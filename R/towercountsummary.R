#'Tower Count Summary
#'@description Reads raw tower count data and computes the abundance estimate
#'  and standard error, incorporating all specified interpolations.
#'
#'  Count days are considered as beginning at 16:00 and ending at 15:59, with
#'  counts occurring between 16:00 and 23:59 on a given calendar date treated as
#'  occurring on the following count day.
#'
#'  Hourly run proportions are calculated by summing the absolute counts (or
#'  counts) across days for each hour block, and dividing by the total of
#'  absolute counts (or counts.)  Summation is performed for days in which no
#'  hourly count is missed, unless otherwise specified.
#'
#'  If all counts were represented during a count day (having values of clarity
#'  recorded as less than or equal to \code{maxclarity}), daily escapement is
#'  estimated as
#'
#'  N_d = (M_d/m_d)*sum(y_dj) for j=1,...,m_d
#'
#'  in which the subscript d denotes day, M_d is the total number of possible
#'  paired 10-minute counting periods in day d, m_d is the number of sampled
#'  10-minute counting periods in day d, and y_dj is the observed count for
#'  period j within day d.
#'
#'  Variance of the daily total is estimated as
#'
#'  V(N_d) = (1-m_d/M_d)*(M_d^2)*(s_d^2)/m_d
#'
#'  in which
#'
#'  s_d^2 = 1/(2*(m_d-1)) * sum((y_dj - y_d(j-1))^2) for j=2,...,m_d
#'
#'  If some counts are missing due to water quality, but more than 25 percent of
#'  the expected daily run (calculated using the hourly run proportions,
#'  expressing the diel migratory pattern) and 25 percent or more of the peak
#'  run periods are represented in the daily count, the quantity
#'
#'  y_dc,interp = y_dc*(1-p_edp)/p_edp
#'
#'  is added to the daily count, in which y_dc denotes daily count, and p_edp
#'  represents the proportion of the expected daily escapement successfully
#'  counted.  If this quantity is greater than zero, the difference is allocated
#'  to the hourly counts such that the hourly proportions are consistent with
#'  the diel migratory pattern.  If observed hourly counts are greater than
#'  interpolated hourly values, observed counts are retained.
#'
#'  In these cases, the variance of the expanded daily escapement is calculated
#'  as given above, but with m_d calculated as
#'
#'  m_d = m_d,observed*p_edp
#'
#'  This inflates the resulting variance according to the diminished effective
#'  daily sample size.
#'
#'  If daily counts represent less than 25 percent of the expected daily
#'  escapement or if less than 25 percent of the peak run periods are
#'  represented in the daily count, daily totals are interpolated according to
#'  the form
#'
#'  N_i = sum(I(day j was sampled)*N_j)/sum(I(day j was sampled)) for
#'  j=(i-k),...,(i+k)
#'
#'  in which k = the number of days missed due to adverse viewing conditions.
#'  The variance of the daily total is estimated as the maximum daily variance
#'  of the days used for interpolation.
#'
#'  The total estimated escapement and associated variance are then calculated
#'  according to the sums
#'
#'  N = sum(N_d) for d = all possible days
#'
#'  V(N) = sum(V(N_d)) for d = all possible days
#'@param date A vector of dates.  This can be in numeric or date format.
#'@param hour A vector of hour blocks corresponding to counts.  Values must be
#'  numeric, between 0 and 23.
#'@param count A vector of counts.
#'@param clarity A vector of clarity values of water clarity, defined as
#'
#'  1: Excellent / 2: Good / 3: Fair / 4: Poor / 4.5: Very poor / 5:
#'  Unobservable
#'@param maxclarity The largest value (poorest visibility) of clarity that will
#'  be accepted as a count, without interpolating.  Defaults to 4.
#'@param hourfrac The proportion of each hour block that counts were conducted.
#'  Defaults to 1/6.
#'@param interpwithin Whether to interpolate hourly counts with poor clarity
#'  within an otherwise good-clarity count day.  Specifying TRUE imputes
#'  interpolated hourly counts, and specifying FALSE accepts the recorded
#'  counts.  In both cases, the bad-visibility counts are treated as
#'  no-counts in the variance calculation.  Defaults to TRUE.
#'@param interpbetween Whether to interpolate daily counts in which less than 25
#'  percent of the estimated daily run is recorded due to poor visibility.
#'  Specifying TRUE imputes interpolated daily counts, and specifying FALSE
#'  accepts the recorded counts.  Defaults to TRUE.
#'@param hourproprows A vector of count-day rows to use for calculating hourly
#'  run proportions.  If the default (NULL) is accepted, all rows not missing
#'  hourly counts will be used.
#'@param abscounts4props Whether to use absolute counts for calculating hourly
#'  run proportions (TRUE) or treat negative counts (fish running down-river) as
#'  negative in the summation (FALSE).  Defaults to TRUE.
#'@return A countsummary object, with elements \itemize{
#'  \item{\code{$maxclarity}: }{The largest value of clarity that was accepted
#'  as a count} \item{\code{$countmatrix}: }{The matrix of hourly counts that
#'  was used for escapement and variance calculations, including any
#'  interpolations} \item{\code{$rawcountmatrix}: }{The matrix of observed
#'  hourly counts} \item{\code{$claritymatrix}: }{The matrix of clarity values,
#'  corresponding to the counts in \code{$countmatrix} and
#'  \code{$rawcountmatrix}} \item{\code{$propbyhour}: }{A numeric vector of
#'  hourly run proportions, defined as the proportion for hours 16 through 15}
#'  \item{\code{$calcmatrix}: }{The matrix of calculations leading up to the
#'  variance calculation (returned as a data frame)} \item{\code{$Nhat}: }{The
#'  total escapement estimate} \item{\code{$Varhat}: }{The variance of the total
#'  escapement estimate} \item{\code{$SE}: }{The standard error of the total
#'  escapement estimate} \item{\code{$CI}: }{An approximate Wald-type 95 percent
#'  confidence interval for the total escapement} \item{\code{$maxpass}: }{A
#'  two-element vector, defined as the first and last hours of the period of
#'  peak passage.  This is calculated as the hours corresponding to the minimum
#'  time interval representing greater than 80 percent of the overall run.  If
#'  multiple intervals of the same length had greater than 80 percent of the
#'  total run, the median interval was selected.} \item{\code{$totcount}: }{The
#'  total net escapement visually counted.} }
#'@author Matt Tyers
#'@examples
#' data(Gulk_2015)
#'
#' chin2015 <- towercountsummary(date=Gulk_2015$date, hour=Gulk_2015$hour,
#'      count=Gulk_2015$chin, clarity=Gulk_2015$clarity)
#'
#' str(chin2015)
#'@export
towercountsummary <- function(date, hour, count, clarity, maxclarity=4,hourfrac=1/6,interpwithin=T,interpbetween=T,hourproprows=NULL,abscounts4props=T) {
  # checking inputs
  if(!any(class(date)==c("Date","numeric"))) stop("date must be date or numeric")
  if(!is.numeric(hour)) stop("hour must be numeric")
  if(any(floor(hour)!=hour) | min(hour)<0 | max(hour)>23) stop("hour must be integer values between 0 and 23")
  if(!is.numeric(clarity)) stop("clarity must be numeric")
  for(i in 1:length(clarity)) {
    if(is.na(clarity[i])) stop("Missing value in clarity argument")
    if(!any(clarity[i]==c(1,2,3,4,4.5,5))) stop("Invalid value in clarity argument")
  }
  countday <- date
  countday[hour>=16] <- countday[hour>=16] + 1   # counts after 1559 occur on the next count day
  rawcount <- count   # saving a vector of raw counts before removing bad-visibility counts
  count[clarity>maxclarity] <- NA  # removing bad-visibility counts

  # setting up the count and clarity matrices
  rownames <- sort(unique(countday))
  colnames <- c(16:23,0:15)
  countmatrix <- matrix(NA,nrow=length(rownames),ncol=24)
  dimnames(countmatrix)[[1]] <- as.character(rownames)
  dimnames(countmatrix)[[2]] <- colnames
  rawcountmatrix <- countmatrix
  claritymatrix <- countmatrix

  # making the count and clarity matrices
  for(i in 1:length(rownames)) {
    for(j in 1:24) {
      if(length(count[countday==rownames[i] & hour==colnames[j]])>0) {
        countmatrix[i,j] <- count[countday==rownames[i] & hour==colnames[j]]
        rawcountmatrix[i,j] <- rawcount[countday==rownames[i] & hour==colnames[j]]
        claritymatrix[i,j] <- clarity[countday==rownames[i] & hour==colnames[j]]
      }
    }
  }

  # proportions by hour
  propbyhour <- NA
  abscountmatrix <- abs(countmatrix)
  if(!is.null(hourproprows)) {
    allclear <- rep(F,length(rownames))
    allclear[hourproprows] <- T
  }
  if(is.null(hourproprows)) {
    allclear <- NA
    for(i in 1:length(rownames)) allclear[i] <- all(!is.na(countmatrix[i,]))
  }
  if(!abscounts4props) {
    for(j in 1:24) propbyhour[j] <- sum(countmatrix[allclear,j],na.rm=T)/sum(countmatrix[allclear,],na.rm=T)
  }
  if(abscounts4props) {
    for(j in 1:24) propbyhour[j] <- sum(abscountmatrix[allclear,j],na.rm=T)/sum(abscountmatrix[allclear,],na.rm=T)  # i like this one
  }

  # determining the hours of maximum (>80%) passage
  width <- 1
  maxpass <- NA
  # starting with a width of 1 hour, checking all possible intervals and increasing the width by one hour
  while(is.na(maxpass[1])) {
    bigenough <- FALSE
    for(i in 1:(24-width+1)) {
      pass <- sum(propbyhour[i:(i+width-1)])
      bigenough[i] <- (pass>=.8)
    }
    if(any(bigenough)) {
      firstone <- floor(median((1:length(bigenough))[bigenough]))
      maxpass <- c(firstone,(firstone+width-1))
    }
    width <- width+1
  }

  # indices of first and last counts
  first <- last <- NA
  j<-1
  while(is.na(first)) {
    if(!all(is.na(countmatrix[1,(1:j)]))) first <- j
    j<-j+1
  }
  j<-24
  while(is.na(last)) {
    if(!all(is.na(countmatrix[(dim(countmatrix)[1]),(j:24)]))) last <- j
    j<-j-1
  }

  # abundance estimation for each row, including interpolations
  sumYdi <- proprun <- md <- NA
  for(i in 1:(dim(countmatrix)[1])) {

    # proportion of expected daily run for each hour
    proprun[i] <- sum(propbyhour[!is.na(countmatrix[i,])])
    if(i==1) proprun[1] <- sum(propbyhour[!is.na(countmatrix[1,])])/sum(propbyhour[first:24])
    if(i==dim(countmatrix)[1]) proprun[dim(countmatrix)[1]] <- sum(propbyhour[!is.na(countmatrix[(dim(countmatrix)[1]),])])/sum(propbyhour[1:last])

    # interpolation for missing counts within a given day
    if(interpwithin) {
      if(i!=1 & i!=dim(countmatrix)[1]) {  # if i is not the first or last row
        if(proprun[i]>=.25 & any(is.na(countmatrix[i,]))) {  # if interpolation is needed
          if(length(!is.na(countmatrix[i,(maxpass[1]:maxpass[2])]))>=0.25*(maxpass[2]-maxpass[1])) { # if more than 25% of peak periods were counted
            Ydi_i <- sum(countmatrix[i,],na.rm=T)/proprun[i]#*(1-proprun[i])
            toallocate <- Ydi_i - sum(countmatrix[i,],na.rm=T)  # difference in actual and expected counts to allocate to missing counts
            # doing the allocation
            # countmatrix[i,][is.na(countmatrix[i,])] <- ifelse(sum(propbyhour[is.na(countmatrix[i,])])==0,rep(0,length(countmatrix[i,][is.na(countmatrix[i,])])),
            #                                                   round(toallocate*propbyhour[is.na(countmatrix[i,])]/sum(propbyhour[is.na(countmatrix[i,])])))
            if(sum(propbyhour[is.na(countmatrix[i,])])<=0) {
              countmatrix[i,][is.na(countmatrix[i,])] <- rep(0,length(countmatrix[i,][is.na(countmatrix[i,])]))
            }
            if(sum(propbyhour[is.na(countmatrix[i,])])>0) {
              countmatrix[i,][is.na(countmatrix[i,])] <- round(toallocate*propbyhour[is.na(countmatrix[i,])]/sum(propbyhour[is.na(countmatrix[i,])]))
            }
          }
          # if raw count is bigger than interpolated count, accepting the raw count
          for(j in 1:24) countmatrix[i,j] <- ifelse(any(!is.na(c(countmatrix[i,j],rawcountmatrix[i,j]))),max(countmatrix[i,j],rawcountmatrix[i,j],na.rm=T),NA)
        }
      }
      if(i==1) {  # first row
        if(proprun[i]>=.25 & any(is.na(countmatrix[i,first:24]))) {
          if(length(!is.na(countmatrix[i,(maxpass[1]:maxpass[2])]))>=0.25*(maxpass[2]-maxpass[1])) { # if more than 25% of peak periods were counted
            Ydi_i <- sum(countmatrix[i,],na.rm=T)/proprun[i]#*(1-proprun[i])
            toallocate <- Ydi_i - sum(countmatrix[i,],na.rm=T)
            # countmatrix[i,first:24][is.na(countmatrix[i,first:24])] <- ifelse(sum(propbyhour[first:24][is.na(countmatrix[i,first:24])])==0,rep(0,length(countmatrix[i,first:24][is.na(countmatrix[i,first:24])])),
            #                                                                   round(toallocate*propbyhour[first:24][is.na(countmatrix[i,first:24])]/sum(propbyhour[first:24][is.na(countmatrix[i,first:24])])))
            if(sum(propbyhour[first:24][is.na(countmatrix[i,first:24])])<=0) {
              countmatrix[i,][is.na(countmatrix[i,])] <- rep(0,length(countmatrix[i,first:24][is.na(countmatrix[i,first:24])]))
            }
            if(sum(propbyhour[first:24][is.na(countmatrix[i,first:24])])>0) {
              countmatrix[i,][is.na(countmatrix[i,])] <- round(toallocate*propbyhour[first:24][is.na(countmatrix[i,first:24])]/sum(propbyhour[first:24][is.na(countmatrix[i,first:24])]))
            }
          }
          for(j in first:24) countmatrix[i,j] <- ifelse(any(!is.na(c(countmatrix[i,j],rawcountmatrix[i,j]))),max(countmatrix[i,j],rawcountmatrix[i,j],na.rm=T),NA)
        }
      }
      if(i==dim(countmatrix)[1]){  # last row
        if(proprun[i]>=.25 & any(is.na(countmatrix[i,1:last]))) {
          if(length(!is.na(countmatrix[i,(maxpass[1]:maxpass[2])]))>=0.25*(maxpass[2]-maxpass[1])) { # if more than 25% of peak periods were counted
            Ydi_i <- sum(countmatrix[i,],na.rm=T)/proprun[i]#*(1-proprun[i])
            toallocate <- Ydi_i - sum(countmatrix[i,],na.rm=T)
            # countmatrix[i,1:last][is.na(countmatrix[i,1:last])] <- ifelse(sum(propbyhour[1:last][is.na(countmatrix[i,1:last])])==0,rep(0,length(countmatrix[i,1:last][is.na(countmatrix[i,1:last])])),
            #                                                               round(toallocate*propbyhour[1:last][is.na(countmatrix[i,1:last])]/sum(propbyhour[1:last][is.na(countmatrix[i,1:last])])))
            if(sum(propbyhour[1:last][is.na(countmatrix[i,1:last])])<=0) {
              countmatrix[i,][is.na(countmatrix[i,])] <- rep(0,length(countmatrix[i,1:last][is.na(countmatrix[i,1:last])]))
            }
            if(sum(propbyhour[1:last][is.na(countmatrix[i,1:last])])>0) {
              countmatrix[i,][is.na(countmatrix[i,])] <- round(toallocate*propbyhour[1:last][is.na(countmatrix[i,1:last])]/sum(propbyhour[1:last][is.na(countmatrix[i,1:last])]))
            }
          }
          for(j in 1:last) countmatrix[i,j] <- ifelse(any(!is.na(c(countmatrix[i,j],rawcountmatrix[i,j]))),max(countmatrix[i,j],rawcountmatrix[i,j],na.rm=T),NA)
        }
      }
    }
    if(!interpwithin) countmatrix[i,] <- rawcountmatrix[i,]  # if told not to interpolate, accepting raw counts for all

    # expansion
    sumYdi[i] <- sum(countmatrix[i,],na.rm=T)
    md[i] <- 24*proprun[i]
  }

  # expansion
  md[1] <- (24-first+1)*proprun[1]  # special case for first row
  md[dim(countmatrix)[1]] <- last*proprun[dim(countmatrix)[1]] # special case for last row
  Yd <- sumYdi/hourfrac

  # interpolation between days if counts were missed for full days
  if(interpbetween) {
    interprow <- proprun<0.25            # which rows to interpolate (less than 25% of the expected daily run)
    for(i in 1:dim(countmatrix)[1]) {    # which rows to interpolate (less than 25% of the peak period represented)
      if(length(!is.na(countmatrix[i,(maxpass[1]:maxpass[2])]))<0.25*(maxpass[2]-maxpass[1])) interprow[i] <- TRUE
    }

    # some weird stuff to determine how many days were missed in each block of days missed
    countup <- 1*interprow[1]
    for(i in 2:length(interprow)) {
      if(!interprow[i]) countup[i] <- 0
      if(interprow[i]) countup[i] <- countup[i-1]+1
    }
    ninterp <- countup
    for(i in (length(interprow)-1):1) {
      if(countup[i]!=0) ninterp[i] <- max(ninterp[i:(i+1)])
      if(countup[i]==0) ninterp[i] <- 0
    }

    # doing the interpolation
    Ydinterp <- Yd
    Ydinterp[interprow] <- NA
    for(i in 1:length(Ydinterp)) {
      if(interprow[i]) {
        Yd[i] <- mean(Ydinterp[(max(c(1,(i-ninterp[i])))):(min(c((i+ninterp[i]),length(Ydinterp))))],na.rm=T)
        countmatrix[i,] <- rawcountmatrix[i,]
      }
    }

  }

  # variance calculations
  RDS <- NA*countmatrix[,1:23]
  for(j in 1:23) RDS[,j] <- (countmatrix[,(j+1)]-countmatrix[,(j)])^2
  sumRDS <- NA
  for(i in 1:(dim(countmatrix)[1])) sumRDS[i] <- sum(RDS[i,],na.rm=T)
  ssqd2d <- sumRDS/(2*(md-1))
  f2d <- md/(24/hourfrac)
  rowvar <- ifelse(md==0, 0, ((24/hourfrac)^2)*(1-f2d)*ssqd2d/md)
  if(interpbetween) { # in cases of interpolations between rows
    rowvarinterp <- rowvar
    rowvarinterp[interprow] <- NA
    for(i in 1:length(Ydinterp)) {
      if(interprow[i]) {
        rowvar[i] <- max(rowvarinterp[(max(c(1,(i-ninterp[i])))):(min(c((i+ninterp[i]),length(Ydinterp))))],na.rm=T)
      }
    }
  }

  # some summary output
  rowse <- sqrt(rowvar)
  rowcv <- rowse/Yd

  calcmatrix <- data.frame(sumYdi,proprun,md,Yd,sumRDS,ssqd2d,f2d,rowvar,rowse,rowcv,interprow,ninterp)
  row.names(calcmatrix) <- row.names(claritymatrix)

  Nhat <- sum(Yd)
  Varhat <- sum(rowvar)
  SE <- sqrt(Varhat)

  # more summary output
  CI <- Nhat+c(-1.96,1.96)*SE
  totcount <- sum(rawcountmatrix,na.rm=TRUE)

  stufftoreturn <- list(maxclarity=maxclarity,countmatrix=countmatrix,rawcountmatrix=rawcountmatrix,claritymatrix=claritymatrix,propbyhour=propbyhour,calcmatrix=calcmatrix,Nhat=Nhat,Varhat=Varhat,SE=SE,CI=CI,maxpass=maxpass,totcount=totcount)
  class(stufftoreturn) <- "countsummary"
  return(stufftoreturn)
}
