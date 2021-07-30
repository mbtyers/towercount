#'Escapement Estimates for River Tower Counts
#' @description Provides escapement estimates with associated variance for river tower counts.  Also included are functions to visualize and write summary tables.
#' @author Matt Tyers
#'
#' Maintainer: Matt Tyers <matt.tyers@@alaska.gov>
#' @details Version: 0.1
#'
#' Date: 2015-10-23
#'
#' License: GPL (>=2)
#' @docType package
#' @name towercount
#' @aliases towercount-package
NULL


#' Dataset: Datetime 2015
#'
#' A vector of decimal datetimes from a 2015 count
#'
#' @docType data
#' @keywords datasets
#' @name datetime2015
#' @usage data(datetime2015)
#' @format numeric
NULL


#' Dataset: Gulkana Counts 2015
#'
#' A data frame of tower counts from 2015.  This was a very good year for water quality.
#'
#' @docType data
#' @keywords datasets
#' @name Gulk_2015
#' @usage data(Gulk_2015)
#' @format A data frame with the following fields:
#' \itemize{
#' \item{\code{$date: }}{Count dates in calendar date format}
#' \item{\code{$hour: }}{Hour blocks for each count}
#' \item{\code{$clarity: }}{Clarity values for each count}
#' \item{\code{$chin: }}{Chinook salmon counts for each hour block}
#' \item{\code{$sock: }}{Sockeye salmon counts for each hour block}}
NULL


#' Dataset: Gulkana Counts 2014
#'
#' A data frame of tower counts from 2014.  There were some substantial high water events this year, compromising counts.
#'
#' @docType data
#' @keywords datasets
#' @name Gulk_2014
#' @usage data(Gulk_2014)
#' @format A data frame with the following fields:
#' \itemize{
#' \item{\code{$date: }}{Count dates in calendar date format}
#' \item{\code{$hour: }}{Hour blocks for each count}
#' \item{\code{$clarity: }}{Clarity values for each count}
#' \item{\code{$chin: }}{Chinook salmon counts for each hour block}
#' \item{\code{$sock: }}{Sockeye salmon counts for each hour block}}
NULL
