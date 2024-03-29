% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/request_data.R
\name{request_data}
\alias{request_data}
\title{Retrieve data from the StreamPULSE database}
\usage{
request_data(
  sitecode,
  startdate = NULL,
  enddate = NULL,
  variables = NULL,
  token = NULL
)
}
\arguments{
\item{sitecode}{underscore-separated region and site code, e.g. 'NC_Eno'.
Full list of regions and site codes available at
\url{https://data.streampulse.org/sitelist}. Or, you can use the
\code{query_available_data} function in this package.}

\item{startdate}{date string formatted 'YYYY-MM-DD', representing the first
day of data to be requested.
If data coverage does not extend this far back in time, the argument
will be adjusted to the first day with available data. Omit this argument to
include all records up to the first available. To see the range of available
dates for a particular site, use the
\code{query_available_data} function in this package.}

\item{enddate}{date string formatted 'YYYY-MM-DD', representing the last
day of data to be requested.
If data coverage does not extend this far in time, the argument
will be adjusted to the last day with available data. Omit this argument to
include all records up to the last available. To see the range of available
dates for a particular site, use the
\code{query_available_data} function in this package.}

\item{variables}{character vector of variable names to request. To see which
variables are available for a given site, use the
\code{query_available_data} function in this package.
Omit this argument to
request all variables potentially useful for metabolism modeling: c('DO_mgL','DOsat_pct','satDO_mgL','WaterPres_kPa',
'Depth_m','WaterTemp_C','Light_PAR','AirPres_kPa','Discharge_m3s').}

\item{token}{a unique alphanumeric string for each registered user of
StreamPULSE. Only necessary for accessing embargoed data.}
}
\value{
returns a \code{list} containing two \code{data.frame}s.
  The first contains all requested
  data. The second contains site metadata. Variable names and corresponding
  values are stacked in two long columns. Use \code{prep_metabolism} to
  format the output of this function.
}
\description{
Uses StreamPULSE API to query data via MySQL. Visit the data portal at
\url{https://data.streampulse.org} to download data via a browser.
}
\details{
Request data for a single region and site within the StreamPULSE database.
All sites are embargoed for 1 year or more from the time they first appear
in the database. Embargoed data are only available to the user who submits
them, and can be accessed via a unique token.
}
\examples{
query_available_data(region='all')

streampulse_data = request_data(sitecode='NC_Eno',
    startdate='2016-06-10', enddate='2016-10-23')
}
\seealso{
\code{\link{prep_metabolism}} for organizing data returned by this
  function.
}
\author{
Mike Vlah, \email{vlahm13@gmail.com}

Aaron Berdanier
}
