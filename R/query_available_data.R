#' View sites, dates, and variables available in the StreamPULSE database
#'
#' Uses StreamPULSE API to query data via MySQL. To download data, see
#' \code{request_data}.
#'
#' Only records fully encompassing the range defined by startdate and enddate
#' will be returned UNLESS variable is also specified. In that case,
#' any sites for which that variable has been measured at any point during the
#' specified date range will be returned.
#'
#' The data returned by this function depend on which of the input parameters
#' are specified and which are omitted. See examples for all possible ways
#' to use this function.
#'
#' @author Mike Vlah, \email{vlahm13@gmail.com}
#' @param region string representing region. Must be a StreamPULSE region
#'   abbreviation like 'NC' or 'AZ'. Set this
#'   argument to 'all' for a list of all available region and site abbreviations.
#' @param site string representing site. Must be a StreamPULSE site
#'   abbreviation like 'Eno' or 'LV'.
#' @param startdate date string formatted 'YYYY-MM-DD'. Must be specified if
#'   enddate is. Only records fully
#'   encompassing the range defined by startdate and enddate will be returned.
#'   See details for a caveat.
#' @param enddate date string formatted 'YYYY-MM-DD'. Must be specified if
#'   enddate is. Only records fully
#'   encompassing the range defined by startdate and enddate will be returned.
#'   See details for a caveat.
#' @param variable string representing a StreamPULSE variable name, such as
#'   'DO_mgL'. Set this argument to 'all' for a list of all variable names.
#' @param powell logical; set to TRUE to query only data that are part of
#'   the Powell Center metabolism synthesis. Set to FALSE (default) to
#'   query only data from StreamPULSE and NEON. This parameter affects
#'   queries that only include startDate, endDate, and variable. If your
#'   query is of any other form, this parameter will be ignored.
#' @return returns a list containing one or more data frames containing
#'   site, variable, or time data. Which of these data frames is/are returned
#'   depends on which of the input parameters are specified. See examples
#' @seealso \code{\link{request_data}} for downloading StreamPULSE data.
#' @export
#' @examples
#' #View all sites in region
#' query_available_data(region='AZ')
#'
#' #View all sites in all regions
#' query_available_data(region='all')
#'
#' #View all variables at site and available time range for site
#' query_available_data(region='AZ', site='LV')
#'
#' #View all sites at which variable has been sampled
#' query_available_data(variable='Depth_m')
#'
#' #View full list of variables
#' query_available_data(variable='all')
#'
#' #View all sites for which temporal coverage encompasses the specified dates
#' query_available_data(startdate='2016-06-13', enddate='2017-11-01')
#'
#' #View all variables available at site, encompassing date range
#' query_available_data(region='NC', site='Eno', startdate='2017-01-08',
#'     enddate='2017-04-09')
#'
#' #View date range available for variable at site
#' query_available_data(region='NC', site='Eno', variable='WaterTemp_C')
#'
#' #View all sites at which variable has been sampled within date range.
#' #This one takes a while to run.
#' query_available_data(startdate='2017-01-08', enddate='2017-04-09',
#'     variable='DO_mgL')
query_available_data = function(region=NULL, site=NULL, startdate=NULL,
    enddate=NULL, variable=NULL, powell=FALSE){

    #basic checks (more in Flask code)
    if(length(region) > 1 ||
        (! is.null(region) && ! is.character(region))){
        stop("region must be a single string.", call.=FALSE)
    }
    if(length(site) > 1 ||
        (! is.null(site) && ! is.character(site))){
        stop("site must be a single string.", call.=FALSE)
    }
    if(length(variable) > 1 ||
        (! is.null(variable) && ! is.character(variable))){
        stop('variable must be a single string.', call.=FALSE)
    }
    if(length(startdate) > 1 ||
        (! is.null(startdate) &&
            ! grepl('^[0-9]{4}-[0-9]{2}-[0-9]{2}$', startdate))){
        stop('startdate must be a single string of the form: YYYY-MM-DD',
            call.=FALSE)
    }
    if(length(enddate) > 1 ||
        (! is.null(enddate) &&
            ! grepl('^[0-9]{4}-[0-9]{2}-[0-9]{2}$', enddate))){
        stop('enddate must be a single string of the form: YYYY-MM-DD',
            call.=FALSE)
    }
    if(! is.null(site) && is.null(region)){
        stop('Cannot specify site without region.', call.=FALSE)
    }
    if(is.null(site) && ! is.null(region) &&
        (! is.null(variable) || ! is.null(startdate))){
        stop('If region is specified but not site, all other arguments must be NULL',
            call.=FALSE)
    }
    if(! is.null(region) && region == 'all' &&
        (! is.null(variable) || ! is.null(startdate) || ! is.null(site))){
        stop('If region is "all", all other arguments must be NULL',
            call.=FALSE)
    }
    if(! is.null(variable) && variable == 'all' &&
        (! is.null(region) || ! is.null(startdate) || ! is.null(site))){
        stop('If variable is "all", all other arguments must be NULL',
            call.=FALSE)
    }
    if(is.null(startdate) && is.null(region) && is.null(variable)){
        stop('At least one argument must be specified.',
            call.=FALSE)
    }
    if(! is.null(startdate) & ! is.null(enddate)){
        if(as.Date(enddate) < as.Date(startdate)){
            stop("Start date is after end date.", call.=FALSE)
        }
    }
    if(! is.null(startdate) && is.null(enddate)){
        stop("If startdate is supplied, enddate must be too.", call.=FALSE)
    }
    if(is.null(startdate) && ! is.null(enddate)){
        stop("If enddate is supplied, startdate must be too.", call.=FALSE)
    }
    if(! is.null(startdate) && ! is.null(variable) && ! is.null(region)){
        stop("Received arguments for dates, region+site, and variable.\n\t",
            "Omit at least one of these.", call.=FALSE)
    }
    if(! is.null(startdate) && ! is.null(variable) && is.null(region) &&
        ! is.logical(powell)){
        stop("Argument to parameter 'powell' must be logical.")
    }
    if(! is.null(startdate) && ! is.null(variable) && is.null(region) &&
        ! powell){
        cat(paste("NOTE: set powell=TRUE to query Powell Center metabolism",
            "\n\tsynthesis data. If you don't see this message, this note does",
            "\n\tnot apply. It's only necessary here because this particular",
            "\n\tquery takes ages. It would take ages^2 if StreamPULSE and Powell",
            "\n\tdata were queried simultaneously.\n"))
    }

    #assemble url based on user input
    #u = "localhost:5000/query_available_data?"
    u = "https://data.streampulse.org/query_available_data?"
    if(!is.null(region)) u = paste0(u, "&region=", region)
    if(!is.null(site)) u = paste0(u, "&site=", site)
    if(!is.null(startdate)) u = paste0(u, "&startdate=", startdate)
    if(!is.null(enddate)) u = paste0(u, "&enddate=", enddate)
    if(!is.null(variable)) u = paste0(u, "&variable=", variable)
    u = paste0(u, "&powell=", powell)

    cat(paste0('API call: ', u, '\n\n'))
    if(is.null(region) && ! is.null(variable) && ! is.null(startdate)){
        message(paste('NOTE: This particular request takes the server a while,',
            '\n\teven if you abort it from R.'))
    }

    #retrieve json; read into r object
    r = httr::GET(u)
    json = httr::content(r, as="text", encoding="UTF-8")
    d = try(jsonlite::fromJSON(json), silent=TRUE)

    #check for errors
    if(class(d) == 'try-error'){
        stop(paste0('Unable to process request. Please check your\n\t',
            'arguments.'), call.=FALSE)
    }
    if(length(d) == 1 && ! is.null(d$error)){
        stop(d$error, call.=FALSE)
    }

    #organize response data into a neat list
    if('sites' %in% names(d)){
        col_order = c('region', 'site', 'name', 'latitude', 'longitude',
            'usgsGage', 'addDate', 'firstRecord', 'lastRecord',
            'embargoDaysLeft', 'contact', 'contactEmail')
        d$sites = d$sites[order(d$sites$region, d$sites$site),
            as.vector(na.omit(match(col_order, colnames(d$sites))))]
        d$sites$addDate = as.POSIXct(d$sites$addDate, tz="UTC",
            format="%a, %d %b %Y %T GMT")
        d$sites$firstRecord = as.POSIXct(d$sites$firstRecord, tz="UTC",
            format="%a, %d %b %Y %T GMT")
        d$sites$lastRecord = as.POSIXct(d$sites$lastRecord, tz="UTC",
            format="%a, %d %b %Y %T GMT")
        d$sites$usgsGage = as.numeric(d$sites$usgsGage)
    }
    if('variables' %in% names(d)){
        d$variables = data.frame('variables'=sort(strsplit(d$variables,
            ',')[[1]]))
    }
    if('datebounds' %in% names(d)){
        dates = gsub('\\[\\[\\\"(.*?Z)\\\",\\\"(.*?Z).*', '\\1,\\2',
            d$datebounds)
        dates = as.POSIXct(strsplit(dates, ',')[[1]], tz="UTC",
            format="%Y-%m-%dT%H:%M:%S.000Z")
        d$datebounds = data.frame(firstRecord=dates[1], lastRecord=dates[2])
    }

    return(d)
}
