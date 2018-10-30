#' View a list of results of streamMetabolizer
#'
#' Uses StreamPULSE API to query data via MySQL.
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
#' @return returns a \code{list} containing one or more data frames containing
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
query_available_results = function(region=NULL, site=NULL, year=NULL){

    #basic checks (more in Flask code)
    if(length(region) > 1 ||
        (! is.null(region) && ! is.character(region))){
        stop("region must be a single string.", call.=FALSE)
    }
    if(length(site) > 1 ||
        (! is.null(site) && ! is.character(site))){
        stop("site must be a single string.", call.=FALSE)
    }
    if(! is.null(year)){
        year = as.character(year)
    }
    if(length(year) > 1){
        stop('year must be a single string.', call.=FALSE)
    }
    if(! is.null(site) && is.null(region)){
        stop('Cannot specify site without region.', call.=FALSE)
    }
    if(! is.null(region) && region == 'all' &&
        (! is.null(site) || ! is.null(year))){
        stop('If region is "all", all other arguments must be NULL',
            call.=FALSE)
    }
    if(is.null(year) && is.null(region) && is.null(site)){
        stop('At least one argument must be specified.',
            call.=FALSE)
    }
    if(! is.null(year) && ! is.null(site) && ! is.null(region)){
        stop("Too many arguments.", call.=FALSE)
    }

    #assemble url based on user input
    u = "localhost:5000/query_available_results?"
    # u = "https://data.streampulse.org/query_available_results?"
    if(!is.null(region)) u = paste0(u, "&region=", region)
    if(!is.null(site)) u = paste0(u, "&site=", site)
    if(!is.null(year)) u = paste0(u, "&year=", year)
    cat(paste0('\nAPI call: ', u, '\n\n'))

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

    #convert results to dataframe in a list
    df = data.frame('region'=d[[1]][,1], 'site'=d[[1]][,2],
        'year'=d[[1]][,3])
    df = df[order(df$region, df$site, df$year),]
    rownames(df) = 1:nrow(df)
    d$available_model_results = df

    return(d)
}
