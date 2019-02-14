#' View metabolism model results available on the StreamPULSE server.
#'
#' Uses StreamPULSE API to query data. To download model results, see
#' \code{request_results}.
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
#' @param year string or numeric representing year, e.g. '2015' or 2015.
#' @return returns a list containing a single data frame with region, site,
#'   and year columns.
#' @seealso \code{\link{request_results}} for downloading StreamPULSE model
#'   outputs.
#' @export
#' @examples
#' View all available model results from region
#' query_available_results(region='NC')
#'
#' View all available model results from all regions
#' query_available_results(region='all')
#'
#' View all available model results from site
#' query_available_results(region='NC', site='Eno')
#'
#' View all available model results from year
#' query_available_results(year=2017)
#'
#' View all available model results from region and year
#' query_available_results(region='NC', year=2017)
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
    #u = "localhost:5000/query_available_results?"
    u = "https://data.streampulse.org/query_available_results?"
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
