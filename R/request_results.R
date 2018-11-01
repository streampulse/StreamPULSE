#' Retrieve model results from the StreamPULSE server
#'
#' Uses StreamPULSE API to retrieve \code{.rds} objects from streamMetabolizer
#' runs for specific siteyears. To view available model results, see
#' \code{query_available_results}. Visit the data portal at
#' \url{http://data.streampulse.org:3838/streampulse_diagnostic_plots/}
#' to visualize model results via a browser.
#'
#' Request data for a single region and site within the StreamPULSE database.
#' All sites are embargoed for 1 year or more from the time they first appear
#' in the database. Embargoed data are only available to the user who submits
#' them, and can be accessed via a unique token.
#'
#' @author Mike Vlah, \email{vlahm13@gmail.com}
#' @author Aaron Berdanier
#' @param sitecode underscore-separated region and site code, e.g. 'NC_Eno'.
#'   Full list of regions and site codes available at
#'   \url{https://data.streampulse.org/sitelist}. Or, you can use the
#'   \code{query_available_data} function in this package.
#' @param startdate date string formatted 'YYYY-MM-DD', representing the first
#'   day of data to be requested.
#'   If data coverage does not extend this far back in time, the argument
#'   will be adjusted to the first day with available data. Omit this argument to
#'   include all records up to the first available. To see the range of available
#'   dates for a particular site, use the
#'   \code{query_available_data} function in this package.
#' @param enddate date string formatted 'YYYY-MM-DD', representing the last
#'   day of data to be requested.
#'   If data coverage does not extend this far in time, the argument
#'   will be adjusted to the last day with available data. Omit this argument to
#'   include all records up to the last available. To see the range of available
#'   dates for a particular site, use the
#'   \code{query_available_data} function in this package.
#' @param variables character vector of variable names to request. To see which
#'   variables are available for a given site, use the
#'   \code{query_available_data} function in this package.
#'   Omit this argument to
#'   request all variables potentially useful for metabolism modeling: c('DO_mgL','DOsat_pct','satDO_mgL','WaterPres_kPa',
#'   'Depth_m','WaterTemp_C','Light_PAR','AirPres_kPa','Discharge_m3s').
#' @param token a unique alphanumeric string for each registered user of
#'   StreamPULSE. Only necessary for accessing embargoed data.
#' @return returns a \code{list} containing two \code{data.frame}s.
#'   The first contains all requested
#'   data. The second contains site metadata. Variable names and corresponding
#'   values are stacked in two long columns. Use \code{prep_metabolism} to
#'   format the output of this function.
#' @seealso \code{\link{prep_metabolism}} for organizing data returned by this
#'   function.
#' @export
#' @examples
#' query_available_data(region='all')
#'
#' streampulse_data = request_data(sitecode='NC_Eno',
#'     startdate='2016-06-10', enddate='2016-10-23')
request_results = function(sitecode, year, token=NULL){

    #basic checks (more in Flask code)
    if(length(sitecode) > 1 || ! is.character(sitecode)){
        stop("sitecode must be a single string.", call.=FALSE)
    }
    if(! grepl('\\w+_\\w+', sitecode)){
        stop("sitecode must be something like REGION_SITE, e.g. 'NC_Eno'.",
            call.=FALSE)
    }
    year = as.character(year)
    if(length(year) > 1){
        stop('year must be a single string.', call.=FALSE)
    }
    if(! grepl('[0-9]{4}', year)){
        stop('year must be 4 numerical digits.', call.=FALSE)
    }

    #assemble url based on user input
    # u = paste0("https://data.streampulse.org/request_results?sitecode=", sitecode)
    u = paste0("localhost:5000/request_results?sitecode=", sitecode,
        "&year=", year)
    cat(paste0('\nAPI call: ', u, '\n\n'))

    #retrieve raw rds binary from response
    if(is.null(token)){
        r = httr::GET(u)
    }else{
        r = httr::GET(u, httr::add_headers(Token = token))
    }
    raw = try(httr::content(r), silent=TRUE)

    #check for errors
    if(class(raw) == 'try-error'){
        stop(paste0('Unable to process request. Please check your\n\t',
            'arguments.'), call.=FALSE)
    }
    if(length(raw) == 1 && ! is.null(raw$error)){
        stop(raw$error, call.=FALSE)
    }

    tmp_filename = tempfile('rdata_temp', fileext='.rds')
    con = file(tmp_filename, "wb") #open temp file
    writeBin(raw, con) #write raw rds binary to temp file
    close(con)
    mod_results = readRDS(tmp_filename)

    return(mod_results)
}
