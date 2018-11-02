#' Retrieve model results from the StreamPULSE server
#'
#' Uses StreamPULSE API to retrieve \code{.rds} objects containing the outputs
#' of streamMetabolizer
#' runs for specific siteyears. To view available model results, see
#' \code{query_available_results}. Visit the data portal at
#' \url{http://data.streampulse.org:3838/streampulse_diagnostic_plots/}
#' to visualize model results via a browser.
#'
#' Request results for a single region, site, and year.
#' Some sites are "embargoed," meaning results from those sites are
#' kept private and can only be accessed by authorized users.
#' If you are authorized, you can access your embargoed results using
#' a unique token, which currently you can only receive by emailing StreamPULSE
#' developer Mike Vlah (\email{vlahm13@gmail.com}).
#'
#' @author Mike Vlah, \email{vlahm13@gmail.com}
#' @param sitecode underscore-separated region and site code, e.g. 'NC_Eno'.
#'   Full list of regions and site codes available at
#'   \url{https://data.streampulse.org/sitelist}. Or, you can use the
#'   \code{query_available_data} function in this package.
#' @param year string or numeric representing year, e.g. '2015' or 2015.
#' @param token a unique alphanumeric string for each registered user of
#'   StreamPULSE. Only necessary for accessing embargoed results. Email
#'   StreamPULSE developer Mike Vlah (\email{vlahm13@gmail.com}) to receive
#'   your token.
#' @return returns a list containing an unfortunately limited collection of
#'   information about the "best" model result we have on record for the
#'   site and year requested. If the output of \code{streamMetabolizer}'s
#'   \code{metab}
#'   function is a variable called \code{x}, this function will give you x@fit,
#'   x@data, and x@data_daily. Likewise, if the output of \code{StreamPULSE}'s
#'   fit_metabolism is a variable called \code{y}, this function will give you
#'   y$fit@fit, y$fit@data, and y$fit@data_daily. In the future, this function
#'   will return a lot more relevant information, but at the moment only these
#'   three list elements are stored on our server for each "best" model run.
#'
#'   Model "bestness" is determined by an automatic
#'   comparison of five criteria each time a new model is fit: 1) proportion of
#'   daily GPP estimates that go negative; 2) proportion of daily ER estimates
#'   that go positive; 3) correlation between daily ER and K600;
#'   4) maximum daily K600; and 5) temporal coverage. A score is assigned
#'   for each of these criteria, and an aggregate score is determined.
#' @seealso \code{\link{query_available_results}} for determining which models
#'  are available for download.
#' @export
#' @examples
#' res = query_available_results(sitecode='FL_NR1000', year=2016)
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
    u = paste0("https://data.streampulse.org/request_results?sitecode=",
        sitecode, "&year=", year)
    # u = paste0("localhost:5000/request_results?sitecode=", sitecode,
    #     "&year=", year)
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
