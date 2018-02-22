#' Prepare StreamPULSE data for metabolism modeling
#'
#' Formats the output of \code{\link{request_data}} for stream metabolism model
#' of choice. Filters flagged data and imputes missing data. Acquires/estimates
#' additional variables if necessary.
#'
#' \code{BASE} and \code{streamMetabolizer}, the two metabolism modeling
#' platforms available via StreamPULSE, require different data input formats.
#' Formatting also varies depending on whether one is using a Bayesian framework
#' or MLE. This function supplements and rearranges the raw output of
#' \code{\link{request_data}} to prepare it for a desired set of model
#' specifications.
#'
#' Both \code{BASE} and \code{streamMetabolizer} require dissolved oxygen (DO)
#' concentration, water temperature, and light (PAR) data. If light is missing,
#' it will automatically be estimated based on solar angle. In addition to these
#' variables,
#' \code{streamMetabolizer} requires DO % saturation and depth, and
#' \code{BASE} requires atmospheric pressure. If DO % saturation is missing,
#' it will be calculated automatically from DO concentration, water temperature,
#' and atmospheric pressure. In turn, atmospheric pressure estimates will
#' be automatically retrieved from NOAA (NCEP), if missing,
#' for sites within the U.S.A.
#'
#' If \code{streamMetabolizer} is being used and
#' \code{type='bayes'}, discharge data are also required. These can be estimated
#' from depth using a rating curve (specifications for which can be made via
#' \code{zq_curve}.
#'
#' All single-station models assume that, where applicable, variables
#' represent averages
#' throughout an area delineated by the width of the stream and the approximate
#' oxygen turnover distance. More on this and other considerations can be
#' found by clicking the "Before modeling stream metabolism..." button on
#' \url{https://data.streampulse.org}.
#'
#' The between-sample interval specified by the \code{interval} parameter is
#' also determined programmatically
#' for each variable within \code{d}, and is assumed to be the mode
#' if the between-sample interval varies within a series. If the
#' between-sample interval varies across variables, the longest interval is
#' used for the whole dataset. If the user-specified \code{interval}
#' (i.e. the argument to this parameter) is a multiple of the programmatically
#' determined interval, the dataset will be quietly coerced to the desired
#' interval.
#' This is useful for thinning extremely long datasets in order to
#' avoid out-of-memory errors while running models.
#' If the desired interval is not a multiple of the true sample interval
#' (i.e. the programmatically determined interval), gaps will be introduced to
#' the dataset and the user will be warned. If user-specified and
#' programmatically-determined intervals are identical, nothing will be
#' changed.
#'
#' @author Mike Vlah, \email{vlahm13@gmail.com}
#' @author Aaron Berdanier
#' @param d the output of \code{request_data},
#'   or a \code{list} of \code{data.frame}s so organized.
#' @param model either 'streamMetabolizer' (the default) or 'BASE'. If 'BASE',
#'   \code{type} must be set to \code{'bayes'}.
#' @param type either 'mle' or 'bayes'. If \code{model='BASE'}, \code{type}
#'   must be set to \code{'bayes'}.
#' @param interval a string specifying the between-sample time interval of the
#'   dataset, or the interval to which the dataset should be coerced. Must be
#'   of the form '<number> <unit>', as in '15 min' (the default). Unit can be
#'   'min' or 'hour'. Non-integer hours are tolerated, but minutes must be
#'   specified as integers. See details.
#' @param rm_flagged a list containing any of 'Interesting', 'Questionable',
#'   and 'Bad Data'. Any data points flagged with these specified tags will be
#'   removed (replaced with NA), and then imputed according to \code{fillgaps}.
#'   If data for a selected site and timespan have been cleaned using
#'   \url{https://data.streampulse.org/clean}, it is a good idea to remove any
#'   data points flagged as "Questionable" or "Bad Data". Set this argument to
#'   'none' to keep all flagged data points. Defaults to
#'   \code{list('Questionable', 'Bad Data')}.
#' @param fillgaps a string specifying one of the imputation methods available
#'   to imputeTS::na.seasplit, namely: 'interpolation', 'locf', 'mean',
#'   'random', 'kalman', or 'ma'. May also be 'none', though gaps are not
#'   tolerated by either \code{BASE} or \code{streamMetabolizer}. The specified
#'   imputation method will be attempted after seasonal decomposition.
#'   Periodicity depends on the between-sample interval, and is determined
#'   programmatically (see details for \code{interval}). If the desired
#'   imputation method
#'   fails, which sometimes occurs when series consist largely of NAs, basic
#'   linear interpolation will be performed instead and the user will be
#'   notified.
#' @return returns a list containing two \code{data.frame}s.
#' @seealso \code{\link{prep_metabolism}} for organizing data returned by this
#' @export
#' @examples
#' streampulse_data = request_data(sitecode='NC_Eno',
#'     startdate='2016-06-10', enddate='2016-10-23')
