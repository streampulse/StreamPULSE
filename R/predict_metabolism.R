#' Predict metabolism from a fitted model
#'
#' A thin wrapper for \code{streamMetabolizer}'s \code{predict_metab} function.
#' Generates estimates of GPP, ER, and K600.
#'
#' See documentation for \code{streamMetabolizer::predict_metab}.
#'
#' @author Mike Vlah, \email{vlahm13@gmail.com}
#' @author Aaron Berdanier
#' @param model_fit the output of \link{fit_metabolism}.
#' @return returns the output of \code{streamMetabolizer}'s \code{predict_metab}
#'   function.
#' @seealso \code{\link{request_data}} for acquiring StreamPULSE data;
#'   \code{\link{prep_metabolism}} for organizing data and acquiring additional
#'   variables; \code{\link{fit_metabolism}} for fitting models.
#' @export
#' @examples
#' streampulse_data = request_data(sitecode='NC_Eno',
#'     startdate='2016-06-10', enddate='2016-10-23')
#'
#' fitdata = prep_metabolism(d=streampulse_data, type='bayes',
#'     model='streamMetabolizer', interval='15 min',
#'     rm_flagged=list('Bad Data', 'Questionable'), fillgaps=fillgaps,
#'     zq_curve=list(sensor_height=NULL, Z=Z_data, Q=Q_data,
#'     fit='power', plot=TRUE), estimate_areal_depth=TRUE)
#'
#' modelfit = fit_metabolism(fitdata)
#'
#' predictions = predict_metabolism(modelfit)
predict_metabolism = function(model_fit){
    if(class(model_fit)=="BASE"){
        directory = model_fit$output_directory
        read.csv(paste0(directory,"/output/BASE_results.csv")) %>%
            separate(File, c("fileX", "date", "extX"), "_|\\.") %>%
            select(-fileX, -extX) %>% mutate(date=as.Date(date))
    }else{
        predictions = predict_metab(model_fit)
        return(predictions)
    }
}
