#' Fit stream metabolism models
#'
#' Passes data and model parameters to \code{streamMetabolizer} or
#' \code{BASE} according to arguments passed to \link{prep_metabolism}.
#' NOTE: support for modeling with
#' \code{BASE} is currently in development. Please use \code{streamMetabolizer}
#' in the meantime.
#'
#' This function is a wrapper for \code{streamMetabolizer}'s model specification,
#' (\code{mm_name, specs}) fitting (\code{metab}), and prediction
#' (\code{predict_metab}) functions, and all parameters except for \code{d}
#' are \code{streamMetabolizer} parameters.
#' Default arguments passed to this function
#' specify the recommended starting point for fitting a metabolism model to
#' data from a typical stream or river. In addition to the parameters
#' documented here, this function calculates \code{K600_lnQ_nodes_centers}
#' (parameter of \code{streamMetabolizer}'s \code{specs} function) as
#' a sequence of length 7 from
#' the minimum to the maximum natural log of daily discharge.
#'
#' To access the full model specification interface of \code{streamMetabolizer},
#' please call its corresponding functions
#' (\code{mm_name, specs, metab, predict_metab}) directly.
#' See the example section below for more.
#'
#' @author Mike Vlah, \email{vlahm13@gmail.com}
#' @author Aaron Berdanier
#' @param d the output of \link{prep_metabolism}.
#' @return returns a list containing the output of \code{streamMetabolizer}'s
#'   \code{metab} function (the model fit object), the output of
#'   \code{streamMetabolizer}'s
#'   \code{predict_metab} function (metabolism predictions), and additional
#'   information about model performance and specifications.
#' @param pool_K600 character. Should the model pool information among
#'   days to get more consistent daily estimates for K600? Options (see Details
#'   section of \code{streamMetabolizer}'s \code{mm_name} function for more):
#'   \itemize{
#'     \item \code{none}: no pooling of K600
#'     \item \code{binned}: \eqn{K600 ~ N(B[Q_bin], sigma)} where \eqn{mu ~
#'     N(mu_mu, mu_sigma)} and \eqn{sigma ~ N(sigma_mu, sigma_sigma)}
#'   }
#' @param err_obs_iid logical. Should IID observation error be included? If not,
#'   the model will be fit to the differences in successive DO measurements,
#'   rather than to the DO measurements themselves.
#' @param err_proc_acor logical. Should autocorrelated process error (with the
#'   autocorrelation term phi fitted) be included?
#' @param err_proc_iid logical. Should IID process error be included?
#' @param ode_method character. The method to use in solving the ordinary
#'   differential equation for DO. Options:
#'   \itemize{
#'     \item \code{euler}, formerly \code{Euler}: the final change in DO from
#'     t=1 to t=2 is solely a function of GPP, ER, DO, etc. at t=1
#'     \item \code{trapezoid}, formerly \code{pairmeans}: the final change in DO
#'     from t=1 to t=2 is a function of the mean values of GPP, ER, etc. across
#'     t=1 and t=2.
#'     \item for \code{type='mle'}, options also include \code{rk2} and any
#'     character method accepted by \code{\link[deSolve]{ode}} in the
#'     \code{deSolve} package (\code{lsoda}, \code{lsode}, \code{lsodes},
#'     \code{lsodar}, \code{vode}, \code{daspk}, \code{rk4}, \code{ode23},
#'     \code{ode45}, \code{radau}, \code{bdf}, \code{bdf_d}, \code{adams},
#'     \code{impAdams}, and \code{impAdams_d}; note that many of these have not
#'     been well tested in the context of \code{streamMetabolizer} models)
#'   }
#' @param deficit_src character. From what DO estimate (observed or modeled)
#'   should the DO deficit be computed? Options:
#'   \itemize{
#'     \item \code{DO_mod}: the DO deficit at time t will be {(DO.sat(t) -
#'     DO_mod(t))}, the difference between the equilibrium-saturation value and
#'     the current best estimate of the true DO concentration at that time
#'     \item \code{DO_obs}: the DO deficit at time t will be {(DO.sat(t) -
#'     DO.obs(t))}, the difference between the equilibrium-saturation value and
#'     the measured DO concentration at that time
#'     \item \code{DO_obs_filter}: applicable only to \code{type='night'}: a
#'     smoothing filter is applied over the measured DO.obs values before
#'     applying nighttime regression
#'     \item \code{NA}: applicable only to \code{type='Kmodel'}, for which DO deficit is not estimated
#'   }
#' @seealso \code{\link{request_data}} for acquiring StreamPULSE data;
#'   \code{\link{prep_metabolism}} for organizing data and acquiring additional
#'   variables.
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
#' ## fit model with default parameters
#' modelfit = fit_metabolism(fitdata)
#'
#' ## fit model with custom parameters
#' class(fitdata) = 'data.frame' #just a formality, execute and disregard
#' modname = mm_name(type='bayes', pool_K600='binned',
#'     err_obs_iid=TRUE, err_proc_acor=FALSE, err_proc_iid=TRUE,
#'     ode_method = 'trapezoid', deficit_src='DO_mod', engine='stan')
#' modspecs = specs(modname)
#'
#' #overwrite default ln(Q) node centers
#' addis = tapply(log(fitdata$discharge), substr(fitdata$solar.time,1,10), mean)
#' modspecs$K600_lnQ_nodes_centers = seq(from=min(addis),
#'     to=max(addis), length.out=7)
#'
#' modelfit = metab(specs=modspecs, data=fitdata)
#' predictions = predict_metab(modelfit)
fit_metabolism = function(d, pool_K600='binned', err_obs_iid=TRUE,
    err_proc_acor=FALSE, err_proc_iid=TRUE, ode_method='trapezoid',
    deficit_src='DO_mod'){

    fitdata = d$data

    # check class of fitdata to determine which model to fit
    model = class(fitdata)
    if(model=="streamMetabolizer"){
        model_type = fitdata@type
    } else if(model=="BASE"){
        stop(paste("Modeling with BASE is not currently supported.\n\t",
            "Please rerun prep_metabolism with model='streamMetabolizer'."),
            call.=FALSE)
    }
    # then, reset class of fitdata to data.frame, may not be necessary?
    class(fitdata) = "data.frame"

    if(model=="streamMetabolizer"){

        #choose appropriate model specifications based on model type
        if(model_type == 'bayes'){
            # if(pool_K600 == 'none'){
            #     message(paste0("pool_K600 can't be 'none' in Bayesian ",
            #         "framework.\n\tSetting it to 'binned'."))
            #     pool_K600 = 'binned'
            # }
            engine = 'stan'#; pool_K600 = 'binned'; proc_err = TRUE
        } else {
            if(pool_K600 == 'binned'){
                message(paste0("pool_K600 can't be 'binned' in MLE ",
                    "framework.\n\tSetting it to 'none'."))
                pool_K600 = 'none'
            }
            if(proc_err == TRUE){
                message(paste0("proc_err can't be TRUE in MLE ",
                    "framework.\n\tSetting it to FALSE."))
                proc_err = FALSE
            }
            engine = 'nlm'#; pool_K600 = 'none'; proc_err = FALSE
        }

        # set most model specs
        modname = streamMetabolizer::mm_name(type=model_type, pool_K600=pool_K600,
            err_obs_iid=err_obs_iid, err_proc_acor=err_proc_acor,
            err_proc_iid=err_proc_iid,
            ode_method=ode_method, deficit_src=deficit_src, engine=engine)
        # modname = mm_name(type=model_type, pool_K600=pool_K600,
        #     err_obs_iid=TRUE, err_proc_acor=FALSE, err_proc_iid=proc_err,
        #     ode_method = "trapezoid", deficit_src="DO_mod", engine=engine)
        modspecs = streamMetabolizer::specs(modname)

        # fitdata2 <<- fitdata
        # stop('a')
        # fitdata = fitdata2
        #get average log daily discharge and use it to parameterize k v. Q curve
        if(engine == 'stan' && pool_K600 == 'binned'){
            addis = tapply(log(fitdata$discharge),
                substr(fitdata$solar.time,1,10), mean)
            # sum(is.infinite(log(fitdata$discharge))
            modspecs$K600_lnQ_nodes_centers = seq(from=min(addis, na.rm=TRUE),
                to=max(addis, na.rm=TRUE), length.out=7)
        }

        #fit model
        model_fit = streamMetabolizer::metab(specs=modspecs, data=fitdata)
        # return(model_fit)

    }else if(model=="BASE"){
        tmp = tempdir() # the temp dir for the data and model
        if(!dir.exists(tmp)) dir.create(tmp) # create if does not exist
        # create BASE directory
        directory = file.path(tmp,"BASE")
        if(!dir.exists(directory)){
            dir.create(directory)
            # add input folder
            dir.create(file.path(directory,"input"))
            # add output folder
            dir.create(file.path(directory,"output"))
            # - add instantaneous rates folder
            dir.create(file.path(directory,"output","instantaneous rates"))
            # - add validation plots folder
            dir.create(file.path(directory,"output","validation plots"))
        }
        # download BASE_metab_model_v2.2.txt
        download.file("https://raw.githubusercontent.com/streampulse/BASE/master/BASE/BASE_metab_model_v2.2.txt",
            file.path(directory,"BASE_metab_model_v2.2.txt"), quiet=TRUE)
        file.remove(list.files(file.path(directory,"input"),full.names=T)) # clear out input folder
        fitdata = split(fitdata, fitdata$Date)
        # write individual date base files
        lapply(fitdata, function(xx){
            if(nrow(xx)==96 && all(complete.cases(xx))){ # only full days
                date = unique(xx$Date)[1]
                write.csv(xx, file=paste0(directory,"/input/BASE_",date,".csv"), row.names=F)
            }
        })
        cat("Estimating metabolism with BASE.\n")
        # found in BASE_functions.R
        fit_BASE(directory=directory, interval=900, n.iter=30000, n.burnin=15000)
        # structure(list(output_directory = directory), class="BASE") # return the BASE directory with class BASE
        model_fit = list(output_directory=directory)
        attr(model_fit, 'class') = 'BASE'
    }

    #extract predictions from fit object (this block was formerly a
    # separate predict_metabolism function)
    if(class(model_fit)=="BASE"){
        directory = model_fit$output_directory
        preds = read.csv(paste0(directory,"/output/BASE_results.csv")) %>%
            tidyr::separate(File, c("fileX", "date", "extX"), "_|\\.") %>%
            select(-fileX, -extX) %>% mutate(date=as.Date(date))
        #develop stuff here if we ever use BASE again
    }else{
        predictions = streamMetabolizer::predict_metab(model_fit)
        output = list(predictions=predictions, fit=model_fit)
    }

    # ppp <<- predictions
    # fff <<- model_fit

    #extract data related to current model run and performance
    deets = extract_model_details(model_fit, predictions, d$specs)
    output$details = deets

    #if model covers <= one calendar year, it's eligible for inclusion
    #on the data portal; otherwise just return results to user
    mod_startyr = substr(d$specs$startdate, 1, 4)
    mod_endyr = substr(d$specs$enddate, 1, 4)
    if(mod_startyr != mod_endyr){
        cat(paste('Done. Returning your model fit and predictions.\n'))
        output$details$current_best = NULL
        return(output)
    }

    #assemble retrieval API call for the current best model's details
    cat(paste0('Checking StreamPULSE database for model results to\n\t',
        'compare with the model you just fit.\n'))

    u = paste0("https://data.streampulse.org/api/model_details_download?region=",
    # u = paste0("localhost:5000/api/model_details_download?region=",
        d$specs$region, "&site=", d$specs$site, "&year=", mod_startyr)

    #retrieve details for the current best model
    if(d$specs$token == 'none'){
        # r = httr::GET(u)
        tryCatch(R.utils::withTimeout(r <- httr::GET(u),
            timeout=12, onTimeout='error'),
            error=function(e){
                warning(paste('Could not reach StreamPULSE.',
                    'Check your internet connection.'), call.=FALSE)
                })
    } else {
        # r = httr::GET(u, httr::add_headers(Token=d$specs$token))
        tryCatch(R.utils::withTimeout(r <- httr::GET(u,
            httr::add_headers(Token=d$specs$token)), timeout=12,
            onTimeout='error'),
            error=function(e){
                warning(paste('Could not reach StreamPULSE.',
                    'Check your internet connection.'), call.=FALSE)
                })
    }

    if(exists('r')){
        json = httr::content(r, as="text", encoding="UTF-8")
        modspec = try(jsonlite::fromJSON(json), silent=TRUE)
    }

    #check for errors
    if(!exists('modspec') || class(modspec) == 'try-error'){
        cat(paste0('Failed to retrieve data from StreamPULSE.\n\t',
            'Returning your model fit and predictions.\n'))
        output$details$current_best = NULL
        return(output)
    }
    # if(length(modspec) == 1 && ! is.null(modspec$error)){
    #     return(output) #this line should never run; just here for later
    # }

    #if no model data on server, push this model up and return
    if(length(modspec$specs) == 0){
        cat(paste0("No model fits detected for this site and calendar year",
            ",\n\tso yours is the best fit by default!\n\t",
            "Pushing your results to the StreamPULSE database\n\t",
            "and returning model fit and predictions.\n"))
        push_model_to_server(output=output, deets=deets)
        output$details$current_best = NULL
        return(output)
    }

    #compare this model to the current best
    this_model_pen = get_model_penalty(deets)
    best_model_pen = get_model_penalty(modspec$specs)

    mods_equal = this_model_pen == best_model_pen
    pen_dif = best_model_pen - this_model_pen
    coverage_dif = deets$coverage - modspec$specs$coverage

    #if mod penalties equal, compare coverage. if better, push. return.
    if(mods_equal){
        if(coverage_dif > 0){ #this model has better coverage
            cat(paste0("Your model outperformed the best one on file!\n\t",
                "Pushing your results to the StreamPULSE database\n\t",
                "and returning your model fit and predictions.\n"))
            push_model_to_server(output=output, deets=deets)
        }

        output$details$current_best = NULL
        return(output)
    }

    #COVERAGE SHOULD NOT BE BASED SOLELY ON START AND END DATE, BUT ALSO NA!

    #if penalties differ, evaluate penaty and coverage differences
    accept = compare_models(pen_dif, coverage_dif)

    if(accept){
        cat(paste0("Your model outperformed the best one on file!\n\t",
            "Pushing your results to the StreamPULSE database\n\t",
            "and returning your model fit and predictions.\n"))
        push_model_to_server(output=output, deets=deets)
        output$details$current_best = NULL
        return(output)
    }

    cat(paste('Done. Returning your model fit and predictions.\n'))
    output$details$current_best = NULL
    return(output)

}
