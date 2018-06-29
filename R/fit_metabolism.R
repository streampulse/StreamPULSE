#' Fit stream metabolism models
#'
#' Passes data and model parameters to \code{streamMetabolizer} or
#' \code{BASE} according to arguments passed to \link{prep_metabolism}.
#' NOTE: support for modeling with
#' \code{BASE} is currently in development. Please use \code{streamMetabolizer}
#' in the meantime.
#'
#' This function supports only default streamMetabolizer parameters at the
#' moment. Namely, \code{err_obs_iid=TRUE, err_proc_acor=FALSE,}
#' \code{ode_method='trapezoid', deficit_src='DO_mod'}.
#' If \code{prep_metabolism}'s
#' \code{model_type} argument is set to \code{'mle'}, then the following
#' defaults will also be set: \code{engine='nlm', pool_K600='none',}
#' \code{proc_err=FALSE}. If \code{model_type} is set to \code{'bayes'},
#' these will
#' instead be set to \code{engine='stan', pool_K600='binned', proc_err=TRUE},
#' and \code{K600_lnQ_nodes_centers} will be a sequence of length 7 from
#' the minimum to the maximum natural log of daily discharge.
#'
#' Until this function is made more flexible, you can adjust the above
#' parameters by using \code{streamMetabolizer}'s functions directly.
#' See the example section below for more.
#'
#' @author Mike Vlah, \email{vlahm13@gmail.com}
#' @author Aaron Berdanier
#' @param d the output of \link{prep_metabolism}.
#' @return returns a list containing the output of \code{streamMetabolizer}'s
#'   \code{metab} function (fit) and of \code{streamMetabolizer}'s
#'   \code{predict_metab} function (predictions).
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
fit_metabolism = function(d){

    fitdata = d$data

    # check class of fitdata to determine which model to fit
    model = class(fitdata)
    if(model=="streamMetabolizer") model_type = fitdata@type
    # then, reset class of fitdata to data.frame, may not be necessary?
    class(fitdata) = "data.frame"

    if(model=="streamMetabolizer"){

        #choose appropriate model specifications based on model type
        if(model_type=='bayes'){
            engine = 'stan'; pool_K600='binned'; proc_err = TRUE
        } else {
            engine = 'nlm'; pool_K600='none'; proc_err = FALSE
        }

        # set most model specs
        modname = mm_name(type=model_type, pool_K600=pool_K600,
            err_obs_iid=TRUE, err_proc_acor=FALSE, err_proc_iid=proc_err,
            ode_method = "trapezoid", deficit_src="DO_mod", engine=engine)
        modspecs = specs(modname)

        # fitdata2 <<- fitdata
        # stop('a')
        # fitdata = fitdata2
        #get average log daily discharge and use it to parameterize k v. Q curve
        if(engine == 'stan'){
            addis = tapply(log(fitdata$discharge),
                substr(fitdata$solar.time,1,10), mean)
            # sum(is.infinite(log(fitdata$discharge))
            modspecs$K600_lnQ_nodes_centers = seq(from=min(addis),
                to=max(addis), length.out=7)
        }

        #fit model
        model_fit = metab(specs=modspecs, data=fitdata)
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
            separate(File, c("fileX", "date", "extX"), "_|\\.") %>%
            select(-fileX, -extX) %>% mutate(date=as.Date(date))
        #develop stuff here if we ever use BASE again
    }else{
        predictions = predict_metab(model_fit)
        output = list(predictions=predictions, fit=model_fit)
    }

    #if model covers <= one calendar year...
    # ppp <<- predictions
    # fff <<- model_fit
    mod_startyr = substr(d$specs$startdate, 1, 4)
    mod_endyr = substr(d$specs$enddate, 1, 4)
    if(mod_startyr == mod_endyr){

        #assemble model spec retrieval API call
        cat(paste0('Checking StreamPULSE database for model results to\n\t',
            'compare with the model you just fit.'))

        #LOCALHOST
        u = paste0("localhost:5000/api/model_details_download?region=",
            d$specs$region, "&site=", d$specs$site, "&year=", mod_startyr)
        # u = paste0("localhost:5000/api/model_details?region=NC",
        #     "&site=Eno&year=2019")

        #retrieve specs for the current best model
        if(d$specs$token == 'none'){
            r = httr::GET(u)
        }else{
            r = httr::GET(u, httr::add_headers(Token=d$specs$token))
        }
        json = httr::content(r, as="text", encoding="UTF-8")
        modspec = try(jsonlite::fromJSON(json), silent=TRUE)

        #check for errors
        if(class(modspec) == 'try-error'){
            cat(paste0('Failed to retrieve data from StreamPULSE.\n\t',
                'Returning your model fit and predictions.'))
            return(output)
        }
        if(length(modspec) == 1 && ! is.null(modspec$error)){
            return(output) #this line should never run; just here for later
        }

        #extract results related to model run and performance
        deets = extract_model_details(model_fit, predictions, d$specs,
            mod_startyr)
        # ddd <<- deets

        #if no model data on server, push this model up
        if(length(modespec$specs) == 0){
            cat(paste0("No model fits detected for this site and calendar year",
                ",\n\tso this is the best fit by default!\n\t",
                "Pushing your results to the StreamPULSE database.\n\t"))

            #first push model details to database
            u2 = paste0("localhost:5000/api/model_details_upload")
            o = httr::POST(url=u2, body=deets, encode='form') #send data

            jsono = httr::content(o, as="text", encoding="UTF-8") #get response
            oo = try(jsonlite::fromJSON(jsono), silent=TRUE)

            #check for errors
            if(class(oo) == 'try-error' || oo$callback != 'success'){
                cat(paste0('Failed to push data to StreamPULSE.\n\t',
                    'Returning your model fit and predictions.'))
                return(output)
            }

            #then save model output and predictions to temp files as RDS
            data_daily = model_fit@data_daily
            data = model_fit@data
            fit = model_fit@fit
            mod_out = list('data_daily'=data_daily, 'data'=data, 'fit'=fit)
            tmp1 = tempfile()
            saveRDS(mod_out, file=tmp1)
            tmp2 = tempfile()
            saveRDS(predictions, file=tmp2)

            #then push those RDS files to server
            file_id = paste(deets$region, deets$site, '2019', sep='_')
            # file_id = paste(deets$region, deets$site, deets$year, sep='_')
            u3 = paste0("localhost:5000/api/model_upload")
            o = httr::POST(url=u3,
                body=list(modOut=httr::upload_file(tmp1),
                    predictions=httr::upload_file(tmp2)),
                httr::add_headers(file_id=file_id))

            jsono = httr::content(o, as="text", encoding="UTF-8") #get response
            oo = try(jsonlite::fromJSON(jsono), silent=TRUE)

            if(class(oo) == 'try-error' || oo$callback != 'success'){
                cat(paste0('Failed to push data to StreamPULSE.\n\t',
                    'Returning your model fit and predictions.'))
                return(output)
            }

            cat('Done!')
            return(output)
        }

        # mmm <<- modspec




    }

}
