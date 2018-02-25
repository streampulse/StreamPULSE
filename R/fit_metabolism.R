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
#' @param fitdata the output of \link{prep_metabolism}.
#' @return returns the output of \code{streamMetabolizer}'s \code{metab}
#'   function.
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
fit_metabolism = function(fitdata){
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
        modfit = metab(specs=modspecs, data=fitdata)
        return(modfit)

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
        structure(list(output_directory = directory), class="BASE") # return the BASE directory with class BASE
    }
}
