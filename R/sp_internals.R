#' Internal functions
#'
#' Not intended to be called directly by the user.
#'
#' Not intended to be called directly by the user.
#'
#' @keywords internal
#'
FindandCollect_airpres = function(lat, long, start_datetime, end_datetime) {
    tf = tempfile()
    download.file("ftp://ftp.ncdc.noaa.gov/pub/data/noaa/isd-history.txt",tf,mode="wb")
    noaa.sites <- read.fwf(tf, skip = 22, header = F, widths = c(7,6,30, 5, 3, 6, 8, 9, 8, 9, 8), comment.char = "", col.names = c("USAF", "WBAN", "STATION NAME", "CTRY", "ST", "CALL", "LAT", "LON", "ELEV(M)", "BEGIN", "END"), flush = TRUE)
    noaa.sites <- na.omit(noaa.sites)
    noaa.sites <- noaa.sites %>%
        mutate(LAT = as.numeric(as.character(LAT))) %>%
        mutate(LON = as.numeric(as.character(LON))) %>%
        filter(LAT < (lat + 5) & LAT > (lat - 5) & LON < (long + 5) & LON > (long - 5))
    pt1 <- cbind(rep(long, length.out = length(noaa.sites$LAT)), rep(lat, length.out = length(noaa.sites$LAT)))
    pt2 <- cbind(noaa.sites$LON, noaa.sites$LAT)
    dist <- diag(distm(pt1, pt2, fun = distHaversine))/1000
    noaa.sites$dist <- dist
    tmp <- which((as.numeric(substr(noaa.sites$END,1,4)) >= as.numeric(substr(end_datetime, 1, 4))) & as.numeric(substr(noaa.sites$BEGIN,1,4)) <= as.numeric(substr(start_datetime, 1, 4)))
    noaa.sites <- noaa.sites[tmp,]
    noaa.sites <- noaa.sites[with(noaa.sites, order(dist)),]

    yrs <- seq(as.numeric(substr(start_datetime, 1, 4)),as.numeric(substr(end_datetime, 1, 4)), by = 1)
    for (i in 1:length(noaa.sites$dist)) {
        k <- i
        available <- vector(mode = 'logical', length = length(yrs))
        USAF <- as.character(noaa.sites$USAF[i])
        if(nchar(as.character(noaa.sites$WBAN[i])) == 5){
            WBAN <- as.character(noaa.sites$WBAN[i])
        } else {
            WBAN <- paste0(0,as.character(noaa.sites$WBAN[i]))
        }
        y <- as.data.frame(matrix(NA, nrow = 1, ncol = 12))
        for(j in 1:length(yrs)){
            tf = tempfile()
            download.file(paste0("ftp://ftp.ncdc.noaa.gov/pub/data/noaa/isd-lite/",
                yrs[j],"/",USAF,"-",WBAN,"-",yrs[j],".gz"),tf,mode="wb")
            x = read.table(tf)
            x[x==-9999] = NA
            if(length(which(!is.na(x$V7))) >= 0.9 * length(x$V7)) {
                available[j] <- TRUE
                y <- rbind(x,y)
            }else {
                break
            }
        }
        if(length(yrs) == length(which(available))){
            break
        }
    }
    y <- y[!is.na(y$V1),]
    colnames(y) = c("y","m","d","h","air_temp","dewtemp","air_kPa","winddir","sindspeed","skycover","precip1h","precip6h")
    y$air_kPa = y$air_kPa/100
    y$air_temp = y$air_temp/10
    y$DateTime_UTC = parse_datetime(paste0(y$y,"-",sprintf("%02d",y$m),"-",sprintf("%02d",y$d)," ",sprintf("%02d",y$h),":00:00 0"), "%F %T %Z")
    y <- y[with(y, order(DateTime_UTC)),]
    y = as_tibble(y) %>% select(DateTime_UTC,air_temp,air_kPa)
    ss = tibble(DateTime_UTC=seq(y$DateTime_UTC[1], y$DateTime_UTC[nrow(y)], by=900))
    xx = left_join(ss, y, by = "DateTime_UTC")
    xx = mutate(xx, air_temp=na.approx(air_temp), air_kPa=na.approx(air_kPa))
    daterng = c(start_datetime, end_datetime)
    xtmp = xx %>% filter(DateTime_UTC>=daterng[1] & DateTime_UTC<=daterng[2])
    # select(xtmp, DateTime_UTC, air_kPa, air_temp)
    # print(noaa.sites[k,])
    return(select(xtmp, DateTime_UTC, air_kPa, air_temp))
}

retrieve_air_pressure = function(md, dd){

    # md2 <<- md
    # dd2 <<- dd
    # stop('a')
    # md = md2
    # dd = dd2

    lat = md$lat
    long = md$lon
    start_datetime = dd$DateTime_UTC[1]
    end_datetime = dd$DateTime_UTC[nrow(dd)]

    cat('Collecting air pressure data from NCDC (NOAA).\n')
    capture.output(df <- suppressWarnings(FindandCollect_airpres(lat, long,
        start_datetime, end_datetime)))
    cat('\n')

    df_out = df %>% mutate(AirPres_kPa = air_kPa) %>%
        select(DateTime_UTC, AirPres_kPa) %>% as.data.frame()

    return(df_out)
}

retrieve_air_pressure2 = function(md, dd){

    # sites2 <<- sites
    # dd2 <<- dd
    sites = md

    #format site data for use with geoknife package
    station = as.data.frame(t(sites[,c('lon','lat')]))
    station = simplegeom(station)

    years = unique(substr(dd$DateTime_UTC, 1, 4))
    cat('Acquiring air pressure',
        'data for', length(years),
        'year(s).\n\tEach year may take a few minutes.\n')

    #retrieve air pressure data from noaa
    pres = data.frame(datetime=.POSIXct(character()), pres=numeric())
    for(i in 1:length(years)){

        fabric = webdata(url=paste0('https://www.esrl.noaa.gov/psd/th',
            'redds/dodsC/Datasets/ncep.reanalysis/surface/pres.sfc.',
            years[i], '.nc'), variables='pres')
        noaa_job = geoknife(stencil=station, fabric=fabric, wait=TRUE)
        noaa_data = result(noaa_job, with.units=TRUE)

        pres = rbind(pres, noaa_data[,c('DateTime','1')])

        cat('Year', i, 'complete.\n')
    }

    pres = data.frame(pres)

    df_out = pres %>% mutate(AirPres_kPa = X1 / 1000,
        DateTime_UTC=DateTime) %>% select(AirPres_kPa, DateTime_UTC)

    return(df_out)
}

estimate_discharge = function(Z=NULL, Q=NULL, a=NULL, b=NULL,
    sh=NULL, dd=NULL, fit, ignore_oob_Z, plot){

    # dd2 <<- dd
    # stop('a')
    # dd = dd2
    # if(is.numeric(sh)){ #then need to calculate depth. based on:
    #https://web.archive.org/web/20170617070623/http://www.onsetcomp.com/files/support/tech-notes/onsetBCAguide.pdf

    if(is.null(plot)) plot = TRUE
    if(is.null(fit)) fit = 'power'
    if(plot){
        defpar = par()
        if(is.null(Z)){
            par(mfrow=c(1,1), mar=c(4,4,1,1), oma=c(0,0,0,0))
        } else {
            par(mfrow=c(2,1), mar=c(4,4,1,1), oma=c(0,0,0,0))
        }
    }

    if(! 'Depth_m' %in% colnames(dd)){

        cat(paste0('No depth data detected. Estimating level (AKA stage) from ',
            'water pressure.\n'))
        if(any(! c('WaterPres_kPa','AirPres_kPa') %in% colnames(dd))){
            stop(paste0('Air and/or water pressure not detected.',
                '\n\tNot enough information to proceed.'),
                call.=FALSE)
        }

        #remove air pressure so all pressure is from hydraulic head
        hyd_pres = dd$WaterPres_kPa - dd$AirPres_kPa

        if(! 'WaterTemp_C' %in% colnames(dd)){
            warning(paste0('Water temperature not detected.',
                '\n\tAssuming density of water is 1 g/cm^3.'),
                call.=FALSE)

            fl_dens = 1000 #g/m^3

        } else {

            #compute fluid density
            wat_temp = dd$WaterTemp_C
            T1 = 16.945176 * wat_temp
            T1b = 16.879850e-03 * wat_temp
            T2 = 7.9870401e-03 * wat_temp^2
            T3 = 46.170461e-06 * wat_temp^3
            T4 = 105.56302e-09 * wat_temp^4
            T5 = 280.54253e-12 * wat_temp^5
            fl_dens = (999.83952 + T1 - T2 - T3 + T4 - T5) / (1 + T1b)
        }

        #convert fluid density to lb/ft^3
        fl_dens = 0.0624279606 * fl_dens

        #convert hydraulic pressure to fluid depth
        ft_to_m = 0.3048
        kPa_to_psi = 0.1450377
        psi_to_psf = 144.0
        fl_depth = ft_to_m * (kPa_to_psi * psi_to_psf * hyd_pres) / fl_dens

    } else { #else we have depth already
        fl_depth = dd$Depth_m
    }

    #correct for sensor height above bed if sensor_height supplied
    if(!is.null(sh)){

        depth = fl_depth + sh #sh needs to be in meters

        message(paste0('Computing depth as level (AKA stage) plus ',
            'sensor height.\n\tMake sure other parameters supplied to ',
            'zq_curve are also based on depth,\n\trather than level,',
            ' or else omit sensor_height argument.'))

        dep_or_lvl = 'Depth'

        # cat('Quantiles of computed depth (m):\n')

    } else {

        depth = fl_depth

        message(paste0('Without sensor_height argument, ZQ rating curve will ',
            'be based on level (AKA stage),\n\trather than depth. Make sure ',
            'other parameters supplied to zq_curve are\n\talso based on level,',
            ' or else include sensor_height argument.'))

        dep_or_lvl = 'Level'

        # cat('Quantiles of computed level (m):\n')
    }
    # print(round(quantile(depth, na.rm=TRUE), 4))

    #generate rating curve if Z and Q sample data supplied
    if(!is.null(Z)){

        Q = Q[order(Z)]
        Z = Z[order(Z)] #just making sure there's no funny business

        # Q2 <<- Q
        # Z2 <<- Z
        # a2 <<- a
        # b2 <<- b
        # stop()
        # Q = Q2
        # Z = Z2
        # a = a2
        # b = b2

        #try to fit power model
        if(fit == 'power'){
            # mod = tryCatch(nls(Q ~ (a * Z^b), start=list(a=0.1, b=1)),
            #     error=function(e){
            #         stop(paste0('Failed to fit rating curve.\n\tThis is worth ',
            #             'mentioning to Mike: vlahm13@gmail.com.\n\t',
            #             'Note that you can fit your own curve and then supply\n\t',
            #             'a and b (of Q=aZ^b) directly.'), call.=FALSE)
            #     })
            mod = try(nls(Q ~ (a * Z^b), start=list(a=0.1, b=1)), silent=TRUE)
            if(class(mod) == 'try-error'){
                mod = try(nls(Q ~ (a * Z^b), start=list(a=1, b=1)), silent=TRUE)
                if(class(mod) == 'try-error'){
                    stop(paste0('Failed to fit rating curve.\n\tThis is worth ',
                        'mentioning to Mike: vlahm13@gmail.com.\n\t',
                        'Note that you can fit your own curve and then supply',
                        '\n\ta and b (of Q=aZ^b) directly.'), call.=FALSE)
                }
            }
        } else { #try to fit exponential
            if(fit == 'exponential'){
                mod = tryCatch(nls(Q ~ (a * exp(b * Z)),
                    start=list(a=0.01, b=1)),
                    error=function(e){
                        stop(paste0('Failed to fit rating curve. Try using ',
                            'fit="power" instead.\n\tAlso ',
                            'note that you can fit your own curve and then ',
                            'supply\n\ta and b (of Q=ae^(bZ)) directly.'),
                            call.=FALSE)
                    })
            } else { #fit linear
                mod = nls(Q ~ a * Z + b, start=list(a=1, b=0))
            }
        }

        if(plot == TRUE){
            plot(Z, Q, xlab='', ylab='Q samp. data (cms)',
                las=1, main=paste0('Rating curve fit (', fit, ')'))
            mtext('Z sample data (m)', 1, 2)
            plotseq = round(seq(Z[1], Z[length(Z)], diff(range(Z))/50), 2)
            lines(plotseq, predict(mod, list(Z=plotseq)))
        }

        params = summary(mod)$parameters
        a = params[1,1]
        b = params[2,1]

        #display info about curve
        eqns = list('power'='Q=aZ^b', 'exponential'='Q=ae^(bZ)',
            'linear'='Q=aZ+b')
        cat(paste0('Rating curve summary\n\tFit: ', fit, '\n\tEquation: ',
            eqns[fit], '\n\tParameters: a=', round(a, 3), ', b=',
            round(b, 3), '\n'))
        maxZ = round(max(Z, na.rm=TRUE), 2)
        maxD = round(max(depth, na.rm=TRUE), 2)
        prop_Z_oob = round(sum(depth > maxZ, na.rm=TRUE) / length(depth), 1)

        if(ignore_oob_Z){
            depth[depth > maxZ] = NA
            message(paste0(prop_Z_oob * 100, '% of sensor ',
                tolower(dep_or_lvl), ' values ',
                'exceeded Z in the rating curve.\n\tThese have been replaced',
                'with NA because ignore_oob_Z is set to TRUE.'))
        } else {
            if(maxD > maxZ){
                warning(paste0('Max observed ', tolower(dep_or_lvl), ' = ',
                    maxD, '. Max observed input Z = ', maxZ, '.\n\tDischarge ',
                    'estimates for ', tolower(dep_or_lvl), ' > ', maxZ,
                    ' (', prop_Z_oob * 100, '% of observations)',
                    ' may be untrustworthy if using power or exponential ',
                    'fit.\n\tSet ',
                    'ignore_oob_Z=TRUE if this bothers you.'), call.=FALSE)
            }
        }
    } #else a and b have been supplied directly

    #estimate discharge using a and b params from rating curve
    if(fit == 'power'){
        discharge = a * depth^b
    } else {
        if(fit == 'exponential'){
            discharge = a * exp(b * depth)
        } else {
            discharge = a * depth + b
        }
    }

    if(plot){
        plot(depth, discharge, xlab='',
            xlim=c(max(min(depth, na.rm=TRUE), 0), max(depth, na.rm=TRUE)),
            ylab='Est. Q (cms)', main='Rating curve prediction', las=1)
        mtext(paste(dep_or_lvl, 'series data (m)'), 1, 2)
    }
    suppressWarnings(par(defpar))

    return(discharge)
}

create_sm_class = function(){

    streamMetabolizer = setClass("streamMetabolizer",
        representation(type="character"), contains="data.frame")

    return(streamMetabolizer)
}
streamMetabolizer = create_sm_class()
