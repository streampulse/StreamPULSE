#' Visualize metabolism model results
#'
#' Generates a Shiny R visualization similar to those found on the
#' \href{http://data.streampulse.org:3838/streampulse_diagnostic_plots/}
#' {StreamPULSE data portal}.
#'
#' @author Mike Vlah, \email{vlahm13@gmail.com}
#' @param model_out the output of \link{fit_metabolism}.
#' @return Opens a shiny app in default browser via localhost.
#' @seealso \code{\link{request_data}} for acquiring StreamPULSE data;
#'   \code{\link{prep_metabolism}} for organizing data and acquiring additional
#'   variables, and \code{\link{fit_metabolism}} for fitting models to data and
#'   predicting metabolism.
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
#' plot_output(modelfit)
plot_output = function(model_out){

    options(shiny.usecairo=TRUE)

    #create data objects
    modOut = list(data=model_out$fit@data,
        data_daily=model_out$fit@data_daily, fit=model_out$fit@fit)
    predictions = model_out$predictions
    fitpred = list('mod_out'=modOut, 'predictions'=predictions)

    #create mapping of input data fields and their pretty equivalents
    varmap = list('DO.sat'=list('DO sat', 'DO sat (%)'),
        'depth'=list('Depth', 'Depth (m)'),
        'temp.water'=list('Water temp',
            expression(paste('Water temp (', degree, 'C)'))),
        'light'=list('PAR', 'Light (PAR)'),
        'discharge'=list('Discharge',
            expression(paste('Discharge (m'^3, 's'^-1, ')'))))

    #populate overlay selector
    vars_etc = colnames(fitpred$mod_out$data)
    varinds = ! vars_etc %in% c('date', 'solar.time', 'DO.obs', 'DO.mod')
    vars = vars_etc[varinds]

    select_vars = vector('character', length=length(vars))
    for(i in 1:length(vars)){
        if(vars[i] %in% names(varmap)){
            select_vars[i] = varmap[[vars[i]]][[1]]
        }
    }

    #initiate time slider:
    #convert POSIX time to DOY and UNIX time
    DOY = as.numeric(gsub('^0+', '',
        strftime(fitpred$mod_out$data$solar.time, format="%j")))

    #get DOY bounds for slider
    DOYmin = ifelse(DOY[1] %in% 365:366, 1, DOY[1])
    DOYmax = DOY[length(DOY)]

    ui = fluidPage(

        #screen shouldn't go gray when plots are updating.
        tags$style(type="text/css", ".recalculating { opacity: 1.0; }" ),

        #hide error messages on screen
        tags$style(type="text/css",
            ".shiny-output-error { visibility: hidden; }",
            ".shiny-output-error:before { visibility: hidden; }"
        ),

        navbarPage(title=p(strong(a('StreamPULSE',
            href='https://data.streampulse.org/'))), inverse=TRUE,
            windowTitle='StreamPULSE Diagnostics',
            tabPanel('Model Performance',
                sidebarLayout(
                    sidebarPanel(
                        p(strong('Select DOY range:')),
                        p('Drag blue bar to move fixed range',
                            style=paste0('color:gray; font-size:80%;',
                                'padding:0; margin:0')),
                        p('Press play to autoscroll*',
                            style='color:gray; font-size:80%'),
                        htmlOutput('MPtime_slider'),
                        p('Click any point to view its date.',
                            style=paste0('color:gray; font-size:80%;')),
                        p(paste('*Residuals based on linear relationship',
                            'between daily mean K600 and log daily mean Q.'),
                            style=paste0('color:gray; font-size:80%;')),
                        width = 3
                    ),
                    mainPanel(
                        fluidRow(
                            column(6, align='left',
                                plotOutput('KvER', height='300px', width='auto',
                                    click='KvER_click'),
                                plotOutput('KvGPP', height='300px', width='auto',
                                    click='KvGPP_click')
                            ),
                            column(6, align='left',
                                plotOutput('KvQ', height='300px', width='auto',
                                    click='KvQ_click'),
                                plotOutput('QvKres', height='300px', width='auto',
                                    click='QvKres_click')
                            )
                        )
                    )
                )
            ),
            tabPanel(HTML('O<sub>2</sub> and Metabolism'),
                fluidRow(
                    column(12, align='left',
                        div(align='center', style=paste0(
                            'display: inline-block;',
                            'vertical-align:middle;',
                            'margin-right:1em'),
                            div(align='center', style=paste0(
                                'display: inline-block;',
                                'vertical-align:middle;',
                                'margin-right:1em'),
                                p(strong('Select DOY range:')),
                                p('Drag blue bar to move fixed range',
                                    style=paste0(
                                        'color:gray; font-size:80%;',
                                        'padding:0; margin:0')),
                                p('Press play to autoscroll*',
                                    style='color:gray; font-size:80%')
                            ),
                            div(align='left', style=paste0(
                                'display: inline-block;',
                                'vertical-align:middle;'),
                                htmlOutput('time_slider')
                            )
                        ),
                        hr()
                    )
                ),
                fluidRow(
                    column(9, align='center',
                        plotOutput('metab_legend', height='20px',
                            width='auto'),
                        plotOutput('metab_plot', height='300px', width='auto'),
                        plotOutput('O2_legend', height='20px',
                            width='auto'),
                        plotOutput('O2_plot', brush='O2_brush',
                            height='300px', width='auto')

                    ),
                    column(3, align='center',
                        p(strong(HTML('Cumulative O<sub>2</sub> (gm<sup>-2</sup>d<sup>-1</sup>)'))),
                        tableOutput('cumul_metab'),
                        br(),
                        # plotOutput('cumul_legend', height='20px',
                        #     width='auto'),
                        # plotOutput('cumul_plot', height='200px', width='auto'),
                        plotOutput('kernel_legend', height='20px',
                            width='auto'),
                        plotOutput('kernel_plot', height='200px', width='auto'),
                        br(),
                        selectInput('metab_overlay', 'Model param overlay',
                            list('None', 'mean daily K600'), selected='None'),
                        selectInput('O2_overlay', 'Input data overlay',
                            as.list(c('None', select_vars)), selected='None'),
                        radioButtons('xformat', 'Series x-axis', inline=TRUE,
                            list('DOY', 'Date'), selected='DOY')
                    )
                ),
                br(),
                p(paste("*At this time, plotting may fail if your data",
                    "do not all occur within the same calendar year.",
                    "If plots do not load, try toggling between tabs."),
                    style='color:gray; font-size:100%')
            )
        )
    )

    server = function(input, output, session){

        output$time_slider = renderUI({

            if(!is.null(fitpred)){

                #convert POSIX time to DOY and UNIX time
                DOY = as.numeric(gsub('^0+', '',
                    strftime(fitpred$mod_out$data$solar.time, format="%j")))

                #get DOY bounds for slider
                DOYmin = ifelse(DOY[1] %in% 365:366, 1, DOY[1])
                DOYmax = DOY[length(DOY)]

                sliderInput("range", label=NULL,
                    min=DOYmin, max=DOYmax, value=c(DOYmin, DOYmax),
                    ticks=TRUE, step=3,
                    animate=animationOptions(interval=2000)
                )
            }
        })

        output$MPtime_slider = renderUI({

            if(!is.null(fitpred)){

                #convert POSIX time to DOY and UNIX time
                MPDOY = as.numeric(gsub('^0+', '',
                    strftime(fitpred$mod_out$data$solar.time, format="%j")))

                #get DOY bounds for slider
                MPDOYmin = ifelse(MPDOY[1] %in% 365:366, 1, MPDOY[1])
                MPDOYmax = MPDOY[length(MPDOY)]

                sliderInput("MPrange", label=NULL,
                    min=MPDOYmin, max=MPDOYmax, value=c(MPDOYmin, MPDOYmax),
                    ticks=TRUE, step=3,
                    animate=animationOptions(interval=2000)
                )
            }
        })

        #update model performance data frames based on time range selection
        get_slices = eventReactive({
            input$MPrange
        }, {

            if(! is.null(fitpred)){
                mod_out = fitpred$mod_out

                MPstart = input$MPrange[1]
                MPend = input$MPrange[2]

                #convert POSIX time to DOY and UNIX time
                DOY = as.numeric(gsub('^0+', '', strftime(mod_out$data$solar.time,
                    format="%j")))
                date = as.Date(gsub('^0+', '', strftime(mod_out$data$solar.time,
                    format="%Y-%m-%d")))

                # replace initial DOYs of 365 or 366 (solar date in previous calendar year) with 1
                if(DOY[1] %in% 365:366){
                    DOY[DOY %in% 365:366 & 1:length(DOY) < length(DOY)/2] = 1
                }

                #filter data by date bounds specified in time slider
                xmin_ind = match(MPstart, DOY)
                if(is.na(xmin_ind)) xmin_ind = 1
                xmin = date[xmin_ind]

                xmax_ind = length(DOY) - match(MPend, rev(DOY)) + 1
                if(is.na(xmax_ind)) xmax_ind = nrow(mod_out$data)
                xmax = date[xmax_ind]

                daily_slice = mod_out$fit$daily[mod_out$fit$daily$date <= xmax &
                        mod_out$fit$daily$date >= xmin,]
                data_daily_slice = mod_out$data_daily[mod_out$data_daily$date <= xmax &
                        mod_out$data_daily$date >= xmin,]

                out = list(daily_slice=daily_slice, data_daily_slice=data_daily_slice,
                    mod_out=mod_out)
            }
        })

        observeEvent({
            input$range
        }, {
            start = input$range[1]
            end = input$range[2]

            if(!is.null(start) && !is.null(end)){
                output$metab_legend = renderPlot({
                    if(input$metab_overlay != 'None'){
                        metab_legend(show_K600=TRUE)
                    } else {
                        metab_legend(show_K600=FALSE)
                    }
                })

                output$metab_plot = renderPlot({
                    ts_full = processing_func(fitpred$predictions, st=start,
                        en=end)
                    par(mar=c(1,4,0,4), oma=rep(0,4))
                    daily = fitpred$mod_out$fit$daily
                    daily$doy = as.numeric(gsub('^0+', '',
                        strftime(daily$date, format="%j")))
                    daily = daily[daily$doy > start & daily$doy < end,]
                    season_ts_func(ts_full, daily, st=start, en=end,
                        input$metab_overlay)
                })

                output$kernel_legend = renderPlot({
                    kernel_legend()
                })

                output$kernel_plot = renderPlot({
                    ts_full = processing_func(fitpred$predictions, st=start,
                        en=end)
                    par(mar=c(3,3.5,0,.5), oma=rep(0,4))
                    kernel_func(ts_full, 'Name and Year')
                })

                output$O2_legend = renderPlot({
                    O2_legend(overlay=input$O2_overlay, varmap=varmap)
                })

                output$O2_plot = renderPlot({
                    par(mar=c(3,4,0,4), oma=rep(0,4))
                    O2_plot(mod_out=fitpred$mod_out, st=start, en=end,
                        brush=input$O2_brush, overlay=input$O2_overlay,
                        xformat=input$xformat, varmap=varmap)
                })

                output$cumul_metab = renderTable({
                    ts_full = processing_func(fitpred$predictions, st=start,
                        en=end)
                    na_rm = na.omit(ts_full)
                    gppsum = sum(na_rm$GPP, na.rm=TRUE)
                    ersum = sum(na_rm$ER, na.rm=TRUE)
                    nepsum = sum(na_rm$NPP, na.rm=TRUE)
                    return(data.frame('GPP'=gppsum, 'ER'=ersum, 'NEP'=nepsum))
                }, striped=TRUE)
            }
        })

        observeEvent({
            input$MPrange
        }, {

            slices = get_slices()
            mod_out = slices$mod_out
            MPstart = input$MPrange[1]
            MPend = input$MPrange[2]

            if(!is.null(MPstart) && !is.null(MPend)){

                output$KvER = renderPlot({
                    if(!is.null(fitpred$mod_out)){
                        KvER_plot(mod_out=fitpred$mod_out,
                            slice=slices$daily_slice)
                    }
                })

                output$KvQ = renderPlot({
                    if(!is.null(fitpred$mod_out)){
                        KvQ_plot(mod_out=fitpred$mod_out,
                            slicex=slices$data_daily_slice,
                            slicey=slices$daily_slice)
                    }
                })

                output$KvGPP = renderPlot({
                    if(!is.null(fitpred$mod_out)){
                        KvGPP_plot(mod_out=fitpred$mod_out,
                            slice=slices$daily_slice)
                    }
                })

                output$QvKres = renderPlot({
                    if(!is.null(fitpred$mod_out)){
                        QvKres_plot(mod_out=fitpred$mod_out,
                            slicex=slices$data_daily_slice,
                            slicey=slices$daily_slice)
                    }
                })
            }

        })

        #replot when a click is registered; don't react when click handler flushes
        observeEvent({
            if (! is.null(input$KvER_click$x) ||
                    ! is.null(input$KvQ_click$x) ||
                    ! is.null(input$QvKres_click$x) ||
                    ! is.null(input$KvGPP_click$x)) TRUE
            else NULL
        }, {

            slices = get_slices()
            mod_out = slices$mod_out
            MPstart = input$MPrange[1]
            MPend = input$MPrange[2]

            output$KvER = renderPlot({
                KvER_plot(mod_out=mod_out,
                    slice=slices$daily_slice, click=isolate(input$KvER_click))
            })

            output$KvQ = renderPlot({
                KvQ_plot(mod_out=mod_out, slicex=slices$data_daily_slice,
                    slicey=slices$daily_slice, click=isolate(input$KvQ_click))
            })

            output$KvGPP = renderPlot({
                KvGPP_plot(mod_out=mod_out,
                    slice=slices$daily_slice, click=isolate(input$KvGPP_click))
            })

            output$QvKres = renderPlot({
                if(!is.null(mod_out)){
                    QvKres_plot(mod_out=mod_out, slicex=slices$data_daily_slice,
                        slicey=slices$daily_slice, click=isolate(input$QvKres_click))
                }
            })

        }, ignoreNULL=TRUE)

    }

    shinyApp(ui=ui, server=server, options=list(launch.browser=TRUE))
}
