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

    # library(shiny)
    # library(Cairo)
    # library(ks)
    # library(scales)

    options(shiny.usecairo=TRUE)

    #create mapping of input data fields and their pretty equivalents
    varmap = list('DO.sat'=list('DO sat', 'DO sat (%)'),
        'depth'=list('Depth', 'Depth (m)'),
        'temp.water'=list('Water temp',
            expression(paste('Water temp (', degree, 'C)'))),
        'light'=list('PAR', 'Light (PAR)'),
        'discharge'=list('Discharge',
            expression(paste('Discharge (m'^3, 's'^-1, ')'))))

    #create data objects
    modOut = list(data=model_out$fit@data,
        data_daily=model_out$fit@data_daily, fit=model_out$fit@fit)
    predictions = model_out$predictions
    fitpred = list('mod_out'=modOut, 'predictions'=predictions)

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
                                plotOutput('KvER', height='auto',
                                    click='KvER_click'),
                                plotOutput('KvGPP', height='auto',
                                    click='KvGPP_click')
                            ),
                            column(6, align='left',
                                plotOutput('KvQ', height='auto',
                                    click='KvQ_click'),
                                plotOutput('QvKres', height='auto',
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

                        plotOutput('metab_plot', height='200px', width='auto'),
                        plotOutput('O2_legend', height='20px',
                            width='auto'),
                        plotOutput('O2_plot', brush='O2_brush',
                            height='200px', width='auto')

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
                        plotOutput('kernel_plot', height='200px', width='auto')
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

        observeEvent({
            input$range
        }, {
            start = input$range[1]
            end = input$range[2]

            if(!is.null(start) && !is.null(end)){
                output$metab_legend = renderPlot({
                    metab_legend()
                })

                output$cumul_legend = renderPlot({
                    cumul_legend()
                })

                output$metab_plot = renderPlot({
                    ts_full = processing_func(fitpred$predictions, st=start,
                        en=end)
                    par(mar=c(1,4,0,1), oma=rep(0,4))
                    season_ts_func(ts_full, TRUE, st=start, en=end)
                })

                output$cumul_plot = renderPlot({
                    ts_full = processing_func(fitpred$predictions, st=start,
                        en=end)
                    par(mar=c(3,3.5,0.2,0.5), oma=rep(0,4))
                    cumulative_func(ts_full, st=start, en=end)
                })

                output$O2_legend = renderPlot({
                    O2_legend()
                })

                output$kernel_legend = renderPlot({
                    kernel_legend()
                })

                output$O2_plot = renderPlot({
                    par(mar=c(3,4,0,1), oma=rep(0,4))
                    O2_plot(mod_out=fitpred$mod_out, st=start, en=end,
                        input$O2_brush)
                })

                output$kernel_plot = renderPlot({
                    ts_full = processing_func(fitpred$predictions, st=start,
                        en=end)
                    par(mar=c(3,3.5,0,.5), oma=rep(0,4))
                    kernel_func(ts_full, 'Name and Year')
                })
            }
        })

        output$KvQvER = renderPlot({
            try(KvQvER_plot(mod_out=modOut), silent=TRUE)
        })

    }

    shinyApp(ui=ui, server=server, options=list(launch.browser=TRUE))
}
