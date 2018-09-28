#' Visualize metabolism model results
#'
#' Generates a Shiny R visualization similar to those found on the
#' \href{http://data.streampulse.org:3838/streampulse_diagnostic_plots/}
#' {StreamPULSE data portal}.
#'
#' @author Mike Vlah, \email{vlahm13@gmail.com}
#' @param model_out the output of \link{fit_metabolism}.
#' @return Opens a shiny app via localhost, using the default browser.
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

    library(shiny)
    library(Cairo)
    library(ks)
    library(scales)
    library(shinyjs)

    options(shiny.usecairo=TRUE)

    #create data objects
    modOut = list(data=model_out$fit@data,
        data_daily=model_out$fit@data_daily, fit=model_out$fit@fit)
    predictions = model_out$predictions
    fitpred = list('mod_out'=modOut, 'predictions'=predictions)

    # #initiate time slider:
    # #convert POSIX time to DOY and UNIX time
    # DOY = as.numeric(gsub('^0+', '',
    #     strftime(fitpred$mod_out$data$solar.time, format="%j")))
    #
    # #get DOY bounds for slider
    # DOYmin = ifelse(DOY[1] %in% 365:366, 1, DOY[1])
    # DOYmax = DOY[length(DOY)]

    get_plotheight = "
    shinyjs.init = function() {
    $(window).resize(shinyjs.getHeight50);
    }

    //shinyjs.calcHeight = function(propHeight) {
    //  var h = $(window).height() * propHeight;
    //  Shiny.onInputChange('plotHeight', Number(h.toFixed(0)));

    shinyjs.getHeight50 = function() {
    Shiny.onInputChange('height50', $(window).height() * .5);
    }
    shinyjs.getHeight40 = function() {
    Shiny.onInputChange('height40', $(window).height() * .4);
    }
    shinyjs.getHeight35 = function() {
    Shiny.onInputChange('height35', $(window).height() * .35);
    }
    shinyjs.getHeight10 = function() {
    Shiny.onInputChange('height10', $(window).height() * .1);
    }
    shinyjs.getHeight05 = function() {
    Shiny.onInputChange('height05', $(window).height() * .05);
    }
    "

    ui = fluidPage(

        #screen shouldn't go gray when plots are updating.
        tags$style(type="text/css", ".recalculating { opacity: 1.0; }" ),
        shinyjs::useShinyjs(),
        shinyjs::extendShinyjs(text=get_plotheight,
            functions=c('getHeight50', 'getHeight40', 'getHeight35',
                'getHeight10', 'getHeight10', 'getHeight05', 'init')),

        #hide error messages on screen
        tags$style(type="text/css",
            ".shiny-output-error { visibility: hidden; }",
            ".shiny-output-error:before { visibility: hidden; }"
        ),

        navbarPage(title=p(strong(a('StreamPULSE',
            href='https://data.streampulse.org/'))), inverse=TRUE,
            windowTitle='StreamPULSE Diagnostics',
            tabPanel('Model Performance',
                fluidRow(
                    column(10, offset=1, align='center',
                        plotOutput('KvQvER', height='auto')
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
                                # sliderInput("range", label=NULL,
                                #     min=DOYmin, max=DOYmax, value=c(DOYmin, DOYmax),
                                #     ticks=TRUE, step=6,
                                #     animate=animationOptions(interval=2000)
                                # )
                                htmlOutput('time_slider')
                            )
                        ),
                        hr()
                    )
                ),
                fluidRow(
                    column(9, align='center',
                        plotOutput('metab_legend', height='auto',
                            width='auto'),

                        plotOutput('metab_plot', height='auto', width='auto'),
                        plotOutput('O2_legend', height='auto',
                            width='auto'),
                        plotOutput('O2_plot', brush='O2_brush',
                            height='auto', width='auto')

                    ),
                    column(3, align='center',
                        plotOutput('cumul_legend', height='auto',
                            width='auto'),
                        plotOutput('cumul_plot', height='auto', width='auto'),
                        plotOutput('kernel_legend', height='auto',
                            width='auto'),
                        plotOutput('kernel_plot', height='auto', width='auto')
                    )
                ),
                br(),
                p(paste("*At this time, plotting code may fail if your data",
                    "do not all occur within the same calendar year.",
                    "If plots do not load, try toggling between tabs."),
                    style='color:gray; font-size:100%')
            )
        )
    )

    server = function(input, output, session){

        height50 = reactive({
            ifelse(is.null(input$height50), 0, input$height50)
        })
        height40 = reactive({
            ifelse(is.null(input$height40), 0, input$height40)
        })
        height35 = reactive({
            ifelse(is.null(input$height35), 0, input$height35)
        })
        height10 = reactive({
            ifelse(is.null(input$height10), 0, input$height10)
        })
        height05 = reactive({
            ifelse(is.null(input$height05), 0, input$height05)
        })

        js$getHeight50()
        js$getHeight40()
        js$getHeight35()
        js$getHeight10()
        js$getHeight05()

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
                }, height=height05)

                output$cumul_legend = renderPlot({
                    cumul_legend()
                }, height=height05)

                output$metab_plot = renderPlot({
                    ts_full = processing_func(fitpred$predictions, st=start,
                        en=end)
                    par(mar=c(1,4,0,1), oma=rep(0,4))
                    season_ts_func(ts_full, TRUE, st=start, en=end)
                }, height=height35)

                output$cumul_plot = renderPlot({
                    ts_full = processing_func(fitpred$predictions, st=start,
                        en=end)
                    par(mar=c(3,3.5,0.2,0.5), oma=rep(0,4))
                    cumulative_func(ts_full, st=start, en=end)
                }, height=height35)

                output$O2_legend = renderPlot({
                    O2_legend()
                }, height=height05)

                output$kernel_legend = renderPlot({
                    kernel_legend()
                }, height=height05)

                output$O2_plot = renderPlot({
                    par(mar=c(3,4,0,1), oma=rep(0,4))
                    O2_plot(mod_out=fitpred$mod_out, st=start, en=end,
                        input$O2_brush)
                }, height=height35)

                output$kernel_plot = renderPlot({
                    ts_full = processing_func(fitpred$predictions, st=start,
                        en=end)
                    par(mar=c(3,3.5,0,.5), oma=rep(0,4))
                    kernel_func(ts_full, 'Name and Year')
                }, height=height35)
            }
        })

        output$KvQvER = renderPlot({
            try(KvQvER_plot(mod_out=modOut), silent=TRUE)
        }, height=height50)

    }

    shinyApp(ui=ui, server=server)#, options=list(launch.browser=TRUE))
}
