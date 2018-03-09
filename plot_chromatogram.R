# Load packages
library(shiny)
library(shinyWidgets)
library(xcms)
library(plotly)
library(RColorBrewer)
library(stringr)

# Get RSession
load("/srv/shiny-server/data/inputdata.dat")
raw_files <- basename(rownames(xdata@phenoData@data))

# identifier_type name because of filename and not hid
gx_get(raw_files, identifier_type='name')

for (file in raw_files){
  system(sprintf("ln -s /import/%s %s", file, file))
}

# In case of adjusted raws
xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.6))
adjusted <- hasAdjustedRtime(xdata)

# Making a palette
palette <- brewer.pal(length(raw_files), "Dark2")

# Versus
versus_group <- c(which(xdata@phenoData@data$sample_group == "KO"))

ui <- bootstrapPage(
	fluidRow(
		column(2,
			if (adjusted) {
				switchInput(
					inputId = "adjustedTime",
					label = strong("Adjusted Time"),
					value = TRUE
				)
			}
		),
		column(2,
			switchInput(
				inputId = "color_by_group",
				label = strong("Color by Group"),
				value = FALSE,
			)
		),
		column(2,
			switchInput(
				inputId = "versus",
				label = strong("Versus mode"),
				value = FALSE
			)
		)
	),
	fluidRow(
		column(
			6,
			div(style = "border-style: solid; border-width: thin; border-color: #000000",
				plotlyOutput('TIC')
			)
		),
		column(
			6,
			div(style = "border-style: solid; border-width: thin;border-color: #000000",
				plotlyOutput('TIC_chromato')
			)
		)
	),
	fluidRow(
		column(
			6,
			div(style = "border-style: solid; border-width: thin; border-color: #000000",
				plotlyOutput('BPC')
			)
		),
		column(
			6,
			div(style = "border-style: solid; border-width: thin; border-color: #000000",
				plotlyOutput('BPC_chromato')
			)
		)
	)
)

server <- function(input, output){

        # Color by Group
	color <- reactive({
	        if (input$color_by_group){
			color=xdata@phenoData@data$sample_group
                } else {
                        color=xdata@phenoData@data$sample_name
		}
	})

	output$TIC <- renderPlotly({

                # According to Adjusted Time
		if (adjusted) {
			if (input$adjustedTime) {
				title <- "TIC adjusted with tic()"
				rtime <- split(rtime(xdata, adjusted = TRUE)/60, f = fromFile(xdata))
			} else {
				title <- "TIC with tic()"
				rtime <- split(rtime(xdata, adjusted = FALSE)/60, f = fromFile(xdata))
			}
		} else {
			title <- "TIC with tic()"
			rtime <- split(rtime(xdata, adjusted = FALSE)/60, f = fromFile(xdata))
		}
		
                #function tic of MsnBase faster than chromatogram of XCMS (for me, to check)
		intensity <- split(tic(xdata, initial=FALSE), f = fromFile(xdata))

		chromato <- plot_ly(source='alignmentChromato', type='scatter', mode='markers', colors=palette) %>% 
			layout(title=title, xaxis=list(title='Retention time (min)'), yaxis=list(title='Intensity'), showlegend=TRUE) %>% 
			config(scrollZoom=TRUE, showLink=TRUE, displaylogo=FALSE, 
				modeBarButtons=list(list('toImage', 'zoom2d', 'select2d', 'pan2d', 'autoScale2d', 'resetScale2d')))
		
		if(is.null(raw_files)) return(chromato)
		else if(!length(raw_files)) return(chromato)

		for(i in 1:length(raw_files)) {
			index <- which(basename(rownames(phenoData(xdata))) == raw_files[i])
			intens = intensity[[index]]
			if ((index %in% versus_group) & (input$versus)) {
				intens = -intens
			}
			chromato <- chromato %>% add_lines(
				x=rtime[[index]], y=intens, name=basename(raw_files[i]), hoverinfo='text', color=color()[i],
				text=paste('Intensity: ', round(intensity[[index]]), '<br />Retention Time: ', round(rtime[[index]], digits=2))
			)
		}

		return(chromato)
	})

        output$TIC_chromato <- renderPlotly({

                # According to Adjusted Time
		if (adjusted) {
	                if (input$adjustedTime) {
				title <- "TIC adjusted with chromatogram()"
                                rtime <- chromatogram(xdata, aggregationFun = 'sum', adjustedRtime = TRUE)
                        } else {
                                title <- "TIC with chromatogram()"
                                rtime <- chromatogram(xdata, aggregationFun = 'sum', adjustedRtime = FALSE)
                        }
		} else {
                        title <- "TIC with chromatogram()"
			rtime <- chromatogram(xdata, aggregationFun = 'sum', adjustedRtime = FALSE)
                }

                chromato <- plot_ly(source='alignmentChromato', type='scatter', mode='markers', colors=palette) %>%
                        layout(title=title, xaxis=list(title='Retention time (min)'), yaxis=list(title='Intensity'), showlegend=TRUE) %>%
                        config(scrollZoom=TRUE, showLink=TRUE, displaylogo=FALSE,
                                modeBarButtons=list(list('toImage', 'zoom2d', 'select2d', 'pan2d', 'autoScale2d', 'resetScale2d')))

                if(is.null(raw_files)) return(chromato)
                else if(!length(raw_files)) return(chromato)

                for(i in 1:length(raw_files)) {
                        index <- which(basename(rownames(phenoData(xdata))) == raw_files[i])
                        intens = rtime[[index]]@intensity
                        if ((index %in% versus_group) & (input$versus)) {
                                intens = -intens
                        }
                        chromato <- chromato %>% add_lines(
                                x=rtime[[index]]@rtime/60, y=intens, name=basename(raw_files[i]), hoverinfo='text', color=color()[i],
                                text=paste('Intensity: ', round(rtime[[index]]@intensity), '<br />Retention Time: ', round(rtime[[index]]@rtime/60, digits=2))
                        )
                }

                return(chromato)
        })


	
	output$BPC <- renderPlotly({

                # According to Adjusted Time
		if (adjusted) {
                        if (input$adjustedTime) {
                                title <- "BPC adjusted with bpi()"
                                rtime <- split(rtime(xdata, adjusted = TRUE)/60, f = fromFile(xdata))
                        } else {
                                title <- "BPC with bpi()"
                                rtime <- split(rtime(xdata, adjusted = FALSE)/60, f = fromFile(xdata))
                        }
                } else {
                        title <- "BPC with bpi()"
                        rtime <- split(rtime(xdata, adjusted = FALSE)/60, f = fromFile(xdata))
                }

                #function tic of MsnBase faster than chromatogram of XCMS (for me, to check)
                intensity <- split(bpi(xdata, initial = FALSE), f = fromFile(xdata))

		chromato <- plot_ly(source='alignmentChromato', type='scatter', mode='markers', colors=palette) %>% 
			layout(title=title, xaxis=list(title='Retention time (min)'), yaxis=list(title='Intensity'), showlegend=TRUE) %>% 
			config(scrollZoom=TRUE, showLink=TRUE, displaylogo=FALSE, 
				modeBarButtons=list(list('toImage', 'zoom2d', 'select2d', 'pan2d', 'autoScale2d', 'resetScale2d')))

		if(is.null(raw_files)) return(chromato)
		else if(!length(raw_files)) return(chromato)

                for(i in 1:length(raw_files)) {
                        index <- which(basename(rownames(phenoData(xdata))) == raw_files[i])
                        intens = intensity[[index]]
                        if ((index %in% versus_group) & (input$versus)) {
                                intens = -intens
                        }
                        chromato <- chromato %>% add_lines(
                                x=rtime[[index]], y=intens, name=basename(raw_files[i]), hoverinfo='text', color=color()[i],
                                text=paste('Intensity: ', round(intensity[[index]]), '<br />Retention Time: ', round(rtime[[index]], digits=2))
                        )
                }

		return(chromato)
	})

        output$BPC_chromato <- renderPlotly({

                # According to Adjusted Time
		if (adjusted) {
                        if (input$adjustedTime) {
                                title <- "BPC adjusted with chromatogram()"
	                        rtime <- chromatogram(xdata, aggregationFun = 'max', adjustedRtime = TRUE)
                        } else {
                                title <- "BPC with chromatogram()"
	                        rtime <- chromatogram(xdata, aggregationFun = 'max', adjustedRtime = FALSE)
                        }
                } else {
                        title <- "BPC with chromatogram()"
			rtime <- chromatogram(xdata, aggregationFun = 'max', adjustedRtime = FALSE)
                }

                chromato <- plot_ly(source='alignmentChromato', type='scatter', mode='markers', colors=palette) %>%
                        layout(title=title, xaxis=list(title='Retention time (min)'), yaxis=list(title='Intensity'), showlegend=TRUE) %>%
                        config(scrollZoom=TRUE, showLink=TRUE, displaylogo=FALSE,
                                modeBarButtons=list(list('toImage', 'zoom2d', 'select2d', 'pan2d', 'autoScale2d', 'resetScale2d')))

                if(is.null(raw_files)) return(chromato)
                else if(!length(raw_files)) return(chromato)

                for(i in 1:length(raw_files)) {
                        index <- which(basename(rownames(phenoData(xdata))) == raw_files[i])
                        intens = rtime[[index]]@intensity
                        if ((index %in% versus_group) & (input$versus)) {
                                intens = -intens
                        }
                        chromato <- chromato %>% add_lines(
                                x=rtime[[index]]@rtime/60, y=intens, name=basename(raw_files[i]), hoverinfo='text', color=color()[i],
                                text=paste('Intensity: ', round(rtime[[index]]@intensity), '<br />Retention Time: ', round(rtime[[index]]@rtime/60, digits=2))
                        )
                }

                return(chromato)
        })
}

shinyApp(ui, server)
