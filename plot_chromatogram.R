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

#group_colors <- brewer.pal(length(raw_files), "Set1")[1:length(unique(xdata$sample_group))]
#names(group_colors) <- unique(xdata$sample_group)

ui <- bootstrapPage(
	fluidRow(
		column(2,
			switchInput(
				inputId = "color_by_group",
				label = strong("Color_by_Group"),
				value = FALSE
			)
		),
		column(4,
			if (adjusted) {
				switchInput(
					inputId = "adjustedTime",
					label = strong("Adjusted_Time"),
					value = TRUE
				)
			}
		),
		column(6,
			fluidRow(
				tags$head(
					tags$style(
						type="text/css", 
						"label.control-label, 
						.selectize-control.single {display: table-cell; vertical-align: middle; width: 80px;} 
						.form-group {display: table-row;} 
						.selectize-control {margin-bottom: 10px;}"
					)
				),
				column(3,
					switchInput(
						inputId = "versus",
						label = strong("Versus_mode"),
						value = FALSE
					)
				),
				column(3,
					conditionalPanel(
						condition = "input.versus == true",
						selectInput(
							inputId = "group1", 
							label = "1st Group",
							choices = unique(xdata@phenoData@data$sample_group),
							width = "180px"
						)
					)
				),
				column(3,
                                        conditionalPanel(
                                                condition = "input.versus == true",
						uiOutput("versus_group")
                                        )
                                ),
                                column(3,
                                        conditionalPanel(
                                                condition = "input.versus == true",
                                                actionButton(
                                                        inputId = "apply_versus",
                                                        label = "Run Vs"
                                                )
                                        )
                                )
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
			#color=rainbow(length(xdata@phenoData@data$sample_name))
		}
	})

	output$versus_group <- renderUI({
		tagList(
			selectInput(
				inputId = "group2",
				label = "2nd Group",
				choices = subset(unique(xdata@phenoData@data$sample_group), !(unique(xdata@phenoData@data$sample_group) %in% input$group1)),
				width = "180px"
			)
		)
	})

	pos_group <- eventReactive(input$apply_versus, {
		pos_group <- c(input$group1)
	})
	neg_group <- eventReactive(input$apply_versus, {
		neg_group <- c(input$group2)
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
			if (input$versus) {
				if (xdata@phenoData@data$sample_group[index] %in% pos_group()) {
					intens = intensity[[index]]
				} else if (xdata@phenoData@data$sample_group[index] %in% neg_group()) {
					intens = -intensity[[index]]
				} else {
					intens = intensity[[index]]
				}
			} else {
				intens = intensity[[index]]
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
                        if (input$versus) {
                                if (xdata@phenoData@data$sample_group[index] %in% pos_group()) {
		                        intens = rtime[[index]]@intensity
                                } else if (xdata@phenoData@data$sample_group[index] %in% neg_group()) {
                                        intens = -rtime[[index]]@intensity
                                } else {
		                        intens = rtime[[index]]@intensity
                                }
                        } else {
				intens = rtime[[index]]@intensity
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
                        if (input$versus) {
                                if (xdata@phenoData@data$sample_group[index] %in% pos_group()) {
                                        intens = intensity[[index]]
                                } else if (xdata@phenoData@data$sample_group[index] %in% neg_group()) {
                                        intens = -intensity[[index]]
                                } else {
                                        intens = intensity[[index]]
                                }
                        } else {
                                intens = intensity[[index]]
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
                        if (input$versus) {
                                if (xdata@phenoData@data$sample_group[index] %in% pos_group()) {
                                        intens = rtime[[index]]@intensity
                                } else if (xdata@phenoData@data$sample_group[index] %in% neg_group()) {
                                        intens = -rtime[[index]]@intensity
                                } else {
                                        intens = rtime[[index]]@intensity
                                }
                        } else {
                                intens = rtime[[index]]@intensity
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
