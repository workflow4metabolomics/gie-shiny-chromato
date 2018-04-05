#-----------
# DEBUG MODE
file.create("/import/times.log", showWarnings = TRUE)
write("Starting debug mode", file="/import/times.log", append=TRUE)
write(as.character(Sys.time()), file="/import/times.log", append=TRUE)
write("Loading packages", file="/import/times.log", append=TRUE)
#-----------

# Load packages
library(shiny)
library(shinyWidgets)
library(xcms)
library(plotly)
library(RColorBrewer)
library(stringr)

#-----------
# DEBUG MODE
write(as.character(Sys.time()), file="/import/times.log", append=TRUE)
write("Get RSession and Import files (gx_get)", file="/import/times.log", append=TRUE)
#-----------

# Get RSession
load("/srv/shiny-server/samples/chromato_visu/inputdata.dat")
raw_files <- basename(rownames(xdata@phenoData@data))

groups <- unique(xdata@phenoData@data$sample_group)

# Import files by copying (identifier_type name because of filename and not hid)
# gx_get(raw_files, identifier_type='name')

# Symbolic link between /import and APP_PATH
#for (file in raw_files){
#  system(sprintf("ln -s /srv/shiny-server/data/datasets/%s %s", file, file))
#}

# In case of adjusted raws
#xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.6))
adjusted <- hasAdjustedRtime(xdata)

# Making a palette
default_palette <- rainbow(length(raw_files))

#----------
#DEBUG MODE
write(as.character(Sys.time()), file="/import/times.log", append=TRUE)
#----------


ui <- bootstrapPage(
	fluidRow(
		column(2,
			switchInput(
				inputId = "color_by_group",
				label = strong("Color_by_Group"),
				value = FALSE
			)
		),
		column(2,
			if (adjusted) {
				switchInput(
					inputId = "adjustedTime",
					label = strong("Adjusted_Time"),
					value = FALSE
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
				column(4,
					switchInput(
						inputId = "versus",
						label = strong("Versus_mode"),
						value = FALSE
					)
				),
				column(4,
					conditionalPanel(
						condition = "input.versus == true",
						selectInput(
							inputId = "group1", 
							label = "1st Group",
							choices = groups,
							width = "180px"
						)
					)
				),
				column(4,
                                        conditionalPanel(
                                                condition = "input.versus == true",
						uiOutput("versus_group")
                                        )
                                )        
			)
		),
	        column(2,
                        actionButton(
                                inputId = "draw",
                                label = "DRAW"
                        )
		)
	),
	fluidRow(
		plotlyOutput('TIC')
	),
	fluidRow(
		plotlyOutput('BPC')
	)
#	fluidRow(
#		div(style = "border-style: solid; border-width: thin;border-color: #000000",
#			plotlyOutput('TIC_chromato')
#		)
#	),
#	fluidRow(
#		div(style = "border-style: solid; border-width: thin; border-color: #000000",
#			plotlyOutput('BPC_chromato')
#		)
#	)
)

server <- function(input, output){

	# Calculate intensities
        tic_intensity <- split(tic(xdata, initial=FALSE), f = fromFile(xdata))
	bpc_intensity <- split(bpi(xdata, initial=FALSE), f = fromFile(xdata))

	# Display the select input for others groups
	output$versus_group <- renderUI({
		tagList(
			selectInput(
				inputId = "group2",
				label = "2nd Group",
				choices = subset(groups, !(groups %in% input$group1)),
				width = "180px"
			)
		)
	})

	# Get reactive values
	draw_chromato <- reactiveValues(value = 0)
	observeEvent(input$draw, {
		draw_chromato$value <- draw_chromato$value + 1
	}) 

	col_group <- eventReactive(input$draw, {
		col_group <- input$color_by_group
	})

	adjusted_time <- eventReactive(input$draw, {
		adjusted_time <- input$adjustedTime
	})

	versus_mode <- eventReactive(input$draw, {
		versus_mode <- input$versus
	})	

	pos_group <- eventReactive(input$draw, {
		pos_group <- input$group1
	})

	neg_group <- eventReactive(input$draw, {
		neg_group <- input$group2
	})

	color <- eventReactive(input$draw, {
	        if (col_group()){
			color=xdata@phenoData@data$sample_group
                } else {
                        color=xdata@phenoData@data$sample_name
		}
	})

	palette <- eventReactive(input$draw, {
                if (col_group()){
			palette <- brewer.pal(length(xdata@phenoData@data$sample_group), "Set1")
                } else {
			palette <- rainbow(length(raw_files))
                }
        })


	# Building chromatogram function
	build_chromato <- function(data, raws, title, intensity, draw_chromato, adjusted, adjusted_time, versus_mode, color, pos_group, neg_group) {

                if (draw_chromato != 0){

	                chromato <- plot_ly(source='alignmentChromato', type='scatter', mode='markers', colors=palette()) %>%
	                        layout(title=title, xaxis=list(title='Retention time (min)'), yaxis=list(title='Intensity'), showlegend=TRUE) %>%
	                        config(scrollZoom=TRUE, showLink=TRUE, displaylogo=FALSE,
	                                modeBarButtons=list(list('toImage', 'zoom2d', 'select2d', 'pan2d', 'autoScale2d', 'resetScale2d')))

	                if(is.null(raws)) return(chromato)
	                else if(!length(raws)) return(chromato)

                        # According to Adjusted Time
                        if (adjusted) {
                                if (adjusted_time) {
                                        rtime <- split(rtime(data, adjusted = TRUE)/60, f = fromFile(data))
                                } else {
                                        rtime <- split(rtime(data, adjusted = FALSE)/60, f = fromFile(data))
                                }
                        } else {
                                rtime <- split(rtime(data, adjusted = FALSE)/60, f = fromFile(data))
                        }

                        for(i in 1:length(raws)) {
                                index <- which(basename(rownames(phenoData(data))) == raws[i])
                                if (versus_mode) {
                                        if (data@phenoData@data$sample_group[index] %in% pos_group) {
                                                intens = intensity[[index]]
                                        } else if (data@phenoData@data$sample_group[index] %in% neg_group) {
                                                intens = -intensity[[index]]
                                        } else {
                                                intens = 0
                                        }
                                } else {
                                        intens = intensity[[index]]
                                }
                                chromato <- chromato %>% add_lines(
                                        x=rtime[[index]], y=intens, name=basename(raws[i]), hoverinfo='text', color=color()[i],
                                        text=paste('Intensity: ', round(intensity[[index]]), '<br />Retention Time: ', round(rtime[[index]], digits=2))
                                )
                        }
                } else {

                        chromato <- plot_ly(source='alignmentChromato', type='scatter', mode='markers', colors=default_palette) %>%
                                layout(title=title, xaxis=list(title='Retention time (min)'), yaxis=list(title='Intensity'), showlegend=TRUE) %>%
                                config(scrollZoom=TRUE, showLink=TRUE, displaylogo=FALSE,
                                        modeBarButtons=list(list('toImage', 'zoom2d', 'select2d', 'pan2d', 'autoScale2d', 'resetScale2d')))

                        if(is.null(raws)) return(chromato)
                        else if(!length(raws)) return(chromato)

                        rtime <- split(rtime(data, adjusted = FALSE)/60, f = fromFile(data))

                        for(i in 1:length(raws)) {
                                index <- which(basename(rownames(phenoData(data))) == raws[i])
                                intens = intensity[[index]]
                                chromato <- chromato %>% add_lines(
                                        x=rtime[[index]], y=intens, name=basename(raws[i]), hoverinfo='text', color=data@phenoData@data$sample_name[i],
                                        text=paste('Intensity: ', round(intensity[[index]]), '<br />Retention Time: ', round(rtime[[index]], digits=2))
                                )
                        }
               }

		return(chromato)
	}


	output$TIC <- renderPlotly({
		#-----------
		# DEBUG MODE
		write("Building TIC", file="/import/times.log", append=TRUE)
		#-----------

		chromato <- build_chromato(xdata, raw_files, "TIC", tic_intensity, draw_chromato$value, adjusted, adjusted_time(), versus_mode(), color(), pos_group(), neg_group())
		
		#-----------
		# DEBUG MODE  
		write(as.character(Sys.time()), file="/import/times.log", append=TRUE)
                write("Done", file="/import/times.log", append=TRUE)  
		#-----------

		return(chromato)
	})

	output$BPC <- renderPlotly({		
		#-----------
                # DEBUG MODE
                write("Building BPC", file="/import/times.log", append=TRUE)
		#-----------

                chromato <- build_chromato(xdata, raw_files, "BPC", bpc_intensity, draw_chromato$value, adjusted, adjusted_time(), versus_mode(), color(), pos_group(), neg_group())

                #-----------
		# DEBUG MODE
                write(as.character(Sys.time()), file="/import/times.log", append=TRUE)
                write("Done", file="/import/times.log", append=TRUE)
		#-----------

		return(chromato)
	})

	# Not used now
	if(FALSE){
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
}

shinyApp(ui, server)
