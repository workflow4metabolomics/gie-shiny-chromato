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
library(webshot)

#-----------
# DEBUG MODE
write(as.character(Sys.time()), file="/import/times.log", append=TRUE)
write("Get RSession and Import files", file="/import/times.log", append=TRUE)
#-----------

# Get RSession
load("/srv/shiny-server/samples/chromato_visu/inputdata.dat")
raw_names <- xdata@phenoData@data$sample_name
raw_group <- xdata@phenoData@data$sample_group

## Settings
# Graph settings
height <- "600"
# Samples number to display
samples_to_display <- 50

# Get group names sorted
groups <- sort(names(table(raw_group)))

## Import files by copying them, not used because to slow (identifier_type='name' because of filename and not hid)
#gx_get(raw_names, identifier_type='name')

# In case of adjusted raws (retcor)
adjusted <- hasAdjustedRtime(xdata)

# Making a color palette
default_palette <- rainbow(length(raw_names))

#----------
#DEBUG MODE
write(as.character(Sys.time()), file="/import/times.log", append=TRUE)
#----------


ui <- bootstrapPage(

	fluidRow(
		includeCSS("styles.css"),
		column(3,
			h5(strong("Chromatogram displayed :")),
			switchInput(
				inputId = "input_chromatogram",
				label = strong("Chromatogram"),
				value = FALSE,
				offLabel = "BPC",
				onLabel = "TIC"
			)
		),
		column(3,
			h5(strong("Intensity :")),
			switchInput(
				inputId = "intensity",
				label = strong("Intensity"),
				value = FALSE,
				offLabel = "Absolute",
				onLabel = "Relative"
			)
		),
		column(3,
			if (adjusted) {
				uiOutput("display_adjusted")
			}
		),
		column(2),
	        column(1,
			br(),
            actionButton(
                inputId = "draw",
                label = "DRAW",
                class = "btn-primary"
            )
		)
	),
	fluidRow(
		column(3,
			h5(strong("Group displayed :")),
			selectInput(
				inputId = "select_group",
				label = NULL,
				choices = groups,
				selected = groups,
				multiple = TRUE,
				selectize = TRUE,
				width = '200px'
			)
		),
		column(3,
			h5(strong("Color by Group :")),
			switchInput(
				inputId = "color_by_group",
				label = strong("Color_by_Group"),
				value = FALSE
			)
		),
		column(6)
	),
	fluidRow(
		style = "height : 110px",
		column(3,
			h5(strong("Versus mode :")),
			switchInput(
				inputId = "versus",
				label = strong("Versus"),
				value = FALSE
			)
		),
		column(3,
			conditionalPanel(
				condition = "input.versus == true",
				h5(strong("Group or Sample :")),
				switchInput(
					inputId = "versus_by",
					label = strong("Versus"),
					value = FALSE,
					offLabel = "Group",
					onLabel = "Sample"
				)
			)
		),
		column(6,
			fluidRow(
				conditionalPanel(
					condition = "input.versus == true && input.versus_by == false",
					column(5,
						h5(strong("Upper Group :")),
						uiOutput("upper_group")
					)
				),
                conditionalPanel(
                    condition = "input.versus == true && input.versus_by == false",
                    column(5,
						h5(strong("Under Group :")),
						uiOutput("under_group")
                    )
                ),
                conditionalPanel(
					condition = "input.versus == true && input.versus_by == true",					
					column(5,
						h5(strong("Upper Sample(s) :")),
						uiOutput("upper_sample")
					)
				),
				conditionalPanel(
                    condition = "input.versus == true && input.versus_by == true",
					column(5,				
						h5(strong("Under Sample(s) :")),
						uiOutput("under_sample")
                    )
				),
				column(2)
			)
		)
	),
	fluidRow(
		style = paste("height:", height, "px", sep=""), 
		column(2,
			style = "height : 100%",
			h5(strong("Sample List :")),
			uiOutput("sample_list")
		),
		column(10,
			plotlyOutput('CHROM')
		)
	),
	fluidRow(
		column(10),
		column(2,
            actionButton(
            	icon = icon("export", lib = "glyphicon"),
                inputId = "export",
                label = "EXPORT"
            )
		)
	)
)

server <- function(input, output){

	# Get pre-calculate chromatogram from chromTIC and chromBPI objects
	tic_chrom <- chromTIC
	bpc_chrom <- chromBPI
	if (adjusted) {
		tic_chrom_adjusted <- chromTIC_adjusted
		bpc_chrom_adjusted <- chromBPI_adjusted
	}


	# On-click Draw button action
	draw_chromato <- reactiveValues(value = 0)
	observeEvent(input$draw, {
		draw_chromato$value <- draw_chromato$value + 1
	}) 

	# Chromatogram displayed
    which_chromatogram <- eventReactive(input$draw, {
        which_chromatogram <- input$input_chromatogram
    })

	# Intensity
	relative_intensity <- eventReactive(input$draw, {
		relative_intensity <- input$intensity
	})

	# Display Adjusted title
	output$display_adjusted <- renderUI({
		tagList(
			h5(strong("Adjusted RTime :")),
			switchInput(
				inputId = "adjustedTime",
				label = strong("Adjusted"),
				value = FALSE
			)	
		)
	})

	# Adjusted RTime
	adjusted_time <- eventReactive(input$draw, {
		adjusted_time <- input$adjustedTime
	})

	# Group displayed
	selected_groups <- eventReactive(input$draw, {
		selected_groups <- input$select_group
	})

	# Coloration
	col_group <- eventReactive(input$draw, {
		col_group <- input$color_by_group
	})

	palette <- eventReactive(input$draw, {
        if (col_group()){
			palette <- brewer.pal(length(table(raw_group)), "Set1")
        } else {
			palette <- rainbow(length(raw_names))
        }
    })

	#color <- eventReactive(input$draw, {
	    #if (col_group()){
		#	color <- raw_group
        #} else {
        #	color <- raw_names
		#}
	#})

	output$upper_group <- renderUI({
		tagList(
			selectInput(
				inputId = "group1",
				label = NULL,
				choices = subset(groups, (groups %in% input$select_group)),
				width = "180px"
			)	
		)
	})

	output$under_group <- renderUI({
		tagList(
			selectInput(
				inputId = "group2",
				label = NULL,
				choices = subset(input$select_group, !(input$select_group %in% input$group1)),
				width = "180px"
			)
		)
	})

	# TODO : Selectize TRUE with limited nb of sample or FALSE (ctrl + choice) ...
	output$upper_sample <- renderUI({
		tagList(
			selectInput(
				inputId = "sample1",
				label = NULL,
				choices = input$select_sample,
				multiple = TRUE,
				selectize = FALSE,
				width = "180px"
			)	
		)
	})

	output$under_sample <- renderUI({
		tagList(
			selectInput(
				inputId = "sample2",
				label = NULL,
				choices = subset(input$select_sample, !(input$select_sample %in% input$sample1)),
				multiple = TRUE,
				selectize = FALSE,
				width = "180px"
			)
		)
	})

	versus_mode <- eventReactive(input$draw, {
		versus_mode <- input$versus
	})

	versus_by <- eventReactive(input$draw, {
		versus_by <- input$versus_by
	})

	pos_group <- eventReactive(input$draw, {
		pos_group <- input$group1
	})

	neg_group <- eventReactive(input$draw, {
		neg_group <- input$group2
	})

	pos_sample <- eventReactive(input$draw, {
		pos_sample <- input$sample1
	})

	neg_sample <- eventReactive(input$draw, {
		neg_sample <- input$sample2
	})

	output$sample_list <- renderUI({
		tagList(
			checkboxGroupInput(
				inputId = "select_sample",
				label = NULL,
				choices = sort(subset(raw_names, (raw_group %in% input$select_group))),
				selected = sort(subset(raw_names, (raw_group %in% input$select_group))),
				width = '200px'
			)
		)
	})

	# Building chromatogram function
	build_chromato <- function(raw_names, raw_group, which_chromato, draw_chromato, groups_selected, adjusted, adjusted_time, relative_intensity, col_group, versus_mode, versus_by, pos_group, neg_group, pos_sample, neg_sample) {

        if (draw_chromato != 0){

			if (which_chromato) {
				title <- "TIC"
			} else {
				title <- "BPC"
			}

			# Initialize a blank chromatogram
	        displayed_chromatogram <- plot_ly(source='alignmentChromato', type='scatter', mode='markers', colors=palette()) %>%
	            layout(title=title, height=height, xaxis=list(title='Retention time (min)'), yaxis=list(title='Intensity'), showlegend=TRUE) %>%
	            config(scrollZoom=TRUE, showLink=TRUE, displaylogo=FALSE,
	            modeBarButtons=list(list('toImage', 'zoom2d', 'select2d', 'pan2d', 'autoScale2d', 'resetScale2d')))

	        if(is.null(raw_names)) return(displayed_chromatogram)
	        else if(!length(raw_names)) return(displayed_chromatogram)

	        # According to Adjusted Time
			# If retcor has been done
	        if (adjusted) {
				# If Adjusted Option is TRUE
	            if (adjusted_time) {
					if (title=="TIC") {
		                chrom <- tic_chrom_adjusted
					} else if (title=="BPC") {
		                chrom <- bpc_chrom_adjusted
					}
	            } else {
					if (title=="TIC") {
		                chrom <- tic_chrom
					} else if (title=="BPC") {
		                chrom <- bpc_chrom
					}
				}
	        } else {
				if (title=="TIC") {
	                chrom <- tic_chrom
				} else if (title=="BPC") {
	                chrom <- bpc_chrom
				}
	        }

			# Initialize a variable in case of little group of samples
			files_to_add <- 0

            for ( j in 1:length(groups) ) {
				if (groups[j] %in% groups_selected){

					# Get the samples to display
	                files_in_group <- length(raw_names[raw_group == groups[j]])
	                files_to_get <- samples_to_display%/%(length(groups_selected)*2)
	                if ( files_in_group < (files_to_get*2) ) {
						files_to_add <- (files_to_get*2)-files_in_group
	                    files_to_get <- files_in_group
	                } else { 
	                    files_to_get <- files_to_get + files_to_add
						files_to_add <- 0 
					}

	                # Counter initialization
	                group_file_nb <- 1
	
	                for ( i in 1:length(raw_names) ) {
	                    #Check if file is in group[j]
	                    if ( raw_group[i] == groups[j] ) {
			                index <- which(raw_names == raw_names[i])
			                intens <- chrom[[index]]@intensity
							intens_max <- max(chrom[[index]]@intensity)

							# In case of versus mode
			                if (versus_mode) {
			                	if (versus_by) { #sample
			    	            	if (raw_names[index] %in% neg_sample) {
						                intens <- -intens
				                    } else if ( !(raw_names[index] %in% neg_sample) && !(raw_names[index] %in% pos_sample) ) {
				                        intens <- 0
				                	}
			                	} else {
			    	            	if (raw_group[index] %in% neg_group) {
						                intens <- -intens
				                    } else if ( !(raw_group[index] %in% neg_group) && !(raw_group[index] %in% pos_group) ) {
				                        intens <- 0
				                	}	
			                	}          
			                }

							# In case of relative intensity
							if (relative_intensity) {
								intens <- intens*100/intens_max
							}

							# In case of color by group
					        if(col_group){
								# Building the chromatogram
	                	        if ( group_file_nb <= files_to_get ) {
		        	                displayed_chromatogram <- displayed_chromatogram %>% add_lines(
		    	                        x=chrom[[index]]@rtime/60, y=intens, name=raw_names[i], hoverinfo='text', color=raw_group[i],
			                            text=paste('Name: ', raw_names[i], '<br />Intensity: ', round(chrom[[index]]@intensity), '<br />Retention Time: ', round(chrom[[index]]@rtime/60, digits=2))
			                        )
		                        } else if ( group_file_nb > (files_in_group - files_to_get) ) {
		                        	displayed_chromatogram <- displayed_chromatogram %>% add_lines(
		                    	        x=chrom[[index]]@rtime/60, y=intens, name=raw_names[i], hoverinfo='text', color=raw_group[i],
		                	            text=paste('Name: ', raw_names[i], '<br />Intensity: ', round(chrom[[index]]@intensity), '<br />Retention Time: ', round(chrom[[index]]@rtime/60, digits=2))
		            	            )
	            	            } else {
		    	                    displayed_chromatogram <- displayed_chromatogram %>% add_lines(
			                            x=chrom[[index]]@rtime/60, y=intens, name=raw_names[i], hoverinfo='text', visible="legendonly", color=raw_group[i],
			                            text=paste('Name: ', raw_names[i], '<br />Intensity: ', round(chrom[[index]]@intensity), '<br />Retention Time: ', round(chrom[[index]]@rtime/60, digits=2))
			                        )
	                        	}
					        } else {
					        	# Building the chromatogram
		                        if ( group_file_nb <= files_to_get ) {
			                        displayed_chromatogram <- displayed_chromatogram %>% add_lines(
			                            x=chrom[[index]]@rtime/60, y=intens, hoverinfo='text', color=raw_names[i],
			                            text=paste('Name: ', raw_names[i], '<br />Intensity: ', round(chrom[[index]]@intensity), '<br />Retention Time: ', round(chrom[[index]]@rtime/60, digits=2))
			                        )
		                        } else if ( group_file_nb > (files_in_group - files_to_get) ) {
			                        displayed_chromatogram <- displayed_chromatogram %>% add_lines(
			                            x=chrom[[index]]@rtime/60, y=intens, hoverinfo='text', color=raw_names[i],
			                            text=paste('Name: ', raw_names[i], '<br />Intensity: ', round(chrom[[index]]@intensity), '<br />Retention Time: ', round(chrom[[index]]@rtime/60, digits=2))
			                        )
		                        } else {
			                        displayed_chromatogram <- displayed_chromatogram %>% add_lines(
			                            x=chrom[[index]]@rtime/60, y=intens, hoverinfo='text', visible="legendonly", color=raw_names[i],
			                            text=paste('Name: ', raw_names[i], '<br />Intensity: ', round(chrom[[index]]@intensity), '<br />Retention Time: ', round(chrom[[index]]@rtime/60, digits=2))
			                        )
		                        }
					        }

	                        group_file_nb <- group_file_nb + 1
	                    }
	                }
				}
            }

        } else {

			write(as.character(Sys.time()), file="/import/times.log", append=TRUE)
			write("Initialisation", file="/import/times.log", append=TRUE)

			# Initialize a blank chromatogram
            displayed_chromatogram <- plot_ly(source='alignmentChromato', type='scatter', mode='markers', colors=default_palette) %>%
                layout(title="BPC", height=height, xaxis=list(title='Retention time (min)'), yaxis=list(title='Intensity'), showlegend=TRUE) %>%
                config(scrollZoom=TRUE, showLink=TRUE, displaylogo=FALSE,
                modeBarButtons=list(list('toImage', 'zoom2d', 'select2d', 'pan2d', 'autoScale2d', 'resetScale2d')))

            if(is.null(raw_names)) return(displayed_chromatogram)
            else if(!length(raw_names)) return(displayed_chromatogram)

			chrom <- bpc_chrom

			# Initialize a variable in case of little group of samples
			files_to_add <- 0

			write(as.character(Sys.time()), file="/import/times.log", append=TRUE)
			write("Parcours de chaque groupe", file="/import/times.log", append=TRUE)


			for ( j in 1:length(groups) ) {

				# Get the samples to display
				files_in_group <- length(raw_names[raw_group == groups[j]])
				files_to_get <- samples_to_display%/%(length(groups)*2)
				if ( files_in_group < (files_to_get*2) ) {
					files_to_add <- (files_to_get*2)-files_in_group
					files_to_get <- files_in_group
				} else {
					files_to_get <- files_to_get + files_to_add
					files_to_add <- 0
				}

                # Counter initialization
				group_file_nb <- 1

				write("Parcours de chaque echantillon", file="/import/times.log", append=TRUE)

				for ( i in 1:length(raw_names) ) {
					#Check if file is in group[j]
					if ( raw_group[i] == groups[j] ) {
		                index <- which(raw_names == raw_names[i])
						intens <- chrom[[index]]@intensity

						# Building the chromatogram
						if ( group_file_nb <= files_to_get ) {
                            displayed_chromatogram <- displayed_chromatogram %>% add_lines(
                                x=chrom[[index]]@rtime/60, y=intens, hoverinfo='text', color=raw_names[i],
                                text=paste('Name: ', raw_names[i], '<br />Intensity: ', round(chrom[[index]]@intensity), '<br />Retention Time: ', round(chrom[[index]]@rtime/60, digits=2))
                            )
						} else if ( group_file_nb > (files_in_group - files_to_get) ) {
                            displayed_chromatogram <- displayed_chromatogram %>% add_lines(
                                x=chrom[[index]]@rtime/60, y=intens, hoverinfo='text', color=raw_names[i],
                                text=paste('Name: ', raw_names[i], '<br />Intensity: ', round(chrom[[index]]@intensity), '<br />Retention Time: ', round(chrom[[index]]@rtime/60, digits=2))
                            )
						} else {
				            displayed_chromatogram <- displayed_chromatogram %>% add_lines(
				                x=chrom[[index]]@rtime/60, y=intens, hoverinfo='text', visible="legendonly", color=raw_names[i],
				                text=paste('Name: ', raw_names[i], '<br />Intensity: ', round(chrom[[index]]@intensity), '<br />Retention Time: ', round(chrom[[index]]@rtime/60, digits=2))
				            )
						}

						group_file_nb <- group_file_nb + 1
					}
				}

				write(as.character(Sys.time()), file="/import/times.log", append=TRUE)

			}

			write("Fin de parcours des groupes", file="/import/times.log", append=TRUE)
			write(as.character(Sys.time()), file="/import/times.log", append=TRUE)

        }

		return(displayed_chromatogram)
	}

	displayed_chromatogram <- reactive({
		reactive_chromatogram <- build_chromato(raw_names, raw_group, which_chromatogram(), draw_chromato$value, selected_groups(), adjusted, adjusted_time(), relative_intensity(), col_group(), versus_mode(), versus_by(), pos_group(), neg_group(), pos_sample(), neg_sample())
	})		

	output$CHROM <- renderPlotly({
		#-----------
		# DEBUG MODE
		write("------------", file="/import/times.log", append=TRUE)
		write("Building Chromatogram", file="/import/times.log", append=TRUE)
		write("------------", file="/import/times.log", append=TRUE)
		#-----------

		graph <- displayed_chromatogram()

		#-----------
		# DEBUG MODE  
		write("------------", file="/import/times.log", append=TRUE)
		write(as.character(Sys.time()), file="/import/times.log", append=TRUE)
		write("------------", file="/import/times.log", append=TRUE)
		#-----------

		return(graph)
	})

	observeEvent(input$export, {
		file_name <- "test.tsv"
		file_content<-data.frame(a=1:3,b=4:6)
		write.table(file_content, file=file_name, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
		gx_put(file_name)
	})

	## TO EXPORT PNG IMAGE OF THE GRAPH IN THE HISTORY ##
	#observeEvent(input$export, {
	#	tmpFile <- tempfile(pattern = "chrom_", tmpdir = tempdir(), fileext = ".png")
	#	export(displayed_chromatogram(), file = tmpFile)
	#	browseURL(tmpFile)
	#	png(filename="plot.png")
	#	displayed_chromatogram()
	#	dev.off()
	#	path <- getwd()
	#	file <- paste0(path, "/", "plot.png")
	#	gx_put(tmpFile)
	#})

}

shinyApp(ui, server)
