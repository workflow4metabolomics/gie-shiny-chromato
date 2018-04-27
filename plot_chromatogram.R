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
write("Get RSession and Import files", file="/import/times.log", append=TRUE)
#-----------

# Get RSession
load("/srv/shiny-server/samples/chromato_visu/inputdata.dat")
raw_files <- basename(rownames(xdata@phenoData@data))

## Settings
# Graph settings
height <- "600"
# Samples number to display
samples_to_display <- 50


# Get group names sorted by croissant order
groups <- names(sort(table(xdata@phenoData@data$sample_group)))
# AJOUTER GROUP AVEC 1ST LETTRE MAJ POUR LES OPTIONS

## Import files by copying them, not used because to slow (identifier_type='name' because of filename and not hid)
#gx_get(raw_files, identifier_type='name')

# In case of adjusted raws (retcor)
adjusted <- hasAdjustedRtime(xdata)

# Making a color palette
default_palette <- rainbow(length(raw_files))

#----------
#DEBUG MODE
write(as.character(Sys.time()), file="/import/times.log", append=TRUE)
#----------


ui <- bootstrapPage(

# Used to modify the CSS
#	tags$head(
#		tags$style(
#			type="text/css", 
#			"label.control-label, 
#			.selectize-control.single {display: table-cell; vertical-align: middle; width: 80px;} 
#			.form-group {display: table-row;} 
#			.selectize-control {margin-bottom: 10px;}"
#		)
#	),

	fluidRow(
		style = 'margin:0px',
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
				#h5(strong("Adjusted RTime :")), => NOT WORK
				switchInput(
					inputId = "adjustedTime",
					label = strong("Adjusted_Time"),
					value = FALSE
				)
			}
		),
		column(2),
	        column(1,
			br(),
                        actionButton(
                                inputId = "draw",
                                label = "DRAW"
                        )
		)
	),
	br(),
	fluidRow(
		style = 'margin:0px',
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
		column(6,
			fluidRow(
				column(4,
					h5(strong("Switch to Versus mode :")),
					switchInput(
						inputId = "versus",
						label = strong("Versus"),
						value = FALSE
					)
				),
				column(3,
					conditionalPanel(
						condition = "input.versus == true",
						h5(strong("First Group :")),
						selectInput(
							inputId = "group1", 
							label = NULL,
							choices = groups,
							width = "180px"
						)
					)
				),
				column(3,
                                        conditionalPanel(
                                                condition = "input.versus == true",
						h5(strong("Second Group :")),
						uiOutput("versus_group")
                                        )
                                ),
				column(2)
			)
		)
	),
	fluidRow(
		style = paste('margin:0px; height:', height, 'px', sep=''),
		plotlyOutput('CHROM')
	),
	fluidRow(
		column(10),
		column(2,
			br(),
                        actionButton(
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

	color <- eventReactive(input$draw, {
	        if (col_group()){
			color <- xdata@phenoData@data$sample_group
                } else {
                        color <- xdata@phenoData@data$sample_name
		}
	})

	palette <- eventReactive(input$draw, {
                if (col_group()){
			palette <- brewer.pal(length(xdata@phenoData@data$sample_group), "Set1")
                } else {
			palette <- rainbow(length(raw_files))
                }
        })

	# Display the select input for others groups
	output$versus_group <- renderUI({
		tagList(
			selectInput(
				inputId = "group2",
				label = NULL,
				choices = subset(groups, !(groups %in% input$group1)),
				width = "180px"
			)
		)
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


	# Building chromatogram function
	build_chromato <- function(data, raws, which_chromato, draw_chromato, groups_selected, adjusted, adjusted_time, versus_mode, relative_intensity, color, pos_group, neg_group) {

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

	                if(is.null(raws)) return(displayed_chromatogram)
	                else if(!length(raws)) return(displayed_chromatogram)

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
	                                files_in_group <- length(data@phenoData@data$sample_name[data@phenoData@data$sample_group == groups[j]])
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
	
	                                for ( i in 1:length(raws) ) {
	                                        #Check if file is in group[j]
	                                        if ( data@phenoData@data$sample_group[i] == groups[j] ) {
			                                index <- which(basename(rownames(phenoData(data))) == raws[i])
							# In case of versus mode
			                                if (versus_mode) {
			                                        if (data@phenoData@data$sample_group[index] %in% pos_group) {
					                                intens <- chrom[[index]]@intensity
			                                        } else if (data@phenoData@data$sample_group[index] %in% neg_group) {
					                                intens <- -chrom[[index]]@intensity
			                                        } else {
			                                                intens <- 0
			                                        }
			                                } else {
				                                intens <- chrom[[index]]@intensity
			                                }

							# In case of relative intensity
							if (relative_intensity) {
								if(title=="BPC"){
									# Enlever les blancs ou ne pas les afficher?
									bpi_max <- max(data@featureData@data$basePeakIntensity[data@featureData@data$fileIdx==i])
									intens <- (chrom[[index]]@intensity*100)/bpi_max
								} else {
									# MAX NOT WORK
									intens <- chrom[[index]]@intensity
								}
							}

							# Building the chromatogram
	                                                if ( group_file_nb <= files_to_get ) {
		                                                displayed_chromatogram <- displayed_chromatogram %>% add_lines(
		                                                        x=chrom[[index]]@rtime/60, y=intens, name=basename(raw_files[i]), hoverinfo='text', color=color()[i],
		                                                        text=paste('Intensity: ', round(chrom[[index]]@intensity), '<br />Retention Time: ', round(chrom[[index]]@rtime/60, digits=2))
		                                                )
	                                                } else if ( group_file_nb > (files_in_group - files_to_get) ) {
		                                                displayed_chromatogram <- displayed_chromatogram %>% add_lines(
		                                                        x=chrom[[index]]@rtime/60, y=intens, name=basename(raw_files[i]), hoverinfo='text', color=color()[i],
		                                                        text=paste('Intensity: ', round(chrom[[index]]@intensity), '<br />Retention Time: ', round(chrom[[index]]@rtime/60, digits=2))
		                                                )
	                                                } else {
		                                                displayed_chromatogram <- displayed_chromatogram %>% add_lines(
		                                                        x=chrom[[index]]@rtime/60, y=intens, name=basename(raw_files[i]), hoverinfo='text', color=color()[i], visible="legendonly",
		                                                        text=paste('Intensity: ', round(chrom[[index]]@intensity), '<br />Retention Time: ', round(chrom[[index]]@rtime/60, digits=2))
		                                                )
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

                        if(is.null(raws)) return(displayed_chromatogram)
                        else if(!length(raws)) return(displayed_chromatogram)

			chrom <- bpc_chrom

			# Initialize a variable in case of little group of samples
			files_to_add <- 0

			write(as.character(Sys.time()), file="/import/times.log", append=TRUE)
			write("Parcours de chaque groupe", file="/import/times.log", append=TRUE)


			for ( j in 1:length(groups) ) {

				# Get the samples to display
				files_in_group <- length(data@phenoData@data$sample_name[data@phenoData@data$sample_group == groups[j]])
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

				for ( i in 1:length(raws) ) {
					#Check if file is in group[j]
					if ( data@phenoData@data$sample_group[i] == groups[j] ) {
		                                index <- which(basename(rownames(phenoData(data))) == raws[i])
						intens <- chrom[[index]]@intensity

						# Building the chromatogram
						if ( group_file_nb <= files_to_get ) {
                                                        displayed_chromatogram <- displayed_chromatogram %>% add_lines(
                                                                x=chrom[[index]]@rtime/60, y=intens, name=basename(raw_files[i]), hoverinfo='text', color=data@phenoData@data$sample_name[i],
                                                                text=paste('Intensity: ', round(chrom[[index]]@intensity), '<br />Retention Time: ', round(chrom[[index]]@rtime/60, digits=2))
                                                        )
						} else if ( group_file_nb > (files_in_group - files_to_get) ) {
                                                        displayed_chromatogram <- displayed_chromatogram %>% add_lines(
                                                                x=chrom[[index]]@rtime/60, y=intens, name=basename(raw_files[i]), hoverinfo='text', color=data@phenoData@data$sample_name[i],
                                                                text=paste('Intensity: ', round(chrom[[index]]@intensity), '<br />Retention Time: ', round(chrom[[index]]@rtime/60, digits=2))
                                                        )
						} else {
				                        displayed_chromatogram <- displayed_chromatogram %>% add_lines(
				                                x=chrom[[index]]@rtime/60, y=intens, name=basename(raw_files[i]), hoverinfo='text', color=data@phenoData@data$sample_name[i], visible="legendonly",
				                                text=paste('Intensity: ', round(chrom[[index]]@intensity), '<br />Retention Time: ', round(chrom[[index]]@rtime/60, digits=2))
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

	output$CHROM <- renderPlotly({
		#-----------
		# DEBUG MODE
		write("------------", file="/import/times.log", append=TRUE)
		write("Building Chromatogram", file="/import/times.log", append=TRUE)
		write("------------", file="/import/times.log", append=TRUE)
		#-----------

		displayed_chromatogram <- build_chromato(xdata, raw_files, which_chromatogram(), draw_chromato$value, selected_groups(), adjusted, adjusted_time(), versus_mode(), relative_intensity(), color(), pos_group(), neg_group())

		#-----------
		# DEBUG MODE  
		write("------------", file="/import/times.log", append=TRUE)
		write(as.character(Sys.time()), file="/import/times.log", append=TRUE)
		write("------------", file="/import/times.log", append=TRUE)
		#-----------

		return(displayed_chromatogram)
	})

}

shinyApp(ui, server)
