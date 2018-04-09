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
raw_groups <- xdata@phenoData@data$sample_group
groups <- unique(raw_groups)

## Import files by copying them, not used because to slow (identifier_type='name' because of filename and not hid)
#gx_get(raw_files, identifier_type='name')
## Symbolic link between /import and APP_PATH
#for (file in raw_files){
#  system(sprintf("ln -s /srv/shiny-server/data/datasets/%s %s", file, file))
#}

# In case of adjusted raws (retcor)
adjusted <- hasAdjustedRtime(xdata)

# Making a color palette
default_palette <- rainbow(length(raw_files))

# Samples number to display
samples_to_display <- 50

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
)

server <- function(input, output){

	# Get pre-calculate chromatogram from chromTIC and chromBPI objects
	tic_chrom <- chromTIC
	bpc_chrom <- chromBPI
	if (adjusted) {
		tic_chrom_adjusted <- chromTIC_adjusted
		bpc_chrom_adjusted <- chromBPI_adjusted
	}

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
	build_chromato <- function(data, raws, title, chromatogram, chromatogram_adjusted, draw_chromato, adjusted, adjusted_time, versus_mode, color, pos_group, neg_group) {

                if (draw_chromato != 0){

	                displayed_chromatogram <- plot_ly(source='alignmentChromato', type='scatter', mode='markers', colors=palette()) %>%
	                        layout(title=title, xaxis=list(title='Retention time (min)'), yaxis=list(title='Intensity'), showlegend=TRUE) %>%
	                        config(scrollZoom=TRUE, showLink=TRUE, displaylogo=FALSE,
	                                modeBarButtons=list(list('toImage', 'zoom2d', 'select2d', 'pan2d', 'autoScale2d', 'resetScale2d')))

	                if(is.null(raws)) return(displayed_chromatogram)
	                else if(!length(raws)) return(displayed_chromatogram)

	                # According to Adjusted Time
	                if (adjusted) {
	                         if (adjusted_time) {
	                                chrom <- chromatogram_adjusted
	                         } else {
					chrom <- chromatogram
				}
	                } else {
	                         chrom <- chromatogram
	                }

                        for ( j in 1:length(groups) ) {
                                # Nb de fichiers dans le groupe
                                files_in_group <- length(data@phenoData@data$sample_name[data@phenoData@data$sample_group == groups[j]])

                                files_to_get <- (samples_to_display%/%(length(groups)*2))

                                if ( files_in_group < (files_to_get*2) ) {
                                        files_to_get <- files_in_group
                                }

                                #Initialisation d'un compteur
                                group_file_nb <- 1

                                for ( i in 1:length(raws) ) {
                                        #Check if file is in group[j]
                                        if ( data@phenoData@data$sample_group[i] == groups[j] ) {
		                                index <- which(basename(rownames(phenoData(data))) == raws[i])
		                                if (versus_mode) {
		                                        if (data@phenoData@data$sample_group[index] %in% pos_group) {
				                                intens = chrom[[index]]@intensity
		                                        } else if (data@phenoData@data$sample_group[index] %in% neg_group) {
				                                intens = -chrom[[index]]@intensity
		                                        } else {
		                                                intens = 0
		                                        }
		                                } else {
			                                intens = chrom[[index]]@intensity
		                                }

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

                                #samples_to_display <- samples_to_display - files_to_get
                                #files_to_get <- (samples_to_display%/%(length(groups)*2))

                                #-----------
                                # DEBUG MODE
                                write("Groupe :", file="/import/times.log", append=TRUE)
                                write(groups[j], file="/import/times.log", append=TRUE)
                                write("Nb de fichiers dans le groupe :", file="/import/times.log", append=TRUE)
                                write(files_in_group, file="/import/times.log", append=TRUE)
                                write("Nb de fichiers à récupérer dans le groupe :", file="/import/times.log", append=TRUE)
                                write(files_to_get, file="/import/times.log", append=TRUE)
                                write("The Final Countdown :", file="/import/times.log", append=TRUE)
                                write(group_file_nb-1, file="/import/times.log", append=TRUE)
                                #-----------

                        }

                } else {

                        displayed_chromatogram <- plot_ly(source='alignmentChromato', type='scatter', mode='markers', colors=default_palette) %>%
                                layout(title=title, xaxis=list(title='Retention time (min)'), yaxis=list(title='Intensity'), showlegend=TRUE) %>%
                                config(scrollZoom=TRUE, showLink=TRUE, displaylogo=FALSE,
                                        modeBarButtons=list(list('toImage', 'zoom2d', 'select2d', 'pan2d', 'autoScale2d', 'resetScale2d')))

                        if(is.null(raws)) return(displayed_chromatogram)
                        else if(!length(raws)) return(displayed_chromatogram)

			chrom <- chromatogram

			for ( j in 1:length(groups) ) {
				# Nb de fichiers dans le groupe
				files_in_group <- length(data@phenoData@data$sample_name[data@phenoData@data$sample_group == groups[j]])

				files_to_get <- (samples_to_display%/%(length(groups)*2))
				
				if ( files_in_group < (files_to_get*2) ) {
					files_to_get <- files_in_group
				}

				#Initialisation d'un compteur
				group_file_nb <- 1

				for ( i in 1:length(raws) ) {
					#Check if file is in group[j]
					if ( data@phenoData@data$sample_group[i] == groups[j] ) {
		                                index <- which(basename(rownames(phenoData(data))) == raws[i])
						intens = chrom[[index]]@intensity

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

				#samples_to_display <- samples_to_display - files_to_get
				#files_to_get <- (samples_to_display%/%(length(groups)*2))

		                #-----------
		                # DEBUG MODE
                                write("Groupe :", file="/import/times.log", append=TRUE)
		                write(groups[j], file="/import/times.log", append=TRUE)
                                write("Nb de fichiers dans le groupe :", file="/import/times.log", append=TRUE)
                                write(files_in_group, file="/import/times.log", append=TRUE)
                                write("Nb de fichiers à récupérer dans le groupe :", file="/import/times.log", append=TRUE)
                                write(files_to_get, file="/import/times.log", append=TRUE)
                                write("The Final Countdown :", file="/import/times.log", append=TRUE)
		                write(group_file_nb-1, file="/import/times.log", append=TRUE)
		                #-----------

			}
               }

		return(displayed_chromatogram)
	}

	output$TIC <- renderPlotly({
		#-----------
		# DEBUG MODE
		write("Building TIC", file="/import/times.log", append=TRUE)
		#-----------

		displayed_chromatogram <- build_chromato(xdata, raw_files, "TIC", tic_chrom, tic_chrom_adjusted, draw_chromato$value, adjusted, adjusted_time(), versus_mode(), color(), pos_group(), neg_group())
		
		#-----------
		# DEBUG MODE  
		write(as.character(Sys.time()), file="/import/times.log", append=TRUE)
                write("Done", file="/import/times.log", append=TRUE)  
		#-----------

		return(displayed_chromatogram)
	})

	output$BPC <- renderPlotly({		
		#-----------
                # DEBUG MODE
                write("Building BPC", file="/import/times.log", append=TRUE)
		#-----------

                displayed_chromatogram <- build_chromato(xdata, raw_files, "BPC", bpc_chrom, bpc_chrom_adjusted, draw_chromato$value, adjusted, adjusted_time(), versus_mode(), color(), pos_group(), neg_group())

                #-----------
		# DEBUG MODE
                write(as.character(Sys.time()), file="/import/times.log", append=TRUE)
                write("Done", file="/import/times.log", append=TRUE)
		#-----------

		return(displayed_chromatogram)
	})
}

shinyApp(ui, server)
