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
library(shinyBS)
library(xcms)
library(RColorBrewer)
library(stringr)
library(webshot)
library(tools)
library(rlist)
library(Hmisc)

# Get RSession
load("/srv/shiny-server/samples/chromato_visu/inputdata.dat")
raw_names <- xdata@phenoData@data$sample_name
raw_group <- xdata@phenoData@data$sample_group
raw_files <- row.names(xdata@phenoData@data)


## Settings
# Graph settings
height <- "600"

# Samples number to display
samples_to_display <- 200

# Get group names sorted
groups <- sort(names(table(raw_group)))

# In case of adjusted raws (retcor)
post_retcor <- hasAdjustedRtime(xdata)

# Making a color palette
default_palette <- rainbow(length(raw_names))

## Import files by copying them, not used because to slow (identifier_type='name' because of filename and not hid)
#gx_get(raw_names, identifier_type='name')

#----------
#DEBUG MODE
write(as.character(Sys.time()), file="/import/times.log", append=TRUE)
#----------


## User Interface
ui <- bootstrapPage(
	fluidRow(
		includeCSS("styles.css"),
		column(12,
			wellPanel(
				HTML("<h3><b>Filters</h3></b>"),
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
						if (post_retcor) {
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
				)
			)
		)
	),
	br(),
	fluidRow(
		style = paste("height:", height, "px", sep=""), 
		column(2,
			style = "height : 100%",
			uiOutput("sample_list")
		),
		column(10,
			fluidRow(
				uiOutput("hover_info")
			),
			fluidRow(
				plotOutput('CHROM', 
					height = height,
					dblclick = "dblclick",
			        brush = brushOpts(
			        	id = "brush",
			        	resetOnNew = TRUE
			        ),
		        	hover = hoverOpts("plot_hover", delay = 100, delayType = "debounce")
		        )
			)
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

server <- function(input, output, session){

	## Initialize variables
	# Get pre-calculate chromatogram from chromTIC and chromBPI objects
	chrom_tic <- chromTIC
	chrom_bpi <- chromBPI
	if (post_retcor) {
		chrom_tic_adjusted <- chromTIC_adjusted
		chrom_bpi_adjusted <- chromBPI_adjusted
	}

	files_list <- reactiveVal(0)
	merged_data <- reactiveVal(0)


	## Get Filters Values
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

	# Groups displayed
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

	# Versus
	versus_mode <- eventReactive(input$draw, {
		versus_mode <- input$versus
	})
	versus_by <- eventReactive(input$draw, {
		versus_by <- input$versus_by
	})
	#-groups-#
	pos_group <- eventReactive(input$draw, {
		pos_group <- input$p_group
	})
	neg_group <- eventReactive(input$draw, {
		neg_group <- input$n_group
	})
	#-samples-#
	pos_sample <- eventReactive(input$draw, {
		pos_sample <- input$p_sample
	})
	neg_sample <- eventReactive(input$draw, {
		neg_sample <- input$n_sample
	})

    # Samples displayed
	selected_samples <- eventReactive(input$draw, {
		s_list <- lapply(input$select_group, function(group){
			input[[paste0("select_",group)]]
		})
		selected_samples <- sort(unique(c(list.rbind(s_list))))
	})


	## Notification
	showNotification(
		HTML("<center><h3><b>BE CAREFUL.</b></h3></center> <br><br> <center><h4>By default, only few samples are displayed.</h4></center>"),
		duration = 20,
		closeButton = TRUE,
	  	id = NULL,
	  	type = "warning",
		session = getDefaultReactiveDomain()
	)


	## Dynamic User Interface
	# Adjusted filter (if retcor done)
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

	# Versus Upper/Under Filters
	#-groups-#
	output$upper_group <- renderUI({
		tagList(
			selectInput(
				inputId = "p_group",
				label = NULL,
				choices = subset(groups, (groups %in% input$select_group)),
				width = "180px"
			)	
		)
	})
	output$under_group <- renderUI({
		tagList(
			selectInput(
				inputId = "n_group",
				label = NULL,
				choices = subset(input$select_group, !(input$select_group %in% input$p_group)),
				width = "180px"
			)
		)
	})
	#-samples-#
	# TODO : Selectize TRUE with limited nb of sample or FALSE (ctrl + choice) ...
	output$upper_sample <- renderUI({
		tagList(
			selectInput(
				inputId = "p_sample",
				label = NULL,
				choices = sort(subset(raw_names, (raw_group %in% input$select_group) & (raw_names %in% files_list()) )),
				multiple = TRUE,
				selectize = FALSE,
				width = "180px"
			)	
		)
	})
	output$under_sample <- renderUI({
		tagList(
			selectInput(
				inputId = "n_sample",
				label = NULL,
				choices = sort(subset(raw_names, (raw_group %in% input$select_group) & (raw_names %in% files_list()) & !(raw_names %in% input$p_sample))),
				multiple = TRUE,
				selectize = FALSE,
				width = "180px"
			)
		)
	})

	# Samples List
	output$sample_list <- renderUI({
		tagList(
			wellPanel(
				id = "tPanel",
				style = "overflow-y:scroll; max-height: 600px",
				HTML("<h3><b>Samples List</h3></b>"),
				hr(),
				actionButton(
					inputId = "button_all",
					label = "Select all"
				),
				actionButton(
					inputId = "button_no",
					label = "Unselect all"
				),
				br(),br(),
			    lapply(input$select_group, function(group) {
			    	fluidRow(
			    		bsCollapsePanel(
							title = h4(paste0(capitalize(group))),					
							checkboxGroupInput(
			        			inputId = paste0("select_",group),
			        			label = NULL,
			                	choices = sort(subset(raw_names, (raw_group %in% group) & (group==raw_group))),
			                	selected = sort(subset(raw_names, (raw_group %in% group) & (raw_names %in% files_list()) )),
			                	width = "200px"
			            	)
			            )
			      	)
				})
			)
		)
	})

	# Hover Panel
	output$hover_info <- renderUI({

		merged_table<-merged_data()
		hover <- input$plot_hover

	    point <- nearPoints(merged_table, hover, xvar = "rtime", yvar = "intensity", threshold = 10, maxpoints = 1)
	    if (nrow(point) == 0) return(NULL)

	    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
	    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    	left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
	    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
	    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ","left:", left_px + 2, "px; top:", top_px + 2, "px;")

        wellPanel(
	    	style = style,
		    p(HTML(paste0(
		    	"<b>Sample: </b>", point$sample, "<br/>",
		        "<b>Retention time: </b>", point$rtime, "<br/>",
		        "<b>Intensity: </b>", point$intensity
		    )))
        )
	})


	## Update Events
	# Update Sample Versus
	observeEvent(input$draw, {
    	updateSelectInput(
    		session = session, 
    		inputId = "p_sample",
    		choices = selected_samples()
    	)
  	})
	observeEvent(input$draw, {
    	updateSelectInput(
    		session = session, 
    		inputId = "n_sample",
    		choices = subset(selected_samples(), !(selected_samples() %in% input$p_sample)),
    	)
  	})

	# Update Samples List (Un/Select All)
	observeEvent(input$button_all, {
		lapply(input$select_group, function(group) {
			updateCheckboxGroupInput(
				session = session, 
				inputId = paste0("select_", group),
				choices = sort(subset(raw_names, (raw_group %in% group) & (group==raw_group))),
				selected = sort(subset(raw_names, (raw_group %in% group) & (group==raw_group)))
			)
		})
	})
	observeEvent(input$button_no, {
		lapply(input$select_group, function(group) {
			updateCheckboxGroupInput(
				session = session,
				inputId = paste0("select_", group),
				choices = sort(subset(raw_names, (raw_group %in% group) & (group==raw_group))),
				selected = NULL
			)
		})
	})

	## Graph Event
	# Hover
	ranges <- reactiveValues(x = NULL, y = NULL)
	observeEvent(input$dblclick, {
	    brush <- input$brush
	    if (!is.null(brush)) {
	    	ranges$x <- c(brush$xmin, brush$xmax)
	    	ranges$y <- c(brush$ymin, brush$ymax)
	    } else {
	    	ranges$x <- NULL
	    	ranges$y <- NULL
	    }
  	})

	## Export Event
	# Export filter datas in History as RData file => Change to xdata filter by samples display
	observeEvent(input$export, {

		basename <- file_path_sans_ext(Sys.getenv("ORIGIN_FILENAME"))
		extension <- ".chromato.RData"
		filename <- paste0(basename, extension, sep="")
		filetype <- "rdata.xcms.findchrompeaks"

		LL<-c()
		for(i in c(1:length(input$select_sample))){
			raw_id<-which(raw_names==input$select_sample[i])
			LL<-c(LL,raw_id)
		}
		xdata<-filterFile(xdata, file=LL, keepAdjustedRtime = FALSE)
		chromBPI <- chromBPI
		chromTIC <- chromTIC
		if (post_retcor){
			chromBPI_adjusted <- chromBPI_adjusted
			chromTIC_adjusted <- chromTIC_adjusted
		}

		if (exists("zipfile")) { 
			zipfile <- zipfile
		} else {
			zipfile <- NULL
		}

		selected_sample_files<-raw_files[raw_names %in% input$select_sample]
		singlefile <- subset(singlefile, names(singlefile) %in% basename(selected_sample_files))

		md5sumList$origin<-subset(md5sumList$origin, row.names(md5sumList$origin) %in% selected_sample_files)

		sampleNamesList$sampleNamesOrigin <- input$select_sample
		sampleNamesList$sampleNamesMakeNames <- input$select_sample

		#saving R data in .Rdata file to save the variables used in the present tool
		objects2save = c("xdata","zipfile","singlefile","md5sumList","sampleNamesList", "chromBPI", "chromTIC", "chromBPI_adjusted", "chromTIC_adjusted")
		save(list=objects2save[objects2save %in% ls()], file=filename)

		gx_put(filename, file_type=filetype)
	})

	# Export PNG image in history
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


	## Plot
	# Build Dynamic Plot
	displayed_chromatogram <- reactive({
		reactive_chromatogram <- build_chromato(raw_names, raw_group, which_chromatogram(), draw_chromato$value, selected_groups(), selected_samples(), post_retcor, adjusted_time(), relative_intensity(), col_group(), versus_mode(), versus_by(), pos_group(), neg_group(), pos_sample(), neg_sample())
	})		
	# Render Plot
	output$CHROM <- renderPlot({
		plot <- displayed_chromatogram() + coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
		return(plot)
	})


	## Functions
	# Building chromatogram function
	build_chromato <- function(raw_names, raw_group, which_chromato, draw_chromato, groups_selected, samples_selected, post_retcor, adjusted_time, relative_intensity, col_group, versus_mode, versus_by, pos_group, neg_group, pos_sample, neg_sample) {

        if (draw_chromato != 0){

	        if(is.null(raw_names)) return(displayed_chromatogram)
	        else if(!length(raw_names)) return(displayed_chromatogram)

	        # According to Adjusted Time
			# If retcor has been done
	        if (post_retcor) {
				# If Adjusted Option is TRUE
	            if (adjusted_time) {
					if (which_chromato) {
						title <- "TIC adjusted"
		                chrom <- chrom_tic_adjusted
					} else {
						title <- "BPC adjusted"
		                chrom <- chrom_bpi_adjusted
					}
	            } else {
					if (which_chromato) {
						title <- "TIC"
		                chrom <- chrom_tic
					} else {
						title <- "BPC"
		                chrom <- chrom_bpi
					}
				}
	        } else {
				if (which_chromato) {
					title <- "TIC"
	                chrom <- chrom_tic
				} else {
					title <- "BPC"
	                chrom <- chrom_bpi
				}
	        }

	        # Initialize a variable in case of little group of samples
			get_files_list <- lapply(groups, function(group){
				if (group %in% groups_selected){
					files_in_group <- raw_names[raw_group == group & raw_names %in% samples_selected]
				}
			})
			files_list <- list.rbind(get_files_list)
			files_list(files_list)

  			sample_table <- NULL
			data_tables <- lapply(files_list, function(sample){
				if (sample %in% samples_selected){
					index <- which(raw_names==sample)
					sample_group <- raw_group[index]
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

					table <- data.frame(rtime=chrom[[index]]@rtime/60, intensity=intens, sample=as.character(sample), group=as.character(sample_group))
					names <- rownames(table)
					rownames(table) <- NULL
					table <- as.data.frame(cbind(names, table))
					sample_table <- rbind(sample_table,table)
				}
			})
			merged_data <- list.rbind(data_tables)
			merged_data(merged_data)

			# In case of color by group
			if(col_group){
				displayed_chromatogram <- ggplot(data=merged_data, aes(x=c(merged_data[["rtime"]]), y=c(merged_data[["intensity"]]), group=c(as.character(merged_data[["sample"]])), col=c(as.character(merged_data[["group"]])))) + geom_line() + xlab("Retention Time") + ylab("Intensity") + guides(col=guide_legend(title="Samples")) #+ guides(col=FALSE)
			} else {
				displayed_chromatogram <- ggplot(data=merged_data, aes(x=c(merged_data[["rtime"]]), y=c(merged_data[["intensity"]]), col=c(as.character(merged_data[["sample"]])))) + geom_line() + xlab("Retention Time") + ylab("Intensity") + guides(col=guide_legend(title="Samples")) #+ guides(col=FALSE)
			}

        } else {

            if(is.null(raw_names)) return(displayed_chromatogram)
            else if(!length(raw_names)) return(displayed_chromatogram)

			chrom <- chrom_bpi

			# Initialize a variable in case of little group of samples
			get_files_list <- lapply(groups, function(group){
				files_in_group <- raw_names[raw_group == group]
				files_to_get <- samples_to_display%/%length(groups)
				if ( length(files_in_group) < files_to_get ) {
					files_in_group
				} else {
					files_in_group[c(1:(files_to_get/2), (length(files_in_group)-(files_to_get/2)):length(files_in_group))]
				}
			})
			files_list <- list.rbind(get_files_list)
			files_list(files_list)

			sample_table <- NULL
			data_tables <- lapply(files_list, function(sample){
				index <- which(raw_names==sample)
				table <- data.frame(rtime=chrom[[index]]@rtime/60, intensity=chrom[[index]]@intensity, sample=as.character(sample))
				names <- rownames(table)
				rownames(table) <- NULL
				table <- as.data.frame(cbind(names, table))
				sample_table <- rbind(sample_table,table)
			})
			merged_data <- list.rbind(data_tables)
			merged_data(merged_data)

			displayed_chromatogram <- ggplot(data=merged_data, aes(x=c(merged_data[["rtime"]]), y=c(merged_data[["intensity"]]), col=c(as.character(merged_data[["sample"]])))) + geom_line() + xlab("Retention Time") + ylab("Intensity") + guides(col=guide_legend(title="Samples")) + theme(legend.box.margin = unit(c(0,0,0,0), units="mm")) #+ guides(col=FALSE)
        }
		return(displayed_chromatogram)
	}
}

shinyApp(ui, server)