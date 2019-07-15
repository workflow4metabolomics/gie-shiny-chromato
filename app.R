# Load packages
library(shiny)
library(shinyWidgets)
library(shinyBS)
library(shinydashboard)
library(xcms)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(rlist)
library(Hmisc)
library(tools)
library(rio)

# Get RSession
load("/srv/shiny-server/samples/chromato_visu/inputdata.dat")

## Settings
# Graph settings
height <- "650"

# Init raws variables
raw_names <- xdata@phenoData@data$sample_name
raw_group <- xdata@phenoData@data$sample_group
raw_files <- row.names(xdata@phenoData@data)
groups <- sort(names(table(raw_group)))

# Samples number to display
samples_to_display <- 50
random_number <- 50
if (random_number > length(raw_names)) {
	random_number <- length(raw_names)%/%2
}

# In case of adjusted raws (retcor)
post_retcor <- hasAdjustedRtime(xdata)


## User Interface
ui <- dashboardPage(
	dashboardHeader(
		title="Shiny Chromatogram"
	),
	dashboardSidebar(
		h4(strong("Data Filters")),
		h5(strong("Group displayed :")),
		selectInput(
			inputId = "select_group",
			label = NULL,
			choices = groups,
			selected = groups,
			multiple = TRUE,
			selectize = TRUE,
			width = '200px'
		),
		h5(strong("Randomize :")),
		actionButton(
			inputId = "random_samples",
			label = paste0(random_number," random samples")
		),
		h5(strong("Samples List :")),
		fluidRow(
			column(5,
				style = "margin-left: -15px;",
				actionButton(
					inputId = "check_all",
					label = "Check all"
				)
			),
			column(5,
				actionButton(
					inputId = "uncheck",
					label = "Uncheck all"
				)
			)
		),
		uiOutput("sample_list"),
		column(6,
			style = "margin-left: -15px;",
			actionButton(
				inputId = "export_rdata",
				label = "RData",
				icon = icon("export", lib = "glyphicon")
			)
		),
		bsPopover(
			id = "export_rdata",
			title = "",
			content = "Export to history selected samples in new RData file.",
			placement = "bottom",
			trigger = "hover",
			options = NULL
		),
		column(6,
			actionButton(
				inputId = "export_png",
				label = "PNG",
				icon = icon("export", lib = "glyphicon")
			)
		),
		bsPopover(
			id = "export_png",
			title = "",
			content = "Export the graph in PNG file.",
			placement = "bottom",
			trigger = "hover",
			options = NULL
		)
	),
	dashboardBody(
		includeCSS("styles.css"),
		fluidRow(
			column(11,
				style = "margin: 0px 0px 5px -15px",
				bsCollapse(
					id = "options_panel",
					bsCollapsePanel(
						title = HTML("<h3><b>Graph Options</h3></b>"),
						fluidRow(
							column(3,
								h5(strong("Graph height :")),
								numericInput("height", label=NULL, value=500)
							)
						),
						fluidRow(
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
								h5(strong("Color by Group :")),
								switchInput(
									inputId = "color_by_group",
									label = HTML("<b>Color&nbsp;by&nbsp;Group</b>"),
									value = FALSE
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
							)
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
								style = "margin-left: -15px;",
								fluidRow(
									conditionalPanel(
										condition = "input.versus == true && input.versus_by == false",
										column(6,
											h5(strong("Upper Group :")),
											uiOutput("upper_group")
										)
									),
									conditionalPanel(
										condition = "input.versus == true && input.versus_by == false",
										column(6,
											h5(strong("Under Group :")),
											uiOutput("under_group")
										)
									),
									conditionalPanel(
										condition = "input.versus == true && input.versus_by == true",					
										column(6,
											h5(strong("Upper Sample(s) :")),
											uiOutput("upper_sample")
										)
									),
									conditionalPanel(
										condition = "input.versus == true && input.versus_by == true",
										column(6,
											h5(strong("Under Sample(s) :")),
											uiOutput("under_sample")
										)
									)
								)
							)
						)
					)
				)
			),
			column(1,
				fluidRow(
					actionButton(
						inputId = "draw",
						label = "DRAW",
						class = "btn-primary"
					)
				)
			)
		),
		fluidRow(
			HTML("<small>Zoom: Drag and drop and double-click in the box. Double-click to zoom out.</small>")
		),
		fluidRow(
			fluidRow(
				uiOutput("hover_info",
					style = "position: absolute;"
				),
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
	)
)

server <- function(input, output, session){

	## Initialize variables
	files_list <- reactiveVal(0)
	merged_data <- reactiveVal(0)

	palette <- rainbow(length(raw_names))
	names(palette) <- sort(raw_names)
	reac_palette <- reactiveVal(palette)

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

	# Coloration
	col_group <- eventReactive(input$draw, {
		col_group <- input$color_by_group
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
		selected_samples <- do.call(c, s_list)
	})


	## Notification
	showNotification(
		HTML(paste0("<center><h3><b>BEWARE OF</b></h3></center> <br><br> <center><h4>By default, about ", samples_to_display ," samples are displayed for readability.</h4></center>")),
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
				choices = input$select_group,
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
				choices = sort(files_list()),
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
				choices = sort(subset(files_list(), !(files_list() %in% input$p_sample))),
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
				style = paste0("overflow-y:scroll; background-color: transparent; border-color: transparent; max-height:", as.character(as.numeric(height)-250),"px"),
				lapply(input$select_group, function(group) {
					fluidRow(
						bsCollapse(
							open=h4(paste0(capitalize(group))),
							bsCollapsePanel(
								title = h4(paste0(capitalize(group))),	
								coloredCheckboxGroupInput(
									inputId = paste0("select_",group),
									label = NULL,
									choices = subset(raw_names, raw_group == group),
									selected = subset(raw_names, raw_names %in% files_list()),
									colors = reac_palette()
								)
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

		point <- nearPoints(merged_table, hover, xvar = "rtime", yvar = "intensity", threshold = 10, maxpoints = 1) #addDist=TRUE
		if (nrow(point) == 0) return(NULL)

		left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
		top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
		left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
		top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)

		style <- paste0("position:absolute; width:200px; height:80px; z-index:1000; background-color: rgba(245, 245, 245, 0.85); ","left:", left_px-205, "px; top:", top_px-85, "px;")

		wellPanel(
			style = style,
			p(HTML(paste0(
				"<b>Sample: </b>", point$sample, "<br/>",
				"<b>Retention time: </b>", round(point$rtime,3), "<br/>",
				"<b>Intensity: </b>", point$intensity, "<br/>"
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

	# Update Samples List (Un/Check All, Random)
	observeEvent(input$check_all, {
		lapply(input$select_group, function(group) {
			updateCheckboxGroupInput(
				session = session, 
				inputId = paste0("select_", group),
				selected = raw_names[raw_group == group],
			)
		})
	})
	observeEvent(input$uncheck, {
		lapply(input$select_group, function(group) {
			updateCheckboxGroupInput(
				session = session,
				inputId = paste0("select_", group),
				selected = character(0)
			)
		})
	})
	observeEvent(input$random_samples, {
		random_sample_list <- sample(raw_names[raw_group %in% input$select_group], random_number)
		lapply(input$select_group, function(group) {
			updateCheckboxGroupInput(
				session = session, 
				inputId = paste0("select_", group),
				selected = subset(raw_names, (raw_names %in% random_sample_list))
			)
		})
	})


	## Graph Event
	# Zoom
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
	# Display confirm message box before to export RData in history
	observeEvent(input$export_rdata, {
		confirmSweetAlert(
			session = session,
			inputId = "confirmation",
			title = "Export RData file in history containing only displayed samples.",
			text = HTML("To continue the workflow from the exported RData file,<br>please restart from the step <b>groupChromPeaks/adjustRTime</b>."),
			type = "warning",
			btn_labels = c("Cancel", "Confirm"),
			closeOnClickOutside = TRUE,
			html = TRUE
		)
	})
	# Export filter datas in History as RData file => Change to xdata filter by samples display
	observeEvent(input$confirmation, {
		if (isTRUE(input$confirmation)) {
			basename <- file_path_sans_ext(basename(Sys.getenv("ORIGIN_FILENAME")))
			extension <- ".chromato.RData"
			filename <- paste0(basename, extension, sep="")
			#filename <- "xset.merged.chromato.RData"
			filetype <- "rdata.xcms.findchrompeaks"

			# Update xdata
			samples_ids <- do.call(c,lapply(files_list(), function(file){which(raw_names==file)}))
			xdata <- filterFile(xdata, file=samples_ids, keepAdjustedRtime = FALSE)

			# Update chromatograms
			files_names <- rownames(xdata@phenoData@data)[xdata@phenoData@data$sample_name %in% files_list()]
			chromBPI@.Data <- t(as.matrix(chromBPI@.Data[,colnames(chromBPI@.Data) %in% files_names]))
			chromBPI@phenoData <- xdata@phenoData
			chromTIC@.Data <- t(as.matrix(chromTIC@.Data[,colnames(chromBPI@.Data) %in% files_names]))
			chromTIC@phenoData <- xdata@phenoData
			if (post_retcor) {
				chromBPI_adjusted@.Data <- t(as.matrix(chromBPI_adjusted@.Data[,colnames(chromBPI_adjusted@.Data) %in% files_names]))
				chromBPI_adjusted@phenoData <- xdata@phenoData
				chromTIC_adjusted@.Data <- t(as.matrix(chromTIC_adjusted@.Data[,colnames(chromBPI_adjusted@.Data) %in% files_names]))
				chromTIC_adjusted@phenoData <- xdata@phenoData
			}

			if (exists("zipfile")) { 
				zipfile <- zipfile
			} else {
				zipfile <- NULL
			}

			selected_sample_files <- raw_files[raw_names %in% files_list()]
			singlefile <- subset(singlefile, names(singlefile) %in% basename(selected_sample_files))

			md5sumList$origin <- subset(md5sumList$origin, row.names(md5sumList$origin) %in% selected_sample_files)

			sampleNamesList$sampleNamesOrigin <- files_list()
			sampleNamesList$sampleNamesMakeNames <- files_list()

			# Saving R data in .Rdata file to save the variables used in the present tool
			objects2save = c("xdata","zipfile","singlefile","md5sumList","sampleNamesList", "chromBPI", "chromTIC", "chromBPI_adjusted", "chromTIC_adjusted")
			save(list=objects2save[objects2save %in% ls()], file=filename)

			gx_put(filename, file_type=filetype)
		}
	})

	# Export PNG image in history
	observeEvent(input$export_png, {
		png(filename="plot.png", width=2000, height = 1000)
		plot(displayed_chromatogram() + guides(col=guide_legend(title="Samples")))
		dev.off()
		gx_put("plot.png")
	})


	## Plot
	# Build Dynamic Plot
	displayed_chromatogram <- reactive({
		reactive_chromatogram <- build_chromato(raw_names, raw_group, which_chromatogram(), draw_chromato$value, selected_samples(), post_retcor, adjusted_time(), relative_intensity(), col_group(), versus_mode(), versus_by(), pos_group(), neg_group(), pos_sample(), neg_sample())
	})
	# Render Plot
	output$CHROM <- renderPlot({
		plot <- displayed_chromatogram() + coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE) + guides(col=FALSE)
		return(plot)
	})


	## Functions
	# Building chromatogram function
	build_chromato <- function(raw_names, raw_group, which_chromato, draw_chromato, samples_selected, post_retcor, adjusted_time, relative_intensity, col_group, versus_mode, versus_by, pos_group, neg_group, pos_sample, neg_sample) {

		if(is.null(raw_names)) return(displayed_chromatogram)
		else if(!length(raw_names)) return(displayed_chromatogram)

		if (draw_chromato != 0){

			# According to Adjusted Time
			# If retcor has been done
			if (post_retcor) {
				# If Adjusted Option is TRUE
				if (adjusted_time) {
					if (which_chromato) {
						title <- "Total Ion Chromatogram adjusted"
						chrom <- chromTIC_adjusted
					} else {
						title <- "Base Peak Chromatogram adjusted"
						chrom <- chromBPI_adjusted
					}
				} else {
					if (which_chromato) {
						title <- "Total Ion Chromatogram"
						chrom <- chromTIC
					} else {
						title <- "Base Peak Chromatogram"
						chrom <- chromBPI
					}
				}
			} else {
				if (which_chromato) {
					title <- "Total Ion Chromatogram"
					chrom <- chromTIC
				} else {
					title <- "Base Peak Chromatogram"
					chrom <- chromBPI
				}
			}

			files_list(samples_selected)

			data_tables <- lapply(samples_selected, function(sample){
				index <- which(raw_names==sample)
				sample_group <- raw_group[index]
				intens <- chrom[[index]]@intensity
				intens_max <- max(chrom[[index]]@intensity)

				# In case of versus mode
				if (versus_mode) {
				if (versus_by) {
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
			})
			merged_data <- list.rbind(data_tables)
			merged_data(merged_data)

			# In case of color by group
			if (col_group) {
				colors <- rainbow(length(levels(merged_data[["group"]])))
				palette <- lapply(levels(merged_data[["group"]]), function(group){
					palette <- rep(colors[which(group == levels(merged_data[["group"]]))], length(raw_names[group == raw_group]))
					names(palette) <- sort(raw_names[group == raw_group])
					palette
				})
				palette <- do.call(c, palette)
				reac_palette(palette)
			} else {
				reac_palette(palette)
			}
			
			displayed_chromatogram <- ggplot(data=merged_data, aes(x=c(merged_data[["rtime"]]), y=c(merged_data[["intensity"]]), col=c(as.character(merged_data[["sample"]])))) + scale_colour_manual(name = merged_data[["sample"]], values = reac_palette()) + geom_line() + xlab("Retention Time") + ylab("Intensity") + ggtitle(title) + theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

		} else {

			chrom <- chromBPI

			# Limit the number of samples displayed
			get_samples_list <- NULL
			files_to_add <- 0
			for (group in groups) {
				files_in_group <- raw_names[raw_group == group]
				files_to_get <- samples_to_display%/%length(groups)
				if (length(files_in_group) < files_to_get) {
					files_to_add <- files_to_add+files_to_get-length(files_in_group)
					get_samples_list<- c(get_samples_list, files_in_group)
				} else {
					get_samples_list<- c(get_samples_list, files_in_group[c(1:((files_to_get+files_to_add)/2), (length(files_in_group)-((files_to_get+files_to_add)/2)):length(files_in_group))])
					files_to_add <- 0
				}
			}
			samples_selected <- list.rbind(get_samples_list)
			files_list(samples_selected)

			data_tables <- lapply(samples_selected, function(sample){
				index <- which(raw_names==sample)
				table <- data.frame(rtime=chrom[[index]]@rtime/60, intensity=chrom[[index]]@intensity, sample=as.character(sample))
				names <- rownames(table)
				rownames(table) <- NULL
				table <- as.data.frame(cbind(names, table))
			})
			merged_data <- list.rbind(data_tables)
			merged_data(merged_data)

			displayed_chromatogram <- ggplot(data=merged_data, aes(x=c(merged_data[["rtime"]]), y=c(merged_data[["intensity"]]), col=c(as.character(merged_data[["sample"]])))) + scale_colour_manual(name = merged_data[["sample"]], values = reac_palette()) + geom_line() + xlab("Retention Time") + ylab("Intensity") + ggtitle("Base Peak Chromatogram") + theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
		}
		return(displayed_chromatogram)
	}

	# Checkbox Group Input
	coloredCheckboxGroupInput <- function(inputId, label, choices, selected, colors){
		div(
			id=inputId,
			class="form-group shiny-input-checkboxgroup shiny-input-container shiny-bound-input",
			HTML(paste0('<label class="control-label" for="',inputId,'">',label,'</label>')),
			div(
				style="background-color:transparent;",
				class="shiny-options-group",
				HTML(paste0(
					'<div class="checkbox">',
						'<label>',
							'<input type="checkbox" name="', inputId, '" value="', choices, '"', ifelse(choices %in% selected, 'checked="checked"', ''), '/>',
							'<span ', ifelse(choices %in% selected, paste0('style="font-size: 16px; color:', colors[choices],'"'),'style="font-size: 16px;"'), '>',choices,'</span>',
						'</label>',
					'</div>', collapse = " "
				))
			)
		)
	}
}

shinyApp(ui, server)