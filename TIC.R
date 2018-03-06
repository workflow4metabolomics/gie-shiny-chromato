#load packages
library(shiny)
library(xcms)
library(plotly)
library(stringr)

ui <- fluidPage(
	sidebarLayout(
		sidebarPanel(
			uiOutput('uiFiles')
		),
		mainPanel(
			plotlyOutput('TIC'),
			plotlyOutput('BPC')
		)
	)
)

load("/srv/shiny-server/data/inputdata.dat")
raw_files <- basename(rownames(xdata@phenoData@data))

# identifier_type name because of filename and not hid
gx_get(raw_files, identifier_type='name')

for (file in raw_files){
  system(sprintf("ln -s /import/%s %s", file, file))
}

server <- function(input, output){
	output$uiFiles <- renderUI({
		selectInput('files', 'Select file(s)', choices=setNames(raw_files, raw_files), multiple=TRUE)
	})
	
	output$TIC <- renderPlotly({
		chromato <- plot_ly(source='alignmentChromato', type='scatter', mode='markers') %>% 
			layout(title='TIC', xaxis=list(title='Retention time (min)'), yaxis=list(title='Intensity'), showlegend=TRUE) %>% 
			config(scrollZoom=TRUE, showLink=TRUE, displaylogo=FALSE, 
				modeBarButtons=list(list('toImage', 'zoom2d', 'select2d', 'pan2d', 'autoScale2d', 'resetScale2d')))
		
		if(is.null(input$files)) return(chromato)
		else if(!length(input$files)) return(chromato)
		
		#raws <- readMSData(input$files, centroided=TRUE, msLevel=1, mode='onDisk')
		raws <- xdata
		#raws <- xdata[basename(rownames(xdata@phenoData@data)) %in% input$files, ]
		#rtime <- rtime(raws)
		#rtime <- rtime / 60
		#function tic of MsnBase faster than chromatogram of XCMS (for me, to check)
		#intensity <- tic(raws)
		#rtime and intensity contains all the data of all files concatenated
		
		#each value has a name which contain the name of the file and the scan number, separate by a point
		#names(rtime) <- names(intensity) <- sapply(names(rtime), function(x) strsplit(x, '\\.')[[1]][1])
		#rtime <- split(rtime, names(rtime))
		#intensity <- split(intensity, names(intensity))
		rtimes <- chromatogram(raws, aggregationFun = 'sum')
		rtime <- rtimes[ ,basename(colnames(rtimes)) %in% input$files]
		for(i in 1:length(rtime)) {
			chromato <- chromato %>% add_lines(
				x=rtime[[i]]@rtime/60, y=rtime[[i]]@intensity, name=basename(input$files[i]), hoverinfo='text', 
				text=~paste('Intensity: ', round(rtime[[i]]@intensity), '<br />Retention Time: ', 
				round(rtime[[i]]@rtime/60, digits=2))
			)
		}
		return(chromato)
	})
	
	output$BPC <- renderPlotly({
		chromato <- plot_ly(source='alignmentChromato', type='scatter', mode='markers') %>% 
			layout(title='BPC', xaxis=list(title='Retention time (min)'), yaxis=list(title='Intensity'), showlegend=TRUE) %>% 
			config(scrollZoom=TRUE, showLink=TRUE, displaylogo=FALSE, 
				modeBarButtons=list(list('toImage', 'zoom2d', 'select2d', 'pan2d', 'autoScale2d', 'resetScale2d')))
		if(is.null(input$files)) return(chromato)
		else if(!length(input$files)) return(chromato)
		#raws <- readMSData(input$files, centroided=TRUE, msLevel=1, mode='onDisk')
		raws <- xdata
		points <- chromatogram(raws, aggregationFun = 'max')
		#points is a list of chromatogram object
		for(i in 1:length(points)) chromato <- chromato %>% add_lines(
				x=points[[i]]@rtime/60, y=points[[i]]@intensity, name=basename(input$files[i]), hoverinfo='text', 
			text=~paste('Intensity: ', round(points[[i]]@intensity), '<br />Retention Time: ', 
				round(points[[i]]@rtime/60)))
		return(chromato)
	})
}

shinyApp(ui, server)
