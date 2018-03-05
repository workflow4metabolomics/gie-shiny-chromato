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

server <- function(input, output){
	output$uiFiles <- renderUI({
		selectInput('files', 'Select file(s)', choices=setNames(
			c('~/../Documents/GitHub/MaDHaloRS/shiny/mzXMLFiles/negative/20171017_052.mzXML', 
				'~/../Documents/GitHub/MaDHaloRS/shiny/mzXMLFiles/negative/20171017_053.mzXML'), 
			c('20171017_052', '20171017_053')), multiple=TRUE)
	})
	
	output$TIC <- renderPlotly({
		chromato <- plot_ly(source='alignmentChromato', type='scatter', mode='markers') %>% 
			layout(title='TIC', xaxis=list(title='Retention time (min)'), yaxis=list(title='Intensity'), showlegend=TRUE) %>% 
			config(scrollZoom=TRUE, showLink=TRUE, displaylogo=FALSE, 
				modeBarButtons=list(list('toImage', 'zoom2d', 'select2d', 'pan2d', 'autoScale2d', 'resetScale2d')))
		if(is.null(input$files)) return(chromato)
		else if(!length(input$files)) return(chromato)
		raws <- readMSData(input$files, centroided=TRUE, msLevel=1, mode='onDisk')
		rtime <- rtime(raws)
		rtime <- rtime / 60
		#function tic of MsnBase faster than chromatogram of XCMS (for me, to check)
		intensity <- tic(raws)
		#rtime and intensity contains all the data of all files concatenated
		#each value has a name which contain the name of the file and the scan number, separate by a point
		names(rtime) <- names(intensity) <- sapply(names(rtime), function(x) strsplit(x, '\\.')[[1]][1])
		rtime <- split(rtime, names(rtime))
		intensity <- split(intensity, names(intensity))
		for(i in 1:length(rtime)) chromato <- chromato %>% add_lines(
				x=rtime[[i]], y=intensity[[i]], name=basename(input$files[i]), hoverinfo='text', 
			text=~paste('Intensity: ', round(intensity[[i]]), '<br />Retention Time: ', 
				round(rtime[[i]], digits=2)))
		return(chromato)
	})
	
	output$BPC <- renderPlotly({
		chromato <- plot_ly(source='alignmentChromato', type='scatter', mode='markers') %>% 
			layout(title='BPC', xaxis=list(title='Retention time (min)'), yaxis=list(title='Intensity'), showlegend=TRUE) %>% 
			config(scrollZoom=TRUE, showLink=TRUE, displaylogo=FALSE, 
				modeBarButtons=list(list('toImage', 'zoom2d', 'select2d', 'pan2d', 'autoScale2d', 'resetScale2d')))
		if(is.null(input$files)) return(chromato)
		else if(!length(input$files)) return(chromato)
		raws <- readMSData(input$files, centroided=TRUE, msLevel=1, mode='onDisk')
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