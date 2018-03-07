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

		#function tic of MsnBase faster than chromatogram of XCMS (for me, to check)
		rtime <- split(rtime(raws)/60, f = fromFile(raws))
		intensity <- split(tic(raws), f = fromFile(raws))

		for(i in 1:length(input$files)) {
			index <- which(basename(rownames(phenoData(raws))) == input$files[i])
			chromato <- chromato %>% add_lines(
				x=rtime[[index]], y=intensity[[index]], name=basename(input$files[index]), hoverinfo='text', 
				text=~paste('Intensity: ', round(intensity[[index]]), '<br />Retention Time: ', round(rtime[[index]], digits=2))
			)
		}

		#rtimes <- chromatogram(raws, aggregationFun = 'sum')
		#rtime <- rtimes[ ,basename(colnames(rtimes)) %in% input$files]
		#for(i in 1:length(input$files)) {
		#	if (length(input$files)==1) {
                #                chromato <- chromato %>% add_lines(
                #                        x=rtime@rtime/60, y=rtime@intensity, name=basename(input$files), hoverinfo='text',
                #                        text=~paste('Intensity: ', round(rtime@intensity), '<br />Retention Time: ', round(rtime@rtime/60, digits=2))
		#		)
		#	} else {
		#		chromato <- chromato %>% add_lines(
		#			x=rtime[[i]]@rtime/60, y=rtime[[i]]@intensity, name=basename(input$files[i]), hoverinfo='text', 
		#			text=~paste('Intensity: ', round(rtime[[i]]@intensity), '<br />Retention Time: ', round(rtime[[i]]@rtime/60, digits=2))
		#		)
		#	}
		#}

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

		#points is a list of bpc chromatogram object
		points <- chromatogram(raws, aggregationFun = 'max')
		point <- points[ ,basename(colnames(points)) %in% input$files]

		#rtime <- split(rtime(raws)/60, f = fromFile(raws))
                #function tic of MsnBase faster than chromatogram of XCMS (for me, to check)
                #intensity <- split(bpc(raws), f = fromFile(raws))

                for(i in 1:length(input$files)) {
                        index <- which(basename(rownames(phenoData(raws))) == input$files[i])
                        chromato <- chromato %>% add_lines(
                                x=rtime[[index]], y=intensity[[index]], name=basename(input$files[index]), hoverinfo='text',
                                text=~paste('Intensity: ', round(intensity[[index]]), '<br />Retention Time: ', round(rtime[[index]], digits=2))
                        )
                }

		#for(i in 1:length(input$files)) {
		#	if (length(input$files)==1) {
		#		chromato <- chromato %>% add_lines(
		#			x=point@rtime/60, y=point@intensity, name=basename(input$files), hoverinfo='text', 
		#			text=~paste('Intensity: ', round(point@intensity), '<br />Retention Time: ', 
		#			round(point@rtime/60))
		#		)
		#	} else {
		#                 chromato <- chromato %>% add_lines(
                #                        x=point[[i]]@rtime/60, y=point[[i]]@intensity, name=basename(input$files[i]), hoverinfo='text',
	        #                        text=~paste('Intensity: ', round(point[[i]]@intensity), '<br />Retention Time: ',
		#			round(point[[i]]@rtime/60))
		#		)
		#	}
		#}

		return(chromato)
	})
}

shinyApp(ui, server)
