#load packages
library(shiny)
library(xcms)
library(plotly)
library(stringr)

ui <- bootstrapPage(
	checkboxInput(
		inputId = "adjustedTime",
		label = strong("Adjusted Time"),
		value = FALSE
	),
	plotlyOutput('TIC'),
	plotlyOutput('BPC')
)

load("/srv/shiny-server/data/inputdata.dat")
raw_files <- basename(rownames(xdata@phenoData@data))

# identifier_type name because of filename and not hid
gx_get(raw_files, identifier_type='name')

for (file in raw_files){
  system(sprintf("ln -s /import/%s %s", file, file))
}

server <- function(input, output){
	
	output$TIC <- renderPlotly({
		chromato <- plot_ly(source='alignmentChromato', type='scatter', mode='markers') %>% 
			layout(title='TIC', xaxis=list(title='Retention time (min)'), yaxis=list(title='Intensity'), showlegend=TRUE) %>% 
			config(scrollZoom=TRUE, showLink=TRUE, displaylogo=FALSE, 
				modeBarButtons=list(list('toImage', 'zoom2d', 'select2d', 'pan2d', 'autoScale2d', 'resetScale2d')))
		
		if(is.null(raw_files)) return(chromato)
		else if(!length(raw_files)) return(chromato)
		
		raws <- xdata

		#function tic of MsnBase faster than chromatogram of XCMS (for me, to check)
		rtime <- split(rtime(raws)/60, f = fromFile(raws))
		intensity <- split(tic(raws), f = fromFile(raws))

		for(i in 1:length(raw_files)) {
			index <- which(basename(rownames(phenoData(raws))) == raw_files[i])
			chromato <- chromato %>% add_lines(
				x=rtime[[index]], y=intensity[[index]], name=basename(raw_files[i]), hoverinfo='text', 
				text=~paste('Intensity: ', round(intensity[[index]]), '<br />Retention Time: ', round(rtime[[index]], digits=2))
			)
		}

		return(chromato)
	})
	
	output$BPC <- renderPlotly({
		chromato <- plot_ly(source='alignmentChromato', type='scatter', mode='markers') %>% 
			layout(title='BPC', xaxis=list(title='Retention time (min)'), yaxis=list(title='Intensity'), showlegend=TRUE) %>% 
			config(scrollZoom=TRUE, showLink=TRUE, displaylogo=FALSE, 
				modeBarButtons=list(list('toImage', 'zoom2d', 'select2d', 'pan2d', 'autoScale2d', 'resetScale2d')))
		if(is.null(raw_files)) return(chromato)
		else if(!length(raw_files)) return(chromato)

		raws <- xdata

		# With adjusted
		#raws <- adjustRtime(raws, param = ObiwarpParam(binSize = 0.6))

                #function tic of MsnBase faster than chromatogram of XCMS (for me, to check)
		rtime <- split(rtime(raws)/60, f = fromFile(raws))
                intensity <- split(bpi(raws, initial = FALSE), f = fromFile(raws))

                for(i in 1:length(raw_files)) {
                        index <- which(basename(rownames(phenoData(raws))) == raw_files[i])
                        chromato <- chromato %>% add_lines(
                                x=rtime[[index]], y=intensity[[index]], name=basename(raw_files[i]), hoverinfo='text',
                                text=~paste('Intensity: ', round(intensity[[index]]), '<br />Retention Time: ', round(rtime[[index]], digits=2))
                        )
                }

		return(chromato)
	})
}

shinyApp(ui, server)
