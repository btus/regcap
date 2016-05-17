# server.R

library(xts); library(dygraphs); library(data.table); library(lubridate)

#increase max file upload to 160 mb
options(shiny.maxRequestSize=160*1024^2) 

#colSelct takes a data object and a flexible number of column indexes, and returns a new data object including just those columns selected by the user. The dygraph plotting utility requires an entire data object, rather than selecting columns from existing object.

colSelect=function(dat_obj, col1, ...){
	return(dat_obj[,c(col1, ...)])
}

shinyServer(

	function(input, output) {
		
 		output$list_of_files <- renderUI({
 			file_list<-input$file1[,'name']
 			checkboxGroupInput("path", "Choose File", file_list)
		})
		
		output$list_of_variables <- renderUI({
			filepath<-input$file1[input$file1[,'name']==input$path,'datapath']
			#variables_list<-colnames(fread(filepath, header=TRUE, sep='\t'))
			variables_list<-colnames(data.xts())
			checkboxGroupInput("colSelection", "Columns to Plot", variables_list)
		})
 		
 		data.xts<-reactive({
 			filepath<-input$file1[input$file1[,'name']==input$path,'datapath']
 			if(substr(input$path, nchar(input$path)-3, nchar(input$path)) == ".rco"){
 				data<-fread(filepath, header=TRUE, sep="\t")
 			} else{
 				data <- fread(filepath)
 				data <- data[,Index := NULL]
 				}
 			 			
 			if(nrow(data) == 8760){
 				DateTime <- seq(ymd_hm("2015-01-01 00:00", tz='GMT'), ymd_hm("2015-12-31 23:00", tz='GMT'), "hours")
 			} else if(nrow(data) == 525600){
 				DateTime <- seq(ymd_hm("2015-01-01 00:00", tz='GMT'), ymd_hm("2015-12-31 23:59", tz='GMT'), "min")
 			} else if(nrow(data) == 365){
 				DateTime <- seq(ymd_hm("2015-01-01 00:00", tz='GMT'), ymd_hm("2015-12-31 23:00", tz='GMT'), "days")
 			} else{
 				DateTime <- seq(ymd_hm("2015-01-01 00:00", tz='GMT'), ymd_hm("2015-12-31 23:00", tz='GMT'), "months")
 			}
 			
			data_xts <- xts(data, order.by = DateTime)
		})
		
		data.fig<-reactive({
			data_fig <- colSelect(data.xts(), input$colSelection) #input$colSelection
			dates_text <- paste(as.character(input$dateRange), collapse ="/")
			data_fig <- data_fig[dates_text]
 		})
 		
 		# output$filepath_val <- renderText({
 			# input$file1[input$file1[,'name']==input$path,'datapath']
 		# })
 		
 		output$map <- renderDygraph({
			dygraph(data.fig()) %>% 
				dySeries(strokeWidth=2) %>%
				dyRangeSelector(retainDateWindow=TRUE) %>%
				dyRoller(rollPeriod=1) %>%
				dyOptions(useDataTimezone=TRUE)
   })
  }
)