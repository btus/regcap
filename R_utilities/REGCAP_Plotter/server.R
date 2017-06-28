# server.R

library(xts); library(dygraphs); library(data.table); library(lubridate)

#increase max file upload to 160 mb
options(shiny.maxRequestSize=600*1024^2)

#colSelct takes a data object and a flexible number of column indexes, and returns a new data object including just those columns selected by the user. The dygraph plotting utility requires an entire data object, rather than selecting columns from existing object.

colSelect=function(dat_obj, col1, ...){
	return(dat_obj[,c(col1, ...)])
}

add_shades <- function(x, start, end, ...) {
  for(i in 1:length(start)) {
    x <- dyShading(x, from = start[i], to = end[i], ...)
  }
  x
}

shinyServer(

	function(input, output) {
		
 		output$list_of_files <- renderUI({
 			file_list<-input$file1[,'name']
 			checkboxGroupInput("path", "Choose File", file_list)
		})
		
		output$list_of_variables <- renderUI({
			#filepath<-input$file1[input$file1[,'name']==input$path,'datapath']
			#variables_list<-colnames(fread(filepath, header=TRUE, sep='\t'))
			variables_list<-colnames(data.xts())
			checkboxGroupInput("colSelection", "Columns to Plot", variables_list)
		})
		
		output$list_of_indices <- renderUI({
			#filepath<-input$file1[input$file1[,'name']==input$path,'datapath']
			#variables_list<-colnames(fread(filepath))
			#variables_list<-colnames(data.xts())
			Ind_cols<-c('fan1', 'fan2', 'fan3', 'fan4', 'fan5', 'fan6', 'rivecOn', 'occuipied', 'RHind60', 'RHind70', 'None', 'Heat', 'Cool', 'Vent')
			#Ind_cols<-grep('Index', variables_list)
			#Ind_cols2<-grep('illage', variables_list)
			#checkboxGroupInput("colSelection", "Columns to Plot", variables_list[-Ind_cols])
			checkboxGroupInput("IndcolSelection", "Index Values to Plot (Must Select Two Buttons)", Ind_cols)
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
 		
 		starts<-reactive({
  			start_vals<-which(diff(as.numeric(data.xts()[,input$IndcolSelection[[1]]]))==1)
			index(data.xts())[start_vals]
		})
		
		ends<-reactive({
  			end_vals<-which(diff(as.numeric(data.xts()[,input$IndcolSelection[[1]]]))==-1)
			index(data.xts())[end_vals]
		})
		
		starts_2<-reactive({
			#if(length(input$IndcolSelection)>1){
  			start_vals2<-which(diff(as.numeric(data.xts()[,input$IndcolSelection[[2]]]))==1)
			index(data.xts())[start_vals2]
			#}
		})
		
		ends_2<-reactive({
			#if(length(input$IndcolSelection)>1){
  			end_vals2<-which(diff(as.numeric(data.xts()[,input$IndcolSelection[[2]]]))==-1)
			index(data.xts())[end_vals2]
			#}
		})

 		output$map <- renderDygraph({
			dygraph(data.fig()) %>% 
				dySeries(strokeWidth=2) %>%
				dyRangeSelector(retainDateWindow=TRUE) %>%
				dyRoller(rollPeriod=1) %>%
				dyOptions(useDataTimezone=TRUE) %>%
				add_shades(starts(), ends()) %>%
				add_shades(starts_2(), ends_2(), color = "#FFE6E6")
   })
  }
)
