# server.R

library(xts); library(dygraphs); library(data.table); library(lubridate)

#increase max file upload to 160 mb
options(shiny.maxRequestSize=160*1024^2) 

#colSelct takes a data object and a flexible number of column indexes, 
#and returns a new data object including just those columns selected by the user. 
#The dygraph plotting utility requires an entire data object, 
#rather than selecting columns from existing object.

mycol <- rgb(184, 184, 184, max = 255, alpha = 85, names = "grey50")
mycol2 <- rgb(0, 184, 184, max = 255, alpha = 85, names = "50")
mycol3 <- rgb(0, 0, 184, max = 255, alpha = 85, names = "50b")
shade_cols<-c(mycol, mycol2, mycol3)

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
		
		output$Second_yaxis <- renderUI({
			checkboxGroupInput("SecYaxis", "Variables on 2nd y-axis", input$colSelection)
		})
		
		output$list_of_indices <- renderUI({
			#filepath<-input$file1[input$file1[,'name']==input$path,'datapath']
			#variables_list<-colnames(fread(filepath))
			#variables_list<-colnames(data.xts())
			Ind_cols<-c('fan1', 'fan2', 'fan3', 'fan4', 'fan5', 'fan6', 'rivecOn', 'occupied', 'RHind60', 'RHind70', 'None', 'Heat', 'Cool', 'Vent')
			#Ind_cols<-grep('Index', variables_list)
			#Ind_cols2<-grep('illage', variables_list)
			#checkboxGroupInput("colSelection", "Columns to Plot", variables_list[-Ind_cols])
			checkboxGroupInput("IndcolSelection", "Index Values to Plot", Ind_cols)
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
 		
 		#First shading choice.
 		starts<-reactive({
 			dates_text <- paste(as.character(input$dateRange), collapse ="/")
 			data_adj <- as.numeric(data.xts()[,input$IndcolSelection[[1]]][dates_text])
 			data_adj[1] <- 0
 			data_adj[length(data_adj)]<-0
 			diffs<-diff(data_adj)
  			start_vals<-which(diffs==1)
			index(data.xts()[dates_text])[start_vals]
		})
		
		ends<-reactive({
			dates_text <- paste(as.character(input$dateRange), collapse ="/")
			data_adj2 <- as.numeric(data.xts()[,input$IndcolSelection[[1]]][dates_text])
 			data_adj2[1] <- 0
 			data_adj2[length(data_adj2)]<-0
			diffs<-diff(data_adj2)
  			end_vals<-which(diffs==-1)
			index(data.xts()[dates_text])[end_vals]
		})
		
		#Second shading choice.
		 starts_2<-reactive({
 			dates_text <- paste(as.character(input$dateRange), collapse ="/")
 			data_adj3 <- as.numeric(data.xts()[,input$IndcolSelection[[2]]][dates_text])
 			data_adj3[1] <- 0
 			data_adj3[length(data_adj3)]<-0
 			diffs<-diff(data_adj3)
  			start_vals<-which(diffs==1)
			index(data.xts()[dates_text])[start_vals]
		})
		
		ends_2<-reactive({
			dates_text <- paste(as.character(input$dateRange), collapse ="/")
			data_adj4 <- as.numeric(data.xts()[,input$IndcolSelection[[2]]][dates_text])
 			data_adj4[1] <- 0
 			data_adj4[length(data_adj4)]<-0
			diffs<-diff(data_adj4)
  			end_vals<-which(diffs==-1)
			index(data.xts()[dates_text])[end_vals]
		})
		
		output$map <- renderDygraph({
		
			if(length(input$IndcolSelection) == 1){
				dygraph(data.fig()) %>% 
					dySeries(input$SecYaxis, axis = 'y2') %>%
					dySeries(strokeWidth=2) %>%
					dyRangeSelector(retainDateWindow=TRUE) %>%
					dyRoller(rollPeriod=1) %>%
					dyOptions(useDataTimezone=TRUE) %>%
					add_shades(starts(), ends(), color = "#FFE6E6") 
					
			} else if(length(input$IndcolSelection) == 2){
				dygraph(data.fig()) %>% 
					dySeries(input$SecYaxis, axis = 'y2') %>%
					dySeries(strokeWidth=2) %>%
					dyRangeSelector(retainDateWindow=TRUE) %>%
					dyRoller(rollPeriod=1) %>%
					dyOptions(useDataTimezone=TRUE) %>%
					add_shades(starts(), ends(), color = "#FFE6E6") %>%
					add_shades(starts_2(), ends_2(), color = "#CCEBD6")
			} else {
				dygraph(data.fig()) %>% 
					dySeries(input$SecYaxis, axis = 'y2') %>%
					dySeries(strokeWidth=2) %>%
					dyRangeSelector(retainDateWindow=TRUE) %>%
					dyRoller(rollPeriod=1) %>%
					dyOptions(useDataTimezone=TRUE)
			}
		})	
	}
)