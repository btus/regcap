# ui.R
library(xts); library(dygraphs); library(data.table); library(lubridate)
# path<-"/Users/brennanless/GoogleDrive/SmartVentilation_HumidityControl/ExampleRCO"
# setwd(path)
# hour_files<-list.files()[grep(".hours", list.files())]
# example<-fread(hour_files[1], header=TRUE, sep=",")
# example<-example[,Index := NULL]
# col_list<-colnames(example)

shinyUI(fluidPage(
  titlePanel("REGCAP Plotter"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Plot time series data from hourly .rco data files."),
      
      fileInput('file1', 'Choose file to upload',
      			multiple = TRUE,
                accept = c(
                  'text/csv',
                  'text/comma-separated-values',
                  'text/tab-separated-values',
                  'text/plain',
                  '.csv',
                  '.tsv',
                  '.rco',
                  '.hours',
                  '.days',
                  '.hum',
                  '.months'
                )
      ),
      
      uiOutput('list_of_files'),
      
      dateRangeInput('dateRange',
      	label = 'Date range input: yyy-mm-dd',
      	start = '2015-01-01',
      	end = '2015-12-31'
      	),
      
      uiOutput('list_of_variables'),
      
      uiOutput('Second_yaxis'),
    
      uiOutput('list_of_indices')
      


#     
      # selectInput("var", 
        # label = "Choose file to display",
        # choices = input$file1[,'name'],
        # selected = input$file1[1,'name']),
      
      # checkboxGroupInput("colSelection",
        # label = "Columns to Plot",
        # choices = col_list,
        # selected = col_list[1]) #,		        
    ),
    mainPanel(dygraphOutput("map"))
  )
))