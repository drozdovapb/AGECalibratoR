################
###### UI ######
################

library(shiny)
library(shinydashboard)
library(DT)
source("./functions.R")
library(rPeaks)
library(ggplot2)
library(rPeaks)
library(Rcpp)
library(coin)
library(DT)
library(shinydashboard)
library(multcomp)
library(dplyr)
library(tidyr)

##############
### HEADER ###
##############

header <- dashboardHeader(title = 'AGECalibratoR')


###############
### SIDEBAR ###
###############


sidebar <- dashboardSidebar(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "www/custom.css")# css style
  ),# end of tags
  sidebarMenu(
    menuItem(startExpanded = T,
      "Input File",
      fileInput('file1', 'Choose your CSV File',
                accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
      checkboxInput('header', 'Header', TRUE),
      radioButtons('sep', 'Separator',
                   c(Comma=',', Semicolon=';', Tab='\t'))
    ),# end of menuItem
    menuItem(startExpanded = T, 
      "Parameters",
      numericInput(inputId = "t", 
                   label = "Time (min)",
                   value = 240), #end of field
      ##field for user adjustment: Voltage
      numericInput(inputId = "V", 
                   label = "Voltage (V)",
                   value = 30), #end of field
      ##field for user adjustment: The distanse between electrodes (m)
      numericInput(inputId = "d", 
                   label = "The distanse between electrodes (cm)",
                   value = 15), #end of field
      numericInput(inputId = "mppx", 
                   label = "px per cm",
                   value = 10), #end of field
      ##field for user adjustment: The percentage of agarose in gel (%)
      
      numericInput(inputId = "Tgel", 
                   label = "The percentage of agarose in gel (%)",
                   value = 1.5)#, #end of field
    )# end of menuItem
  )# end of sidebarMenu
)# end of dashboardSidebar


############
### BODY ###
############


body <- dashboardBody(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "www/custom.css")# css style
  ),# end tags
  tabBox(title = "AGECalibrator", width = 12,
         tabPanel("1. Data input",
                  fluidRow(
                    box(title = textOutput('message'), width = 12, status = "primary", solidHeader = T,
                        dataTableOutput("trace_table"), style = "height:200px; overflow-y: scroll;overflow-x: scroll;"
                    ),# end of box
                    box(title = 'Preprocessing', width = 12, status = "primary", solidHeader = T, 
                        textOutput('message2'),
                        uiOutput("selectLadder"),
                        uiOutput("selectBackground"),
                        uiOutput("Run")
                    )# end of box
                    
                  )# end of fluidRow
         ),# end of tabPanel
         tabPanel('2. Ladder peaks',
                  fluidRow(
                    column(width = 12,#8,
                           box(width = NULL, status = 'primary', solidHeader = T,
                               textOutput('message3'),
                               ##field for user input: DNA ladder values
                               column(width = 5, textInput(inputId = 'ladder', 
                                         label = "DNA ladder (heavy to light)", 
                                         value = "10,8,6,5,4,3,2.5,2,1.5,1, .75, .5")),
                               column(width = 2, numericInput(inputId = "sigma", 
                                            label = "Sigma",
                                            value = 5)),
                               ##field for user adjustment: threshold(what to consider a peak)
                               column(width = 2, numericInput(inputId = "threshold", 
                                            label = "Threshold",
                                            value = 10))
                           ),# end of box
                           box(title = 'Plot', width = NULL, status = 'primary', solidHeader = T,
                               plotOutput('scatterplot'), style = "height:800px;"
                           )# end of box
                    )# end of column
                  )# end of fluidRow
         ),# end of tabPanel
         tabPanel('3. Aggregate weights',
                  fluidRow(
                    column(width = 4,
                           box(status = 'primary', solidHeader = T, width = NULL,
                               HTML("Before doing this please make sure the DNA ladder table is correct...<br> Ready?"),
                               hr(),
                               actionButton(inputId = "regression", label = "Calculate reg. model")
                           ),# end of box
                           box(status = 'success', solidHeader = T, width = NULL,
                               downloadButton("downloadData", HTML("Download <br/> percentile data"))
                           ), 
                           box(status = 'success', solidHeader = T, width = NULL,
                               downloadButton("downloadPeaks", HTML("Download <br/> peak data"))
                           )# end of box
                    ),# end of column
                    column(width = 8,
                           box(title = 'Percentiles', status = 'primary', solidHeader = T, width = NULL,
                               dataTableOutput(outputId = "percentiles"), style = "height:125px; overflow-y: scroll;overflow-x: scroll;"
                           ),
                           box(title = 'Peaks', status = 'primary', solidHeader = T, width = NULL,
                               dataTableOutput(outputId = "PeakMW"), style = "height:125px; overflow-y: scroll;overflow-x: scroll;"
                           )# end of box
                    )# end of column
                    
                  ),# end of fluidRow
                  fluidRow(
                    box(title = 'MW distribution of aggregates (asterisk for the main peak)', 
                        status = 'primary', solidHeader = T, width = 12,
                        plotOutput(outputId = "boxplot")
                    )# end of box
                  )# end of fluidRow
         ),# end of tabPanel
         tabPanel('(4.) Sample comparison',
                  ####################################
                  ####### interactive grouping #######
                  ####################################
                  fluidRow(
                    column(width = 12,
                           box(width = NULL, status = 'primary', solidHeader = T,
                               uiOutput("groups")
                           )# end of box
                    ), # end of column
                    column(width = 12,
                          box(title = 'test', width = NULL, status = 'primary', solidHeader = T,
                              dataTableOutput(outputId = 'summary'),
                              style = "overflow-y: scroll;overflow-x: scroll;"
                          )# end of box
                    )# end of column
                  ),# end of fluidRow
                  fluidRow(
                    column(width = 6,
                           box(width = NULL, status = 'primary', solidHeader = T,
                               uiOutput('variables'),
                               uiOutput("selectStatMethod"),
                               uiOutput('selectParam'),
                               actionButton(inputId = "compare_groups", label = "Compare groups")
                           )# end of box
                    ),# end of column
                    column(width = 6,
                          box(width = NULL, status = 'primary', solidHeader = T,
                              textOutput('testngroups'),
                              verbatimTextOutput("mt"),
                              style = "overflow-y: scroll;overflow-x: scroll;"
                          )# end of box
                    )# end of column
                  ),# end of fluidRow
                  fluidRow(
                    column(width = 12,
                           box(width = NULL, status = 'primary', solidHeader = T,
                               plotOutput('boxplot_medians')
                           )
                    )# end of column
                  )# end of fluidRow
                  #####################################
                  #### end of interactive grouping ####
                  #####################################
                  
                  
                  ### some alternative code for manual grouping [not ready yet]
                  # conditionalPanel(condition = 'output["percentiles"]',
                  #                  ##Statistical analysis
                  #                  hr(),
                  #                  
                  #                  column(width = 7, 
                  #                         HTML("Do you need to compare groups of samples? <br>
                  #           If yes, please enter some description of your groups as a comma-separated list. <br>
                  #           Example: s,s,s,s,s,s,w,w,w,w,w,w <br>
                  #           If you want to leave out some sample(s), lavel them as NA. <br>
                  #           Example: s,s,s,s,NA,NA,w,w,w,w,w,w"),
                  #                         ##input sample groups
                  #                         textInput(inputId = 'sample_groups', 
                  #                                   label = "",
                  #                                   value = "s,w"), #initial value
                  #                         uiOutput("selectStatMethod")
                  #                  ), #end of column
                  #                  
                  #                  column(width = 5,    
                  #                         ##print the table to allow the user to monitor any changes they make
                  #                         tableOutput("groupsTbl")
                  #                  ), #end of column
                  #                  
                  #                  
                  #                  HTML("Sample groups ready?"),
                  #                  actionButton(inputId = "compare_groups", label = "Compare groups")
                  
                  
         )# end of tabPanel
         )# end of tabBox
                  )# dashboardBody end


### UI CONSRTUCTOR ###

ui <- dashboardPage(header = header,
                    sidebar = sidebar,
                    body = body
)# dashboardPage end

