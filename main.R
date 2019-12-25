#Load necessary packages
if (!'shiny' %in% installed.packages()) install.packages("shiny") 
  library(shiny)
if (!'shinydashboard' %in% installed.packages()) install.packages("shinydashboard") 
  library(shinydashboard)
if (!'ggplot2' %in% installed.packages()) install.packages("ggplot2") 
  library(ggplot2)

if (!'remotes' %in% installed.packages()) install.packages("remotes") 
if (!'rPeaks' %in% installed.packages()) remotes::install_github("jrminter/rPeaks")
library(rPeaks)

if (!'Rcpp' %in% installed.packages()) install.packages("Rcpp") 
  library(Rcpp)
if (!'coin' %in% installed.packages()) install.packages("coin") 
  library(coin)
if (!'DT' %in% installed.packages()) install.packages("DT") 
  library(DT)
if (!'multcomp' %in% installed.packages()) install.packages("multcomp") 
  library(multcomp)
if (!'dplyr' %in% installed.packages()) install.packages("dplyr") 
  library(dplyr)
if (!'tidyr' %in% installed.packages()) install.packages("tidyr") 
  library(tidyr)


##source functions
source("functions.R")
source("ui.R")
source("server.R")

##run app
shinyApp(ui = ui, server = server)
