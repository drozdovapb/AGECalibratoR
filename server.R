################
#### Server ####
################

## the main server function

server <- function(input, output) { 
  
  output$message <- renderText("Your table looks like this: ")
  output$trace_table <- renderDataTable({
    inFile <- input$file1
    if (is.null(inFile)) 
      return(NULL)
    datatable(read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote), options = list(dom = 't', ordering=F), rownames = FALSE)
  })
  output$message2 <- renderText({
    inFile <- input$file1
    if (is.null(inFile)) 
      return(NULL)  })
  output$Run <- renderUI({
    ## take input; return 'NULL', if input file is not specified
    inFile <- input$file1
    if (is.null(inFile)) 
      return(NULL)
    ## read data.frame from the input file and create the drop-down list
    tbl <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)
    selectInput("Run","Select column with y coordinates in px", choices = names(tbl), selected = "X")
  })
  
  ## create drop-down list to choose name of variable, corresponding to 
  ## optical density profile (ODP) of DNA ladder
  output$selectLadder <- renderUI({
    ## take input; return 'NULL', if input file is not specified
    inFile <- input$file1
    if (is.null(inFile)) 
      return(NULL)
    ## read data.frame from the input file and create the drop-down list
    tbl <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)
    selectInput("ladder_name","Select the DNA ladder column", choices = names(tbl), selected = "ladder")
  })
  
  ## create drop-down list to choose name of variable, corresponding to 
  ## optical density profile (ODP) of background signal for western blotting
  output$selectBackground <- renderUI({
    ## take input; return 'NULL', if input file is not specified
    inFile <- input$file1
    if (is.null(inFile)) 
      return(NULL)
    ## read data.frame from the input file and create the drop-down list
    tbl <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)
    selectInput("background_name","Select column with WB background signal", choices = names(tbl), selected = "background")
  })
  output$message3 <- renderText({
    ## take input; return 'NULL', if input file is not specified
    inFile <- input$file1
    if (is.null(inFile)) 
      return('ladder')
    "Please enter a vector (comma delimited) of the MW ladder bands in kb (heavy to light). 
    Adjust the Sigma and threshold values to correctly identify the bands; 
    if you need to exclude particular bands from the downstream analysis, 
    please mark them with NA. 
    FYI: 1kb SibEnzyme ladder is 10, 8, 6, 5, 4, 3, 2.5, 2, 1.5, 1, 0.75, 0.5. \n "
  })
  model <- reactive ({
    
    ## take input; return 'NULL', if input file is not specified
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    
    ## read data.frame from the input file and create the drop-down list
    Data <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)
    
    ## Units transfer
    t <- input$t * 60 # minutes to seconds
    d <- input$d / 100 # cm to meters for distance between electrodes
    mppx <- 1 / (100 * input$mppx) # N of pixels in cm to meters in 1 px

    
    ## calculation of coordinates of DNA ladder bands with the function 'DNALadderRun' (returns numeric vector)
    Run <- DNALadderRun(background = 0, ladder=Data[,input$ladder_name], 
                        sigma = input$sigma, threshold = input$threshold)
    
    ## create the plot of DNA ladder ODP with marked bands positions
    req(input$ladder)
    userladder <- as.numeric(unlist(strsplit(input$ladder,",")))
    p <- DNALadderPlot(background = 0, ladder=Data[,input$ladder_name], peaks = list(pos = Run, names = userladder))    
    
    ## create a vector for values with the MW ladder bands (user input)
    Size <- ConvertKbp2MW(userladder) ##if in kb use Size <- userladder
    ## create the table with the vector of DNA ladder coordinates and NA vectors for the ladder MWs and mobility (mu)
    Ld <- data.frame("px" = Run, "kb" = NA, 'mu' = NA)
    
    ## cut the table and add MW ladder in Da
    n <- length(userladder)
    Ld$"kb"[1:n] <- userladder[1:n]

        ## these are the converted units
    Ld$"mu"[1:n] <- PxToMobility(Run, t = t, V = input$V, d = d, mppx = mppx) 
    
    return(list(Plot = p, ld = Ld, data = Data, t = t, d = d, mppx = mppx))

  })
  ## create output plot of ODP with peaks marked. 
  output$scatterplot <- renderPlot({
    model()$Plot}, width=350, height=750)
  ## create output table of DNA ladder runs and MWs
  output$ld <- renderDataTable({
    datatable(model()$ld[,c('px','kb')], options = list(dom = 't', ordering=F))})
  
  ## create the boxplot for pseudosamples of molecular weights for proreins samples
  TestMSG <- eventReactive(input$regression, {
    ## remove rows with NA values 
    ladder_data <- model()$ld[complete.cases(model()$ld),]
    
    ## convert size from kDa to Da and kbp to bp
    ladder_data$kb <- 1000 * ladder_data$kb
    
    ## calculate model to transform mobilities to lengths (bp)
    RTL_model <- Run2LengthModel(ladder_data, Tgel = input$Tgel)
    
    input_table <- as.data.frame(model()$data)
    
    ## calculate pseudosamples of molecular weights for proreins samples
    MWeights <- MWCalculator(data = input_table, background = input$background_name, 
                             model = RTL_model, excluded = list(input$background_name, input$ladder_name, input$Run),
                             t = model()$t, V = input$V, d = model()$d, mppx = model()$mppx)
    
    ## compute lower and upper whiskers for pseudosamples of molecular weights for proreins samples
    stats <- boxplot.stats(MWeights$MW)
    upper_outlier <- min(stats$out[stats$out > stats$stats[5]])
    
    ## create the table with percentile values (0 to 100) describing destribution of MWs of different samples
    percentiles <- percentilesCalculator(MWeights)
    
    ## create the table with Peak coordinates 
    goodnames = names(input_table)[!(names(input_table) %in% list(input$background_name, input$ladder_name, input$Run))]
    Peaks = lapply(X = input_table[goodnames], FUN = function(x) MainPeak(x, sigma = input$sigma))
    MWPeaks = as.data.frame(lapply(X = Peaks, FUN = function(x){
      RunToMW(x, model = RTL_model, t = model()$t, V = input$V, d = model()$d, mppx = model()$mppx)}))
    
    peaksdf <- data.frame(sample = names(MWPeaks), peaks = t(MWPeaks)[, 1])
    
    
    ## create the boxplot for pseudosamples of molecular weights for proreins samples
    b <- ggplot(data = MWeights, mapping = aes(y = MW, x = sample)) + 
      theme_classic(base_size = 13) +
      geom_boxplot(outlier.shape=NA, notch = FALSE) + 
      ### and add the peak to have some idea on how well they correspond to each other
      geom_point(data=peaksdf, aes(x=sample, y=peaks), col="violetred", shape=8) +
      coord_cartesian(ylim = c(0, upper_outlier)) + 
      labs(y = "Molecular weight (kDa)") + theme(axis.text.x=element_text(angle=90, hjust=1))
    
    return(list(Boxplot = b, MWeights = MWeights, percentiles = percentiles, PeakMW=MWPeaks))
    
  })
  
  ## render boxplot for pseudosamples of molecular weights for proreins samples
  output$boxplot <- renderPlot({TestMSG()$Boxplot})
  ## render the table with percentile values (0 to 100) describing destribution of MWs of different samples
  output$percentiles <- renderDataTable({datatable(TestMSG()$percentiles[c(1,26,51,76,101),], 
                                                   options = list(dom = 't', ordering=F))})
  ## Downloadable csv of percentiles table 
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$file1, "_percentiles.csv", sep = "")
    },
    content = function(file) {
      write.csv(TestMSG()$percentiles, file, row.names = FALSE, quote = FALSE)
    })
  ##and a similar downloadable table with peak values
  output$downloadPeaks <- downloadHandler(
    filename = function() {
      paste(input$file1, "_peaks.csv", sep = "")
    },
    content = function(file) {
      write.csv(TestMSG()$PeakMW, file, row.names = FALSE, quote = FALSE)
    })
  
  
  ## render the table with MW of peak value 
  output$PeakMW <- renderDataTable({datatable(TestMSG()$PeakMW, options = list(dom = 't', ordering=F), rownames = FALSE)})
  
  
  ## Assign samples to groups
  enterGroups <- reactive ({
    ## Process user input (a string contaning all groups)
    tblgroups <- unlist(strsplit(input$sample_groups, ","))
    ## Construct a table with sample names & groups they belong to
    groupsTbl <- data.frame(Names = names(TestMSG()$percentiles)[-1], group = NA)
    ## Populate the table with user's values
    groupsTbl$group[1:length(tblgroups)] <- tblgroups
    #return the result
    return(list(groupsTbl = groupsTbl))
  })
  
  ## render the table with sample names & groups for visual control
  output$groupsTbl <- renderDataTable(datatable(enterGroups()$groupsTbl, options = list(dom = 't', ordering = F), rownames = FALSE))
  
  
  ####################################
  ####### interactive grouping #######
  ####################################
  ##Compare all pairs or everything to the control level

  data_for_IQR <- reactive({
    df <- as_data_frame(TestMSG()$percentiles)
    df <- df %>% select(-percentiles) %>%
      gather(key = "group", value)
    value <- df$value
    df <- df %>% select(-value) %>% separate(group, as.character(0:lengths(regmatches(df$group[1], gregexpr("_", df$group[1])))+1)) %>% 
      rename_all(function(x) paste0('name_part_', x)) %>% 
      mutate_all(funs(factor))
    df$value <- value
    return(df)
  })
  output$groups <- renderUI({
    df <- data_for_IQR()
    ## select which part of the name use to define each sample (default is all)
    selectInput(inputId = "grouper", label = "Group variable", 
                choices = names(df[, names(df) != "value"]), 
                multiple = TRUE, selected = names(df[, names(df) != "value"]))
  })
  summary_data <- reactive({
    req(input$grouper)
    datka <- data_for_IQR() %>%
      group_by(!!!rlang::syms(input$grouper)) %>% 
      summarise(Median = median(value), 
                IQR = IQR(value))#,
                ##it's not gonna work like this
                #mainPeak = MainPeak(value, sigma = input$sigma))
    datka$MainPeak <- as.vector(t(TestMSG()$PeakMW[1, ]))
      print(datka)

    return(datka)
  })
  output$summary <- DT::renderDataTable({
    data <- summary_data()
    datatable(data, options = list(dom = 't', ordering=F), rownames = F)
  })
  output$selectStatMethod <- renderUI({
    selectInput(inputId = "statMethod", label = "How should we perform comparison of groups? \n
                Please note that if you choose Dunnett's test you need to make your reference 
                (control) level alphabetically first.", 
                choices = c("Compare all pairs (Mann-Whitney test)", "Compare all groups to one control (Dunnett's test)"), 
                selected = "Compare all pairs (Mann-Whitney test)")
  })
  output$selectParam <- renderUI({
    selectInput(inputId = "statParam", label = "Select param:", 
                choices = c("Median", "IQR", "MainPeak"), 
                selected = "Median")
  })
  outVar <- reactive({
    vars <- input$grouper
    return(vars)
  })
  output$variables = renderUI({
    selectInput('variables2', 'Variables', outVar())
  })
  StatTest <- eventReactive(input$compare_groups, {
    ## Get data to work with 
    data <- summary_data()
    
    var_n <- input$variables2
    param <- input$statParam
    
    mygroups <- list()
    mygroups$group <- data[[var_n]]
    mygroups$parameter <- data[[param]]
    mygroups <- as_data_frame(mygroups)

    ## The number of groups to user check (will be shown for user's check)
    ngroups <- length(levels(mygroups$group))
    ## Perform test
    if (input$statMethod == "Compare all groups to one control (Dunnett's test)") {
      fit <- aov(parameter ~ group, mygroups)
      test_results <- summary(glht(fit, linfct=mcp(group="Dunnett")))
    } else{
      ##If not Dunnett
      ## Two groups => Manny-Whitney
      if (ngroups == 2 ) {
        test_results <- wilcox_test(mygroups$parameter ~ mygroups$group)
      }
      ## More 2 groups => with correction for multiple testing
      if (ngroups > 2 ) {
        test_results <- pairwise.wilcox.test(mygroups$parameter, mygroups$group, p.adjust.method = "fdr")
      }
    } # end of else for Wilcoxon
    ## Plot group medians
    p <- ggplot(mygroups, aes(x = group, y = parameter)) + 
      geom_boxplot() + geom_point() + 
      expand_limits(y=0) +
      ylab("Aggregate weight (MDa)") + 
      # ggtitle(HTML("Distribution of median electrophoretic mobility values")) +
      theme_classic(base_size = 13)
    
    ## Return the results of testing and the boxplot
    return(list(mg = mygroups, mt = test_results, ngroups = ngroups, bplot = p))
  })# end of eventReactive
  
  output$testngroups <- renderText(
    paste0(StatTest()$ngroups, " groups detected. \n 
              If this is not what you intended please check your input, particularly for extra spaces). \n \n
              If you use Dunnett's test and are not satisfied with the reference factor level, please make 
                     sure you encoded it with a word that would be first it your conditions are sorted 
                     by alphabetical order."))
  ## Render the result of the test
  output$mt <- renderPrint(StatTest()$mt)
  output$boxplot_medians <- renderPlot({StatTest()$bplot})

  
  #####################################
  #### end of interactive grouping ####
  #####################################
  
  }# end of server

