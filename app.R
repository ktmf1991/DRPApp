#shiny app model
require(shiny)
require(drc)
require(tidyverse)
require(DT)
require(plotly)
require(ggplot2)

source("function.R")

#---- Define UI for dataset viewer app ----
ui <- fluidPage(
  
  #---- App title ----
  titlePanel("Drug Response Profiling Quick Look"),
  
  #---- Sidebar layout with a input and output definitions ----
  sidebarLayout(
    
    #---- Sidebar panel for inputs ----
    sidebarPanel(
      
      # Patient name
      textInput("ptid", label = "Analysis ID", placeholder = "YYYYMM_ID"),
      
      
      #---- Input result file, 
      fileInput("file1", "Choose CSV File for patients", accept = ".csv"),
      # Input: Analysis button
      actionButton("do", label = "Analyse!"),
      # Input: Selector for choosing dataset ----
      # Patient
      checkboxGroupInput("pat", "Patients", choices = NULL, selected = NULL),
      selectizeInput(
        'drug', label = "Plot", choices = "", selected = "", multiple = T, options = list(maxItems = 1)
      ),
      #
      selectizeInput(
        'drug2', label = "EC50/IC50", choices = "", selected = "", multiple = T
      ),
      
      # Download button for EC and IC
      h3("Download Results"),
      downloadButton("dl1", label = "Download EC Results!"),
      downloadButton("dl2", label = "Download IC Results!"),
      downloadButton("dl4", label = "Download Raw Data!"),
      # Download button for plots
      downloadButton("dl3", label = HTML("<span style=\"font-size:14px;\">Download Plots!</span> <br/> 
                                         <span style=\"font-size:10px;\">(shaded area = 95% CI)</span>"),
                     style='height:60px'),
      checkboxInput("rib","Include 95% interval?",value = T),
      width = 2
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      h2("Drug response curve"),
      h3(textOutput("drug_title")),
      plotlyOutput("drm"),
      h2("EC/IC table"),
      tabsetPanel(
        id = 'dataset',
        tabPanel("EC50", DT::dataTableOutput("EC")),
        tabPanel("IC50", DT::dataTableOutput("IC")),
        tabPanel("Raw data point",DT::dataTableOutput("raw")),
      ),
      width = 6
      
    )
  )
)

# Define server logic to summarize and view selected dataset ----
server <- function(input, output, session) {
  # Read input file
  up_res <- reactive({
    inFile <- input$file1
    if (is.null(inFile)) {
      return(NULL)      
    }
    read.csv(inFile$datapath)
  })
  

  
  # Analysis using drm of the whole dataset
  observeEvent(input$do, {
    output$drug_title <- renderText({
      "Analysis in progress"
    })
    
    dat <<- up_res()
    
    dat.list <- split(dat,dat$Patient, drop = T)
    
    dat.list <- lapply(dat.list, df.trans)
    
    drm.list <- lapply(dat.list,drm.trans)
    
    require(drc)
    ll4.res <<- lapply(drm.list,function (y) lapply(y,function(x) drm(Percentage ~ Concentration_2, data = x, fct = LL.4(), na.action = na.omit)))
    
    plot.res <- lapply(seq_along(dat.list),function(x) plot.list(ll4.res[[x]],dat.list[[x]])) #each plot is in [[patient]]$Drug[[1 or 2]] or [[patient]][drug][[1 (always)]][[1 or 2]] 
    
    plot.res <- lapply(seq_along(plot.res), 
                       function(y) lapply(names(plot.res[[y]]),
                                          function(x) cbind(plot.res[[y]][[x]],
                                                            Patient = rep(y,nrow(plot.res[[y]][[x]])),
                                                            Treatment_Drug = x)
                       )
    )
    plot.res <- tdlc(plot.res)
    plot.res$Patient <- as.character(plot.res$Patient)
    plot.res <<- plot.res
    
    dat.list <- lapply(seq_along(dat.list),function(x) cbind(dat.list[[x]]))
    dat.list <- odlc(dat.list)
    plot.res$Patient <- as.factor(plot.res$Patient)
    dat.list <<- as.data.frame(dat.list)
    
    df <- split(dat.list,dat.list$Treatment_Drug)
    pr <- split(plot.res,plot.res$Treatment_Drug)
    
    exp.plot <- lapply(names(pr),function(x) plot.drug(df[[x]],pr[[x]]))
    names(exp.plot) <- names(pr)
    exp.plot <<- exp.plot
    
    exp.plot.rib <- lapply(seq_along(exp.plot), function(x) {
      g1 <- exp.plot[[x]]
      fits <- pr[[x]]
      g2 <- g1 + geom_ribbon(data = fits, aes(x=conc, y=p ,ymin=pmin, ymax=pmax, fill = factor(Patient)), alpha=0.2) +
        labs(colour = "Patient", fill = "Patient")
      return(g2)
    })
    names(exp.plot.rib) <- names(pr)
    exp.plot.rib <<- exp.plot.rib
    
    # Create EC and IC at different percentage
    ed <- lapply(ll4.res,function (y) lapply(y,function(x) ED(x,c(10,25,50,75,90))))
    ed <- lapply(ed,function(y) lapply(y,function(x) as.vector(t(x))))
    ed <- lapply(ed,function(y) data.frame(matrix(unlist(y), nrow=length(y), byrow=TRUE)))
    ed2 <- lapply(seq_along(ed), function(y) cbind(names(ll4.res[[y]]),rep(y,nrow(ed[[y]])),ed[[y]]))
    ed2 <- lapply(ed2, function(x) {
      dat <- x
      colnames(dat) <- c("Drug","Patient",paste0(rep(paste0("EC",c(10,25,50,75,90)),each = 2),c("","_SE")))
      return(dat)
    })
    ed2 <- odlc(ed2)
    ed2 <<- as.data.frame(ed2)

    ic <- lapply(ll4.res,function (y) lapply(y,function(x) ED(x,c(10,25,50,75,90),type = "absolute")))
    ic <- lapply(ic,function(y) lapply(y,function(x) as.vector(t(x))))
    ic <- lapply(ic,function(y) data.frame(matrix(unlist(y), nrow=length(y), byrow=TRUE)))
    ic2 <- lapply(seq_along(ic), function(y) cbind(names(ll4.res[[y]]),rep(y,nrow(ic[[y]])),ic[[y]]))
    ic2 <- lapply(ic2, function(x) {
      dat <- x
      colnames(dat) <- c("Drug","Patient",paste0(rep(paste0("IC",c(10,25,50,75,90)),each = 2),c("","_SE")))
      return(dat)
    })
    ic2 <- odlc(ic2)
    ic2 <<- as.data.frame(ic2)
    
    # Create EC50 and IC50 for display
    ed <- ed2 %>% select(Patient, Drug,EC50,EC50_SE)
    ed$EC50 <- round(ed$EC50,2)
    ed$EC50_SE <- round(ed$EC50_SE,2)
    ed <- split(ed,ed$Patient,drop = F)
    ed <- lapply(ed, function(x) x %>% select(Drug,EC50, EC50_SE))
    ed <- merge.list(ed, by = "Drug")
    colnames(ed) <- c("Drug",paste0(rep(c("EC50","EC50_SE"),length(unique(dat$Patient))),"(Patient ",
                                    rep(seq_along(unique(dat$Patient)),each = length(unique(dat$Patient))),")"))
    ed <<- ed
     
    ic <- ic2 %>% select(Patient, Drug,IC50,IC50_SE)
    ic$IC50 <- round(ic$IC50,2)
    ic$IC50_SE <- round(ic$IC50_SE,2)
    ic <- split(ic,ic$Patient,drop = F)
    ic <- lapply(ic, function(x) x %>% select(Drug,IC50, IC50_SE))
    ic <- merge.list(ic, by = "Drug")
    colnames(ic) <- c("Drug",paste0(rep(c("IC50","IC50_SE"),length(unique(dat$Patient))),"(Patient ",
                                    rep(seq_along(unique(dat$Patient)),each = length(unique(dat$Patient))),")"))
    ic <<- ic
    
    # Update drug selection list
    uni.drug <- unique(dat$Treatment_Drug)
    updateSelectInput(session, "drug",
                      choices = c("All",uni.drug), selected = "All"
    )
    updateSelectInput(session, "drug2",
                      choices = c("All",uni.drug), selected = "All"
    )
    
    # Update patient number
    uni.pat <- seq_along(unique(dat$Patient))
    updateCheckboxGroupInput(session, "pat",
                             choices = uni.pat, selected = NULL
    )
  })

  # Update table after selecting drug << change this to type in with suggested list for long drug list
  observeEvent(input$drug,{
    output$drm <- renderPlotly({
      req(input$drug)
      req(input$pat)
      if (input$drug == "All"){
        NULL
      } else{
        ds <- paste0(input$drug)
        ps <- paste0(input$pat,collapse = "|")

        ggfits <- plot.res %>% filter(grepl(ds,Treatment_Drug),grepl(ps,Patient))
        ggdf <- dat.list %>% filter(grepl(ds,Treatment_Drug),grepl(ps,Patient))
        
        g <<- plot.drug(ggdf,ggfits)
      }
    })
  })
  observeEvent(input$drug,{
    output$drug_title <- renderText({
      req(input$drug)
      if (input$drug == "All"){
        "Choose a drug"
      } else{
        input$drug
      }
    })
  })
  # Create EC/IC table
  observeEvent(input$do,{
    output$EC <- DT::renderDataTable({
      req(input$drug2)
      if (input$drug2 == "All"){
        DT::datatable(ed)
      } else{
        DT::datatable(ed %>% filter(Drug %in% input$drug2))
      }
    })
  })
  observeEvent(input$do,{
    output$IC <- DT::renderDataTable({
      req(input$drug2)
      if (input$drug2 == "All"){
        DT::datatable(ic)
      } else{
        DT::datatable(ic %>% filter(Drug %in% input$drug2))
      }
    })
  })
  observeEvent(input$do,{
    output$raw <- DT::renderDataTable({
      req(input$drug)
      if (input$drug == "All"){
        DT::datatable(dat.list[,c("Patient","Well","Treatment_Drug","Concentration_1","Unit_1","Concentration_2",
               "Unit_2","Remark","Count","Percentage")]) %>%
          formatRound(10, 2)
      } else{
        DT::datatable(dat.list[dat.list$Treatment_Drug == input$drug|dat.list$Treatment_Drug == "DMSO",
                               c("Patient","Well","Treatment_Drug","Concentration_1","Unit_1","Concentration_2",
                             "Unit_2","Remark","Count","Percentage")] ) %>%
          formatRound(10, 2)
      }
    })
  })
  # 
  # Click EC/IC table to change plot
  observeEvent(input$EC_cell_clicked, {
    info = input$EC_cell_clicked
    # do nothing if not clicked yet, or the clicked cell is not in the 1st column
    if (is.null(info$value)||info$col != 1) return()
    output$drm <- renderPlotly({
      ds <- paste0(info$value)
      ps <- paste0(input$pat,collapse = "|")
      
      ggfits <- plot.res %>% filter(grepl(ds,Treatment_Drug),grepl(ps,Patient))
      ggdf <- dat.list %>% filter(grepl(ds,Treatment_Drug),grepl(ps,Patient))
      
      g <<- plot.drug(ggdf,ggfits)
    })
    output$drug_title <- renderText({info$value})
    output$raw <- DT::renderDataTable({
      DT::datatable(dat.list[dat.list$Treatment_Drug == info$value|dat.list$Treatment_Drug == "DMSO",
                        c("Patient","Well","Treatment_Drug","Concentration_1","Unit_1","Concentration_2",
                          "Unit_2","Remark","Count","Percentage")] ) %>%
        formatRound(10, 2)
    })
  })
  observeEvent(input$IC_cell_clicked, {
    info = input$IC_cell_clicked
    # do nothing if not clicked yet, or the clicked cell is not in the 1st column
    if (is.null(info$value)||info$col != 1) return()
    output$drm <- renderPlotly({
      ds <- paste0(info$value)
      ps <- paste0(input$pat,collapse = "|")
      
      ggfits <- plot.res %>% filter(grepl(ds,Treatment_Drug),grepl(ps,Patient))
      ggdf <- dat.list %>% filter(grepl(ds,Treatment_Drug),grepl(ps,Patient))
      
      g <<- plot.drug(ggdf,ggfits)
    })
    output$drug_title <- renderText({info$value})
    output$raw <- DT::renderDataTable({
      DT::datatable(dat.list[dat.list$Treatment_Drug == info$value|dat.list$Treatment_Drug == "DMSO",
                             c("Patient","Well","Treatment_Drug","Concentration_1","Unit_1","Concentration_2",
                               "Unit_2","Remark","Count","Percentage")] ) %>%
        formatRound(10, 2)
    })
  })
  

  # Create Downloadable data: ECs and ICs
  output$dl1 <- downloadHandler(
      filename = function() {
        paste0(paste('ED',input$ptid, Sys.Date(), sep='_'), '.csv')
      },
      content = function(con) {
        write.csv(ed2, con)
      }
    )
  output$dl2 <- downloadHandler(
    filename = function() {
      paste0(paste('IC',input$ptid, Sys.Date(), sep='_'), '.csv')
    },
    content = function(con) {
      write.csv(ic2, con)
    }
  )
  output$dl4 <- downloadHandler(
    filename = function() {
      paste0(paste('data',input$ptid, Sys.Date(), sep='_'), '.csv')
    },
    content = function(con) {
      write.csv(dat.list, con)
    }
  )
  output$dl3 <- downloadHandler(
    filename = function() {
      if (isTRUE(input$rib)){
        paste0(paste('plot_ci',input$ptid, Sys.Date(), sep='_'), '.zip')
      } else {
        paste0(paste('plot',input$ptid, Sys.Date(), sep='_'), '.zip')
      }
    },
    content = function(plot) {
      if (isTRUE(input$rib)){
        path <- lapply(names(exp.plot),function(x) plotsave(exp.plot.rib[[x]], x))
      } else {
        path <- lapply(names(exp.plot),function(x) plotsave(exp.plot[[x]], x))
      }
      path <- unlist(path)

      zip(zipfile = plot, files = path, extras = '-j')
    }
  )
    
  # Download displayed plot
  plotInput = function() {
    return(g)
  }
  output$shown <- downloadHandler(
    # file name
    filename <- "Export.png",
    # content
    content = function(file){
      ggsave(file, plot = plotInput(), device = "png", scale = 1,
             width = 1800,
             height = 1000,
             units = "px",
             dpi = 300)
    }
  )
}

shinyApp(ui = ui, server = server)
