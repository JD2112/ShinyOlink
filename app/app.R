library(shiny)
library(bslib)
library(OlinkAnalyze)
library(tidyverse)
library(msigdbr)
library(plotly)
library(viridis)
library(writexl)
library(lme4)
library(lmerTest)
library(shinyjs)
library(DT)
library(broom)
library(purrr)


# Source utility functions
source("utilities.R")

# Source all your existing modules
source("ui/ui_main.R")
source("server/server_main.R")

# Source the new start page modules
source("ui/ui_startpage.R")
source("server/server_startpage.R")

# Source the data input and preview modules
source("ui/ui_data_input.R")
source("ui/ui_data_preview.R")

# Define the main UI
ui <- fluidPage(
  useShinyjs(),
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  uiOutput("page_content"),
  uiOutput("background_style")
)

# Define the main server
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 1000000 * 1024^2)
  # Reactive values to store the merged data and var-key merged data
  merged_data <- reactiveVal(NULL)
  var_key_merged <- reactiveVal(NULL)
  ttest_results <- reactiveVal(NULL)
  
  # Reactive value to store the current page
  current_page <- reactiveVal("start")
  
  # Render the appropriate page
  output$page_content <- renderUI({
    print(paste("Rendering page:", current_page()))
    if (current_page() == "start") {
      start_page_ui()
    } else {
      print("Rendering single panel UI")
      single_ui()
    }
  })
  
  # Render the background style
  output$background_style <- renderUI({
    if (current_page() == "start") {
      tags$style(HTML("
        body {
          background-image: url('ppin.gif');
          background-size: cover;
          background-repeat: no-repeat;
          background-attachment: fixed;
          height: 100vh;
          display: flex;
          justify-content: center;
          align-items: center;
        }
        .start-page-content {
          background-color: rgba(255, 255, 255, 0.8);
          padding: 30px;
          border-radius: 10px;
          text-align: center;
        }
        .radio-buttons-center {
          display: flex;
          flex-direction: column;
          align-items: center;
        }
        .radio-buttons-center .shiny-options-group {
          width: auto;
          padding-left: 80px;
        }
        .radio-buttons-center .radio {
          text-align: left;
        }
        .radio-buttons-center .radio label {
          font-size: 18px;
          padding: 5px 0;
        }
      "))
    } else {
      tags$style(HTML("
        body {
          background-image: none;
          background-color: #ffffff;
        }
      "))
    }
  })
  
  # Start page server logic
  observeEvent(input$start_analysis, {
    print("Start analysis button clicked")
    if (input$analysis_type == "single") {
      print("Switching to main page")
      current_page("main")
    } else {
      showNotification("Multi-panel analysis is not yet implemented.", type = "warning")
    }
  })
  
  # Main analysis server logic
  observeEvent(current_page(), {
    print(paste("Current page changed to:", current_page()))
    if (current_page() == "main") {
      print("Calling server_main function")
      server_main(input, output, session, merged_data, var_key_merged, ttest_results)
    }
  })
  
  # Add the new merge functionality
  observeEvent(input$merge_var_key, {
  req(input$var_file, input$key_file)
  
  tryCatch({
    var_data <- safe_read_csv(input$var_file$datapath)
    key_data <- safe_read_csv(input$key_file$datapath)
    
    print("Var data dimensions:")
    print(dim(var_data))
    print("Key data dimensions:")
    print(dim(key_data))
    
    var_data <- var_data %>% distinct(SUBJID, .keep_all = TRUE)
    key_data <- key_data %>% distinct(SampleID, .keep_all = TRUE)
    
    if("SUBJID" %in% colnames(var_data)) {
      var_data$SUBJID <- gsub("\n", "", var_data$SUBJID)
    }
    
    print("Columns in var_data:")
    print(colnames(var_data))
    print("Columns in key_data:")
    print(colnames(key_data))
    
    if (!"SUBJID" %in% colnames(var_data) || !"SampleID" %in% colnames(key_data)) {
      stop("Required columns 'SUBJID' or 'SampleID' are missing")
    }
    
    merged <- var_data %>% left_join(key_data, by = c("SUBJID" = "SampleID"))
    var_key_merged(merged)
    
    print("Merged data dimensions:")
    print(dim(merged))
    
    show_success("Var and Key data merged successfully")
  }, error = function(e) {
    handle_error(paste("Error merging data:", e$message))
  })
})
  
  # Add the new download handler
  output$download_var_key_data <- downloadHandler(
    filename = function() {
      paste("var_key_merged_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      req(var_key_merged())
      write.csv(var_key_merged(), file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)