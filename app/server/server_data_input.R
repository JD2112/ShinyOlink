data_input_server <- function(input, output, session, merged_data, var_key_merged) {
  observeEvent(input$merge_data, {
    withProgress(message = 'Merging data...', value = 0, {
      req(input$npx_file, input$var_file)
      
      npx_data <- read_NPX(input$npx_file$datapath)
      var_data <- read_csv(input$var_file$datapath, col_types = cols()) %>%
        distinct(SUBJID, .keep_all = TRUE)
      
      if("SUBJID" %in% colnames(var_data)) {
        var_data$SUBJID <- gsub("\n", "", var_data$SUBJID)
      }
      
      if (!is.null(input$key_file)) {
        key_data <- read_csv(input$key_file$datapath, col_types = cols()) %>%
          distinct(SampleID, .keep_all = TRUE)
        key_data$SampleID <- as.character(key_data$SampleID)
        
        merged <- npx_data %>%
          left_join(key_data, by = "SampleID") %>%
          left_join(var_data, by = "SUBJID")
      } else {
        merged <- npx_data %>%
          mutate(SUBJID = SampleID) %>%
          left_join(var_data, by = "SUBJID")
      }
      
      merged_data(merged)
      
      updateSelectInput(session, "pca_var", choices = colnames(merged))
      updateSelectInput(session, "ttest_var", choices = colnames(merged))
      updateSelectInput(session, "anova_var", choices = colnames(merged))
      updateSelectInput(session, "volcano_var", choices = colnames(merged))
      updateSelectInput(session, "violin_group", choices = colnames(merged))
      updateSelectInput(session, "violin_protein", choices = unique(merged$Assay))
      updateSelectInput(session, "normality_protein", choices = unique(merged$Assay))

      incProgress(1)
    })
  })

  observeEvent(input$merge_var_key, {
    withProgress(message = 'Merging var and key data...', value = 0, {
      req(input$var_file)
      
      var_data <- read_csv(input$var_file$datapath, col_types = cols()) %>%
        distinct(SUBJID, .keep_all = TRUE)
      
      if("SUBJID" %in% colnames(var_data)) {
        var_data$SUBJID <- gsub("\n", "", var_data$SUBJID)
      }
      
      if (!is.null(input$key_file)) {
        key_data <- read_csv(input$key_file$datapath, col_types = cols()) %>%
          distinct(SampleID, .keep_all = TRUE)
        key_data$SampleID <- as.character(key_data$SampleID)
        
        var_key_merged(var_data %>% left_join(key_data, by = "SUBJID"))
      } else {
        showNotification("Key file is not provided. Cannot merge.", type = "warning")
        var_key_merged(var_data)
      }
      
      incProgress(1)
    })
  })
}