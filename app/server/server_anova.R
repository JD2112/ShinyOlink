anova_server <- function(input, output, session, merged_data, anova_results) {
  
  # Reactive value to store covariate data
  covariate_data <- reactiveVal(list())
  
  output$covariate_inputs <- renderUI({
    req(input$num_covariates)
    num_covariates <- as.numeric(input$num_covariates)
    
    lapply(1:num_covariates, function(i) {
      tagList(
        selectInput(paste0("covariate", i), paste("Select Covariate", i), 
                    choices = c("None", colnames(merged_data()))),
        radioButtons(paste0("covariate_type", i), paste("Covariate", i, "Type"),
                     choices = c("Character", "Factor", "Numeric"))
      )
    })
  })
  
  # Observer to update covariate_data when inputs change
  observe({
    req(input$num_covariates)
    num_covariates <- as.numeric(input$num_covariates)
    new_covariate_data <- lapply(1:num_covariates, function(i) {
      list(
        name = input[[paste0("covariate", i)]] %||% "None",
        type = input[[paste0("covariate_type", i)]] %||% "Character"
      )
    })
    covariate_data(new_covariate_data)
  })
  
  observeEvent(input$run_anova, {
    withProgress(message = 'Running ANOVA...', value = 0, {
      req(merged_data(), input$anova_var)
      data_for_test <- merged_data()
      
      # Process main ANOVA variable
      data_for_test[[input$anova_var]] <- switch(input$anova_var_type,
        "Factor" = as.factor(data_for_test[[input$anova_var]]),
        "Numeric" = as.numeric(data_for_test[[input$anova_var]]),
        as.character(data_for_test[[input$anova_var]])
      )
      
      covariates <- c()
      covariate_types <- c()
      for(cov in covariate_data()) {
        if(!is.null(cov$name) && cov$name != "None") {
          covariates <- c(covariates, cov$name)
          covariate_types <- c(covariate_types, cov$type)
          # Process covariate
          data_for_test[[cov$name]] <- switch(cov$type,
            "Factor" = as.factor(data_for_test[[cov$name]]),
            "Numeric" = as.numeric(data_for_test[[cov$name]]),
            as.character(data_for_test[[cov$name]])
          )
        }
      }
      
      if(length(covariates) == 0) {
        results <- olink_anova(data_for_test, variable = input$anova_var)
      } else {
        results <- olink_anova(data_for_test, variable = input$anova_var, covariates = covariates)
      }
      
      # Store results and model details in reactive value
      anova_results(list(
        results = results,
        model_details = list(
          variable = input$anova_var,
          variable_type = input$anova_var_type,
          covariates = covariates,
          covariate_types = covariate_types
        )
      ))
      
      output$anova_output <- renderDT({
        datatable(results, options = list(scrollX = TRUE))
      })
      incProgress(1)
    })
  })
  
  output$download_anova <- downloadHandler(
    filename = function() { paste("anova_results_", Sys.Date(), ".xlsx", sep="") },
    content = function(file) {
      req(anova_results())
      write_xlsx(list(anova_results = anova_results()$results), file)
    }
  )
}