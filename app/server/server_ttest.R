ttest_server <- function(input, output, session, merged_data, ttest_results) {
  observeEvent(input$run_ttest, {
    print("T-test button clicked")
    withProgress(message = 'Running T-Test...', value = 0, {
      tryCatch({
        req(merged_data(), input$ttest_var)
        data_for_test <- merged_data()
        
        print("T-test input variable:")
        print(input$ttest_var)
        
        if (input$ttest_var_type == "Factor") {
          data_for_test[[input$ttest_var]] <- as.factor(data_for_test[[input$ttest_var]])
        } else {
          data_for_test[[input$ttest_var]] <- as.character(data_for_test[[input$ttest_var]])
        } 
        
        print("Running olink_ttest function")
        results <- olink_ttest(data_for_test, variable = input$ttest_var)
        
        print("T-test results structure:")
        print(str(results))
        
        ttest_results(results)  # Store results in reactive value
        
        print("ttest_results updated")
        
        output$ttest_output <- renderDT({
          datatable(results, 
            options = list(
                         paging = TRUE,
                         searching = TRUE,
                         fixedColumns = TRUE,
                         autowidth = TRUE,
                         ordering = TRUE,
                         dom = 'Bflrtip',
                         lengthMenu = list(c(10, 25, 50, 100, -1), c('10', '25', '50', '100','All')),
                         scrollX = TRUE,
                         buttons = list(
                           list(extend = "excel", text = "Download current page", 
                                filename = "T-test Analysis",
                                exportOptions = list(
                                  modifier = list(page = "current")
                                ))))
                                )
        })
        
        incProgress(1)
        
        print("T-test completed successfully")
      }, error = function(e) {
        print(paste("Error in T-Test:", e$message))
        showNotification(paste("Error in T-Test:", e$message), type = "error")
      })
    })
  })
  
  output$download_ttest <- downloadHandler(
    filename = function() { paste("ttest_results_", Sys.Date(), ".xlsx", sep="") },
    content = function(file) {
      req(ttest_results())
      write_xlsx(list(ttest_results = ttest_results()), file)
    }
  )
  
  # Add an observer to check ttest_results
  observe({
    print("Checking ttest_results:")
    print(str(ttest_results()))
  })
}