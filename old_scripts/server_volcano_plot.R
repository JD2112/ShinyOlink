volcano_plot_server <- function(input, output, session, merged_data) {
  observeEvent(input$run_volcano, {
    withProgress(message = 'Generating Volcano plot...', value = 0, {
      req(merged_data(), input$volcano_var)
      data_for_test <- merged_data()
      
      if (input$volcano_var_type == "Factor") {
        data_for_test[[input$volcano_var]] <- as.factor(data_for_test[[input$volcano_var]])
      } else {
        data_for_test[[input$volcano_var]] <- as.character(data_for_test[[input$volcano_var]])
      }
      
      if (length(unique(data_for_test[[input$volcano_var]])) == 2) {
        results <- olink_ttest(data_for_test, variable = input$volcano_var)
      } else {
        results <- olink_anova(data_for_test, variable = input$volcano_var)
        results <- results %>% mutate(estimate = statistic)
      }
      
      output$volcano_plot <- renderPlot({
        olink_volcano_plot(results)
      })
      incProgress(1)
    })
  })
  
  output$download_volcano <- downloadHandler(
    filename = function() { paste("volcano_plot_", Sys.Date(), ".png", sep="") },
    content = function(file) {
      ggsave(file, plot = last_plot(), device = "png", width = 10, height = 8)
    }
  )
}