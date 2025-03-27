volcano_plot_server <- function(input, output, session, merged_data, ttest_results, anova_results) {
  # Reactive value to store the current plot
  current_plot <- reactiveVal(NULL)

  observeEvent(input$run_volcano, {
    withProgress(message = 'Generating Volcano plot...', value = 0, {
      if (input$volcano_plot_type == "ttest") {
        req(ttest_results())
        results <- ttest_results()
      } else {
        req(anova_results())
        results <- anova_results()$results
        results <- results %>% mutate(estimate = statistic)
      }
      
      # Generate and store the plot
      plot <- olink_volcano_plot(results)
      current_plot(plot)
      
      output$volcano_plot <- renderPlot({
        current_plot()
      })
      incProgress(1)
    })
  })
  
  output$download_volcano <- downloadHandler(
    filename = function() { 
      paste("volcano_plot_", input$volcano_plot_type, "_", Sys.Date(), ".png", sep="") 
    },
    content = function(file) {
      req(current_plot())
      ggsave(file, plot = current_plot(), device = "png", width = 10, height = 8)
    }
  )
}