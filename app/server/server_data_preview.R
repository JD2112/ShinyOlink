data_preview_server <- function(input, output, session, merged_data, var_key_merged) {
  output$data_preview <- renderDT({
    req(merged_data())
    datatable(merged_data(),
      extensions = c('Buttons'),
      options = list(
        paging = TRUE,
        searching = TRUE,
        fixedColumns = TRUE,
        autowidth = TRUE,
        ordering = TRUE,
        dom = 'Bflrtip',
        lengthMenu = list(c(10, 25, 50, 100, 1000, -1), 
                          c('10', '25', '50', '100','1000', 'All')),
        scrollX = TRUE,
        buttons = list(
          list(
            extend = "excel", 
            text = "Download current page", 
            filename = "Merged data (current page)",
            exportOptions = list(
              modifier = list(page = "current")
            )
          )
        )
      )
    )
  })

  output$download_full_data <- downloadHandler(
    filename = function() {
      paste("full_dataset_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(merged_data(), file, row.names = FALSE)
    }
  )

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