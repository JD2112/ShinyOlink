volcano_plot_ui <- function() {
  tagList(
    selectInput("volcano_var", "Grouping Variable", choices = NULL),
    radioButtons("volcano_var_type", "Variable Type", choices = c("Character", "Factor")),
    actionButton("run_volcano", "Generate Volcano Plot", class = "btn-primary"),
    plotOutput("volcano_plot"),
    downloadButton("download_volcano", "Download Plot", class = "btn-success")
  )
}