anova_ui <- function() {
  tagList(
    selectInput("anova_var", "Grouping Variable", choices = NULL),
    radioButtons("anova_var_type", "Variable Type", choices = c("Character", "Factor", "Numeric")),
    selectInput("num_covariates", "Number of Covariates", choices = 0:4),
    uiOutput("covariate_inputs"),
    actionButton("run_anova", "Run ANOVA", class = "btn-primary"),
    DTOutput("anova_output"),
    downloadButton("download_anova", "Download Results", class = "btn-success")
  )
}