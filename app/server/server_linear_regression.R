linear_regression_server <- function(input, output, session, merged_data) {

  regression_results_rv <- reactiveVal(NULL)

  # Generate covariate inputs dynamically
  output$covariate_inputs <- renderUI({
    req(input$num_covariates)
    num <- input$num_covariates

    lapply(seq_len(num), function(i) {
      tagList(
        selectInput(paste0("covariate", i), paste("Select Covariate", i),
                    choices = colnames(merged_data()), selected = NULL),
        radioButtons(paste0("cov_type", i), paste("Covariate", i, "Type"),
                     choices = c("Character", "Factor", "Numeric"), inline = TRUE)
      )
    })
  })

  # Update dependent variable choices based on merged_data
  output$dependent_var_ui <- renderUI({
    req(merged_data())
    numeric_cols <- names(merged_data())[sapply(merged_data(), is.numeric)]
    selectInput("dependent_var", "Select Dependent Variable", choices = numeric_cols)
  })

  # Update covariate selectInput choices when merged_data or num_covariates changes
  observe({
    req(merged_data())
    cols <- colnames(merged_data())

    num_covs <- as.numeric(input$num_covariates %||% 0)
    for (i in seq_len(num_covs)) {
      updateSelectInput(session, paste0("covariate", i), choices = c("None", cols))
    }
  })

  # Perform linear regression when button clicked
  observeEvent(input$run_regression, {
    withProgress(message = 'Running Linear Regression...', value = 0, {
      tryCatch({
        req(merged_data())
        df <- merged_data()

        incProgress(0.1, detail = "Preparing data")

        # Exclude dependent var and covariates from predictors
        excluded_cols <- c(input$dependent_var, unlist(
          lapply(seq_len(input$num_covariates), function(i) input[[paste0("covariate", i)]])
        ))

        # Identify numeric columns not excluded (potential predictors/biomarkers)
        npx_vars <- colnames(df)[sapply(df, is.numeric) & !(colnames(df) %in% excluded_cols)]

        incProgress(0.2, detail = "Processing NPX variables")

        # Standardize NPX variables to Z-score if selected
        if (input$npx_or_zscore == "Z-score") {
          df <- df %>%
            mutate(across(all_of(npx_vars), ~ scale(.), .names = "{.col}_z"))
          npx_vars <- paste0(npx_vars, "_z")
        }

        # Collect covariates selected by user
        covariates <- unlist(lapply(seq_len(input$num_covariates), function(i) input[[paste0("covariate", i)]]))
        covariates <- covariates[!is.na(covariates) & covariates != "" & covariates != "None"]

        incProgress(0.3, detail = "Preprocessing covariates")

        # Preprocess covariate types (Character/Factor/Numeric)
        for (i in seq_along(covariates)) {
          cov <- covariates[i]
          cov_type <- input[[paste0("cov_type", i)]]
          if (!is.null(cov_type)) {
            df[[cov]] <- switch(cov_type,
                                "Character" = as.character(df[[cov]]),
                                "Factor" = as.factor(df[[cov]]),
                                "Numeric" = as.numeric(df[[cov]]),
                                df[[cov]]) # default to original if no match
          }
        }

        dep_var <- input$dependent_var
        if (is.null(dep_var) || dep_var == "") {
          showNotification("Please select a dependent variable.", type = "error")
          return(NULL)
        }

        incProgress(0.5, detail = "Fitting models")

        # Run linear regression models for each NPX variable
        results <- purrr::map_dfr(npx_vars, function(var) {
          safe_dep_var <- paste0("`", dep_var, "`")
          safe_var <- paste0("`", var, "`")
          safe_covariates <- if(length(covariates) > 0) {
            paste0("`", covariates, "`", collapse = " + ")
          } else {
            NULL
          }

          formula_str <- if (!is.null(safe_covariates)) {
            paste(safe_dep_var, "~", paste(c(safe_var, safe_covariates), collapse = " + "))
          } else {
            paste(safe_dep_var, "~", safe_var)
          }
          form <- as.formula(formula_str)

          model <- tryCatch(lm(form, data = df), error = function(e) NULL)
          if (is.null(model)) return(NULL)

          broom::tidy(model) %>%
            filter(term == var) %>%
            mutate(biomarker = gsub("_z$", "", var))
        })

        incProgress(0.9, detail = "Processing results")

        if (nrow(results) > 0) {
          results_clean <- results %>%
            mutate(adj.p.value = p.adjust(p.value, method = "BH")) %>%
            select(biomarker, estimate, std.error, p.value, adj.p.value) %>%
            arrange(adj.p.value)

          regression_results_rv(results_clean)

          output$regression_results <- DT::renderDataTable({
            DT::datatable(results_clean, options = list(pageLength = 15), rownames = FALSE)
          })
        } else {
          regression_results_rv(NULL)
          output$regression_results <- DT::renderDataTable({
            DT::datatable(data.frame(Message = "No valid models or significant results."), rownames = FALSE)
          })
        }

        incProgress(1, detail = "Done")
      }, error = function(e) {
        print(paste("Error in Linear Regression:", e$message))
        showNotification(paste("Error in Linear Regression:", e$message), type = "error")
      })
    })
  })

  # Download handler for regression results CSV
  output$download_regression <- downloadHandler(
    filename = function() {
      paste0("linear_regression_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      results <- regression_results_rv()
      if (!is.null(results)) {
        readr::write_csv(results, file)
      }
    }
  )

  # Observer to print regression results reactive structure (for debugging)
  observe({
    print("Checking regression_results_rv():")
    print(str(regression_results_rv()))
  })
}
