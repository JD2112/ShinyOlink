linear_regression_server <- function(input, output, session, merged_data) {

  regression_results_rv <- reactiveVal(NULL)

  # Dependent variable UI
  output$dependent_var_ui <- renderUI({
    req(merged_data())
    numeric_cols <- names(merged_data())[sapply(merged_data(), is.numeric)]
    selectInput("dependent_var", "Select Dependent Variable", choices = numeric_cols)
  })

  # Covariate UI
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

  # Update covariates when num changes
  observe({
    req(merged_data())
    cols <- colnames(merged_data())
    num_covs <- as.numeric(input$num_covariates %||% 0)
    for (i in seq_len(num_covs)) {
      updateSelectInput(session, paste0("covariate", i), choices = c("None", cols))
    }
  })

  # Run regression
  observeEvent(input$run_regression, {
    withProgress(message = 'Running Linear Regression...', value = 0, {
      tryCatch({
        req(merged_data())
        df <- merged_data()

        incProgress(0.1, detail = "Preparing data")
        dep_var <- input$dependent_var
        if (is.null(dep_var) || dep_var == "") {
          showNotification("Please select a dependent variable.", type = "error")
          return(NULL)
        }

        # Collect covariates
        covariates <- unlist(lapply(seq_len(input$num_covariates), function(i) input[[paste0("covariate", i)]]))
        covariates <- covariates[!is.na(covariates) & covariates != "" & covariates != "None"]

        # Set covariate types
        incProgress(0.2, detail = "Processing covariates")
        for (i in seq_along(covariates)) {
          cov_type <- input[[paste0("cov_type", i)]]
          if (!is.null(cov_type)) {
            df[[covariates[i]]] <- switch(cov_type,
              "Character" = as.character(df[[covariates[i]]]),
              "Factor" = as.factor(df[[covariates[i]]]),
              "Numeric" = as.numeric(df[[covariates[i]]]),
              df[[covariates[i]]]
            )
          }
        }

        incProgress(0.4, detail = "Fitting models per biomarker")
        # Build covariate string
        cov_formula <- if (length(covariates) > 0) paste("+", paste(covariates, collapse = " + ")) else ""

        # Run per-assay regressions
        # results <- df %>%
        #   dplyr::group_by(Assay, OlinkID, UniProt, Panel) %>%
        #   tidyr::nest() %>%
        #   dplyr::mutate(
        #     model = purrr::map(
        #       data,
        #       ~ lm(as.formula(paste("NPX ~", dep_var, cov_formula)), data = .x)
        #     ),
        #     tidy_res = purrr::map(model, ~ broom::tidy(.x, conf.int = TRUE))
        #   ) %>%
        #   tidyr::unnest(tidy_res) %>%
        #   dplyr::filter(term == dep_var) %>%
        #   dplyr::ungroup()
        results <- df %>%
  dplyr::group_by(Assay, OlinkID, UniProt, Panel) %>%
  tidyr::nest() %>%
  dplyr::mutate(
    model = purrr::map(
      data,
      ~ {
        # Keep only complete cases for relevant columns
        vars_needed <- c("NPX", dep_var, covariates)
        dat <- .x %>% dplyr::filter(stats::complete.cases(dplyr::select(., all_of(vars_needed))))
        
        # Skip if too few rows
        if (nrow(dat) < 3) return(NULL)
        
        # Build safe formula with backticks
        formula_str <- paste0("`NPX` ~ `", dep_var, "`",
                              if (length(covariates) > 0) paste0(" + `", covariates, "`", collapse = ""))
        
        # Try fitting the model
        tryCatch(
          lm(as.formula(formula_str), data = dat),
          error = function(e) NULL
        )
      }
    ),
    tidy_res = purrr::map(model, ~ if (!is.null(.x)) broom::tidy(.x, conf.int = TRUE) else NULL)
  ) %>%
  tidyr::unnest(tidy_res) %>%
  dplyr::filter(term == dep_var) %>%
  dplyr::ungroup()


        incProgress(0.8, detail = "Processing results")
        if (nrow(results) > 0) {
          results <- results %>%
            dplyr::mutate(Adjusted_pval = p.adjust(p.value, method = "BH")) %>%
            dplyr::select(Assay, OlinkID, UniProt, Panel,
                          estimate, conf.low, conf.high, statistic, p.value, Adjusted_pval)

          regression_results_rv(results)

          output$regression_results <- DT::renderDataTable({
            DT::datatable(results, options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE)
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

  # Download
  output$download_regression <- downloadHandler(
    filename = function() paste0("linear_regression_results_", Sys.Date(), ".csv"),
    content = function(file) {
      results <- regression_results_rv()
      if (!is.null(results)) readr::write_csv(results, file)
    }
  )

  # Debug observer
  observe({
    print("Checking regression_results_rv():")
    print(str(regression_results_rv()))
  })
}
