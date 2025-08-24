library(shiny)
library(bslib)

# Source UI components
source("ui/ui_data_input.R")
source("ui/ui_data_preview.R")
source("ui/ui_descriptive_stats.R")
source("ui/ui_normality_test.R")
source("ui/ui_pca_plot.R")
source("ui/ui_ttest.R")
source("ui/ui_violin_plot.R")
source("ui/ui_volcano_plot.R")
source("ui/ui_anova.R")
source("ui/ui_outlier_detection.R")
source("ui/ui_enhanced_visualization.R")
source("ui/ui_normalization.R")
source("ui/ui_pathway_enrichment.R")
source("ui/ui_lod_integration.R")
source("ui/ui_heatmap.R")
source("ui/ui_umap.R")
source("ui/ui_bridge_sample.R")
source("ui/ui_wilcox.R")
source("ui/ui_anova_posthoc.R")
source("ui/ui_lme.R")
source("ui/ui_lme_posthoc.R")
#source("ui/ui_additional_stats.R")
source("ui/ui_lme_stats.R")
source("ui/ui_boxplot.R")
source("ui/ui_distribution_plot.R")
source("ui/ui_lme_plot.R")
source("ui/ui_pathway_heatmap.R")
source("ui/ui_qc_plot.R")
#source("ui/ui_plate_randomization.R")
source("ui/ui_linear_regression.R")


single_ui <- function() {
  page_sidebar(
    theme = bs_theme(version = 5, bootswatch = "flatly"),
    
    tags$head(
      tags$style(HTML("
        body {
          display: flex;
          flex-direction: column;
          min-height: 100vh;
          background-image: url('https://laboratorytalk.com/files/assets/image/3017/microarray.gif');
        background-size: cover;
        background-repeat: no-repeat;
        background-position: center;
        background-attachment: fixed;
        }
        .footer {
          position: relative;
          bottom: 0;
          left: 0;
          right: 0;
          background-color:rgb(212, 233, 252);
          padding: 10px 0;
        }
        .footer-content {
          display: flex;
          justify-content: space-between;
          align-items: center;
          max-width: 1200px;
          margin: 0 auto;
          padding: 0 15px;
        }
        .footer-section {
          flex: 1;
          text-align: center;
        }
        .footer-left {
          text-align: left;
        }
        .footer-right {
          text-align: right;
        }
        /* Custom styles for main tabs */
        .nav-tabs .nav-link {
          background-color: #f8f9fa;
          color: #495057;
        }
        .nav-tabs .nav-link.active {
          background-color: #28a745;
          color: #ffffff;
        }
        /* Custom styles for sub-tabs */
        .nav-pills .nav-link {
          background-color: #e9ecef;
          color: #495057;
        }
        .nav-pills .nav-link.active {
          background-color: #007bff;
          color: #ffffff;
        }
      "))
    ),
    
    title = "olinkWrappeR: Shiny Interface for Olink Data Analysis",
    
    h1("olinkWrappeR: Shiny Interface for Olink Data Analysis", class = "text-center mb-4"),

    sidebar = data_input_ui(),
    
    navset_tab(
      nav_panel("A. Data Preview", data_preview_ui()),
      nav_panel("B. Preprocessing", 
        navset_pill(
          nav_panel("1. Bridge Selector", bridge_sample_ui()),
          nav_panel("2. Normalization", normalization_ui()),
          nav_panel("3. LOD", lod_integration_ui()),
          nav_panel("4. Outlier Detection", outlier_detection_ui())
        )
      ),
      nav_panel("C. Statistical Analysis",
        navset_pill(
          nav_panel("1. Normality Test", normality_test_ui()),
          nav_panel("2. T-Test", ttest_ui()),
          nav_panel("3. Wilcoxon Test", wilcox_ui()),
          nav_panel("4. ANOVA", anova_ui()),
          nav_panel("5. ANOVA Post-hoc", anova_posthoc_ui()),
          nav_panel("6. Linear Mixed Effects", lme_stats_ui()),
          nav_panel("7. LME Post-hoc", lme_posthoc_ui())
        )
      ),
      nav_panel("D. Exploratory Analysis",
        navset_pill(
          nav_panel("1. PCA Plot", pca_plot_ui()),
          nav_panel("2. UMAP Plot", umap_ui())
        )
      ),
      nav_panel("E. Visualization",
        navset_pill(
          nav_panel("1. Box Plot", boxplot_ui()),
          nav_panel("2. Distribution Plot", distribution_plot_ui()),
          nav_panel("3. LME Plot", lme_plot_ui()),
          nav_panel("4. Pathway Heatmap", pathway_heatmap_ui()),
          nav_panel("5. QC Plot", qc_plot_ui()),
          nav_panel("6. Heatmap Plot", heatmap_ui()),
          nav_panel("7. Volcano Plot", volcano_plot_ui()),
          nav_panel("8. Violin Plot", violin_plot_ui())
        )
      ),
      nav_panel("F. Pathway Enrichment Analysis", pathway_enrichment_ui()
      ),
      nav_panel("G. Linear Regression", linear_regression_ui())
    ),
    
    tags$footer(
      class = "footer mt-auto",
      div(
        class = "container footer-content",
        div(
          class = "footer-section footer-left",
          paste0("©", format(Sys.Date(), "%Y"), ", Developed & maintained by Jyotirmoy Das, Bioinformactics Unit, Core Facility & Clinical Genomics Linköping, Linköping University")
        ),
        div(
          class = "footer-section",
          a("GitHub repo: olinkWrappeR", href = "https://github.com/JD2112/shinyolink", target = "_blank")
        ),
        div(
          class = "footer-section footer-right",
          "Version 1.3"
        )
      )
    )
  )
}

ui <- single_ui()
server <- function(input, output, session) {}

shinyApp(ui, server)