data_input_ui <- function() {
  sidebar(
    title = "Data Input",
    fileInput("npx_file", "Upload NPX Data (CSV)", accept = c(".csv")),
    fileInput("key_file", "Upload Key File (CSV, Optional)", accept = c(".csv")),
    fileInput("var_file", "Upload Variables File (CSV)", accept = c(".csv")),
    actionButton("merge_data", "Merge Data", class = "btn-primary"),
    actionButton("merge_var_key", "Merge Var and Key Data", class = "btn-secondary"),
    HTML("<hr> <a href='ShinyOlink.txt' target='_blank'> <i class='fa fa-download'> </i> HOW TO USE</a>"),
    br(),
    br(),
    br(),
    HTML(
      "<small><strong>Disclaimer:</strong> Please remember this Shiny app is based on the 
      <a href='https://github.com/Olink-Proteomics/OlinkRPackage' target='_blank'>OlinkAnalyze package</a> 
      from the Olink Developers. Please cite/use the original package.</small>"
    )
  )
}