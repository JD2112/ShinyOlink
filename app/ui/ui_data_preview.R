data_preview_ui <- function() {
  tagList(
    downloadButton("download_full_data", "Download Full Dataset"),
    downloadButton("download_var_key_data", "Download Var-Key Merged Data"),
    br(),
    br(),
    DTOutput("data_preview")
  )
}