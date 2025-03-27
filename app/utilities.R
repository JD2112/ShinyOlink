# utilities.R

# Function to handle errors
handle_error <- function(error_message) {
  showNotification(
    HTML(paste("Error:", error_message)),
    type = "error",
    duration = NULL
  )
}

# Function to show success messages
show_success <- function(message) {
  showNotification(
    message,
    type = "message",
    duration = 5
  )
}

# Function to validate input
validate_input <- function(input, required_fields) {
  for (field in required_fields) {
    if (is.null(input[[field]]) || input[[field]] == "") {
      stop(paste("Please provide", field))
    }
  }
}

# Function to safely read CSV files
safe_read_csv <- function(file_path) {
  tryCatch({
    data <- read_csv(file_path, col_types = cols(.default = "c"))
    if (nrow(data) == 0) {
      stop("The CSV file is empty")
    }
    return(data)
  }, error = function(e) {
    stop(paste("Error reading CSV file:", e$message))
  })
}

# Function to safely execute server modules
safe_call <- function(module_func, ...) {
  tryCatch({
    module_func(...)
  }, error = function(e) {
    message <- paste("Error in", deparse(substitute(module_func)), ":", e$message)
    showNotification(message, type = "error", duration = NULL)
    # Log the error
    print(message)
  })
} 