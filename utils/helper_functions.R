
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) {
    stop("No input parameters file provided")
  }
  return(args[1])
}

.save_results <- function(results, filename) {
  # Automatically collapse list-columns
  is_list_col <- sapply(results, is.list)
  if (any(is_list_col)) {
    results[is_list_col] <- lapply(results[is_list_col], function(col) {
      sapply(col, function(x) {
        if (is.null(x)) return(NA)
        if (length(x) == 0) return("")
        paste(as.character(x), collapse = ";")
      })
    })
  }
  
  # Define file path
  filepath <- file.path(output_dir, filename)
  
  # Save the file
  if (!file.exists(filepath)) {
    write.csv(results, file = filepath, row.names = FALSE)
  } else {
    write.table(results, file = filepath, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  }
}

quietly <- function(expr) {
  val <- NULL
  capture.output(
    val <- force(expr),
    file = NULL
  )
  val
}