# Run CytoSpatio Analysis on Remaining Datasets
# This script runs CytoSpatio on datasets that haven't been fully processed yet

# Set your own working directory here
#setwd("")

# Load CytoSpatio function
source("cytospatio.R")

# Define parameters (same for all runs)
TR <- 500
IR <- 100
HR <- 1

# Define datasets to process (50, 60, 70, 80, 90)
datasets <- list(
  list(file = "example/cell_data_percentile_50.csv", output = "output_percentile_50"),
  list(file = "example/cell_data_percentile_60.csv", output = "output_percentile_60"),
  list(file = "example/cell_data_percentile_70.csv", output = "output_percentile_70"),
  list(file = "example/cell_data_percentile_80.csv", output = "output_percentile_80"),
  list(file = "example/cell_data_percentile_90.csv", output = "output_percentile_90")
)

# Run CytoSpatio on each dataset
for (i in seq_along(datasets)) {
  dataset <- datasets[[i]]

  cat("\n")
  cat("========================================\n")
  cat(sprintf("Processing dataset %d of %d\n", i, length(datasets)))
  cat(sprintf("Input: %s\n", dataset$file))
  cat(sprintf("Output: %s\n", dataset$output))
  cat("========================================\n")

  # Create output directory if it doesn't exist
  if (!dir.exists(dataset$output)) {
    dir.create(dataset$output)
  }

  # Run CytoSpatio
  tryCatch({
    cytospatio(
      input_file = dataset$file,
      output_dir = dataset$output,
      TR = TR,
      IR = IR,
      HR = HR
    )
    cat(sprintf("✓ Successfully completed %s\n", dataset$file))
  }, error = function(e) {
    cat(sprintf("✗ Error processing %s: %s\n", dataset$file, e$message))
  })
}

cat("\n")
cat("========================================\n")
cat("ALL ANALYSES COMPLETE!\n")
cat("========================================\n")
