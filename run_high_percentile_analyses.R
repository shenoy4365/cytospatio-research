# Run CytoSpatio Analysis on High Percentile Datasets
# This script runs CytoSpatio on 70, 80, 90 percentile datasets
# (skipping 50 and 60 which are too large for memory)

# Set your own working directory here
#setwd("")

# Load CytoSpatio function
source("cytospatio.R")

# Define parameters (same for all runs)
TR <- 500
IR <- 100
HR <- 1

# Define datasets to process (70, 80, 90 only - smaller datasets)
datasets <- list(
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
cat("\nNote: Percentiles 50 and 60 were skipped due to memory constraints.\n")
cat("These datasets are too large for CytoSpatio to process.\n")
