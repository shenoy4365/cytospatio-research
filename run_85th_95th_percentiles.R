# This script runs CytoSpatio on the newly generated 85th and 95th percentile files

# Set your own working directory here
#setwd("")

# Load CytoSpatio function
source("cytospatio.R")

# Define parameters (same as other analyses)
TR <- 500
IR <- 100
HR <- 1

# Define the new percentile datasets to process
datasets <- list(
  list(file = "example/cell_data_percentile_85.csv", output = "output_percentile_85"),
  list(file = "example/cell_data_percentile_95.csv", output = "output_percentile_95")
)

# Run CytoSpatio on each dataset
for (i in seq_along(datasets)) {
  dataset <- datasets[[i]]

  cat("\n")
  cat("========================================\n")
  cat(sprintf("Processing dataset %d of %d\n", i, length(datasets)))
  cat(sprintf("Input: %s\n", dataset$file))
  cat(sprintf("Output: %s\n", dataset$output))
  cat("========================================\n\n")

  # Run CytoSpatio
  cytospatio(dataset$file, dataset$output, TR, IR, HR)

  cat("\n")
  cat(sprintf("✓ Completed: %s\n", dataset$file))
  cat("========================================\n")
}

cat("\n")
cat("========================================\n")
cat("ALL ANALYSES COMPLETE!\n")
cat("========================================\n")
cat("\nOutput directories created:\n")
cat("  1. output_percentile_85/\n")
cat("  2. output_percentile_95/\n")
cat("\nYou can now run the Python analysis script to compare all models.\n")
cat("========================================\n\n")
