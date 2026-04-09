# CytoSpatio Example Run Script
# This script runs the example analysis on the provided cell_data.csv

# Set working directory to the cytospatio folder
#setwd("")

# Load the main CytoSpatio function
source("cytospatio.R")

# Create output directory for results
output_dir <- "example_output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Run CytoSpatio on the example data
# Parameters:
#   TR = 500: Look for interactions up to 500 pixels away
#   IR = 100: Check every 100 pixels (creates 5 ranges: 100, 200, 300, 400, 500)
#   HR = 1: Minimum cell separation of 1 pixel

cat("Starting CytoSpatio analysis...\n")
cat("Input file: example/cell_data.csv\n")
cat("Output directory:", output_dir, "\n")
cat("Parameters: TR=500, IR=100, HR=1\n\n")

cytospatio(
  input_file = "example/cell_data.csv",
  output_dir = output_dir,
  TR = 500,
  IR = 100,
  HR = 1
)

cat("\n===== Analysis Complete! =====\n")
cat("Check the", output_dir, "folder for results:\n")
cat("  - Interaction coefficient CSV files\n")
cat("  - Network visualization PNG images\n")
cat("  - Synthetic tissue images\n")
cat("  - Model and data R objects (.Rda files)\n")
