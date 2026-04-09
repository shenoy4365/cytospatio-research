# cytospatio_nosimulation.R

# Source the function files
source("quadrature_scheme_generation.R")
source("quadrature_scheme_concatenation.R")
source("multirange_multitype_modeling.R")
source("visualization.R")

packages <- c("spatstat", "spatstat.utils", "spatstat.data", "ggplot2", 
              "dplyr", "permute", "data.table", "igraph", "deldir", "sf")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Main function: cytospatio
cytospatio_nosimulation <- function(input_file, output_dir, TR = 500, IR = 100, HR = 1) {

  print("CytoSpatio v0.2.0")

  if (!file.exists(output_dir)) {
    dir.create(file.path(".",output_dir))
  }
   
  # Generate quadrature schemes for each SR value based on TR and IR
  schemes_list <- list()
  for (SR in seq(IR, TR, by = IR)) {
    quadrature_scheme_generation(input_file, output_dir, SR, HR)
  }

  # Concatenate all quadrature schemes
  input_dir = basename(input_file)
  quadrature_scheme_concatenation(input_dir, output_dir, TR, IR, HR)

  # Train the multirange multitype point process model
  multirange_multitype_modeling(input_file, output_dir, TR, IR, HR)

  # Visualize the spatial relationships between cell types
  visualization(input_file, output_dir, TR, IR, HR)
}

