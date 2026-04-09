quadrature_scheme_concatenation <- function(input_file, output_dir, tr, ir, hr) {
  cat("Concatenating quadrature schemes...\n")
  input_file_name <- basename(tools::file_path_sans_ext(input_file))
  for (sr in seq(ir, tr, by = ir)) {
    cat(paste0("Interaction range = ", sr, "\n"))
    quad_file <-file.path(output_dir, paste0(input_file_name, "_quad_SR_", sr, "_HR_", hr, ".Rda"))
    load(quad_file)
    current_quad <- Quad_all$moadf
    rm(Quad_all)
    gc()
    for (c in (7:(dim(current_quad)[2]-2))){
      colnames(current_quad)[c] <- paste(colnames(current_quad)[c], sr, sep = 'x')
    }
    if (sr == ir){
      last_current_quad <- current_quad[,7:(dim(current_quad)[2]-2)]
      concatenated_quad <- current_quad[,1:(dim(current_quad)[2]-2)]
      pattern_ID <- current_quad[,(dim(current_quad)[2]-1)]
      caseweight <- current_quad[,dim(current_quad)[2]]
    } else{
      current_quad <- current_quad[,7:(dim(current_quad)[2]-2)]
      concatenated_quad <- cbind(concatenated_quad, current_quad-last_current_quad)
      last_current_quad <- current_quad
    }
  }
  rm(current_quad)
  gc()
  concatenated_quad <- cbind(concatenated_quad, pattern_ID)
  concatenated_quad <- cbind(concatenated_quad, caseweight)
  load(quad_file)
  Quad_all_all <- Quad_all
  Quad_all_all$moadf <- concatenated_quad
  filename <- file.path(output_dir, paste0(input_file_name, "_concat_quad_TR_", tr, "_IR_", ir, "_HR_", hr, ".Rda"))
  save(Quad_all_all, file = filename)
  cat("Quadrature schemes concatenated!\n")
}