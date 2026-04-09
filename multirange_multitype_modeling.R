multirange_multitype_modeling <- function(input_file, output_dir, tr, ir, hr) {
  script_dir = setwd(getwd())
  source(paste0(script_dir, '/modeling/get_devi.R'))
  source(paste0(script_dir, '/modeling/get_formula.R'))
  source(paste0(script_dir, '/modeling/mppm.fit.ppp.R'))
  
  
  cat("Training multirange multitype point process model...\n")
  
  input_file_name <- basename(tools::file_path_sans_ext(input_file))
  concat_quad_file <- file.path(output_dir, paste0(input_file_name, "_concat_quad_TR_", tr, "_IR_", ir, "_HR_", hr, ".Rda"))
  load(concat_quad_file)
  training_original <- Quad_all_all$moadf
  training <- Quad_all_all$moadf
  
  formula_interaction <- colnames(Quad_all_all$moadf)[7:(ncol(Quad_all_all$moadf)-1)]
  formula_interaction <- formula_interaction[-which(formula_interaction == "pattern_ID")]
  formula <- get_formula(formula_interaction)

  marks_real <- Quad_all_all$moadf$marks[which(Quad_all_all$moadf$.mpl.Y != 0)]
  marks_all <- Quad_all_all$moadf$marks
    
  model_train <- mppm.fit.ppp(Data=Quad_all_all, formula) 
  coef <- model_train$Fit$FIT$coefficients
  confid <- model_train$confidence
  resid <- model_train$residuals

  # commented out to save disk space
  #filename <- file.path(output_dir, paste0(input_file_name, "_wholemodel_TR_", tr, "_IR_", ir, "_HR_", hr, ".Rda"))
  #save(model_train, file = filename)

  filename <- file.path(output_dir, paste0(input_file_name, "_model_TR_", tr, "_IR_", ir, "_HR_", hr, ".Rda"))
  save(coef, family, formula, confid, file = filename)
  
  #confidence returns coeff + SE + Confidence Intervals so save "confid" rather than "coef"
  filename <- file.path(output_dir, paste0(input_file_name, "_coef_TR_", tr, "_IR_", ir, "_HR_", hr, ".csv"))
  write.csv(confid, file = filename)
  filename <- file.path(output_dir, paste0(input_file_name, "_resid_TR_", tr, "_IR_", ir, "_HR_", hr, ".csv"))
  write.csv(resid, file = filename)
  cat("Training completed and model saved!\n")

  # avg_devi_per_cell <- get_devi(training_original, coef, fmla, family)
  # 
  # filename <- file.path(output_dir, paste0(input_file_name, "_avg_devi_per_cell_TR_", tr, "_IR_", ir, "_HR_", hr, ".Rda"))
  # save(avg_devi_per_cell, file = filename)
  
}
