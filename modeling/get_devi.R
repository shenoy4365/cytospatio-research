get_devi <- function(individual_data, coef, fmla, family) {
  
  individual_data <- individual_data[individual_data$.mpl.SUBSET, ]
  .mpl.W <- individual_data$.mpl.W
  
  ctrl <- glm.control(maxit = 50)
  devi_list_all <- glm.prep(fmla, family = quasi(link = "log", variance = "mu"), weights = .mpl.W, data = individual_data, start = coef, control = ctrl)
  
  devi_all_avg <- mean(devi_list_all, na.rm = TRUE)
  devi_real_avg <- mean(devi_list_all[individual_data$.mpl.Y != 0], na.rm = TRUE)
  devi_dummy_avg <- mean(devi_list_all[individual_data$.mpl.Y == 0], na.rm = TRUE)
  
  devi_list_cell_type_all_list <- numeric(5)
  devi_list_cell_type_real_list <- numeric(5)
  devi_list_cell_type_dummy_list <- numeric(5)
  
  for (c in 0:4) {
    subset_all <- devi_list_all[individual_data$marks == c]
    subset_real <- subset_all[individual_data$.mpl.Y != 0]
    subset_dummy <- subset_all[individual_data$.mpl.Y == 0]
    
    devi_list_cell_type_all_list[c + 1] <- mean(subset_all, na.rm = TRUE)
    devi_list_cell_type_real_list[c + 1] <- mean(subset_real, na.rm = TRUE)
    devi_list_cell_type_dummy_list[c + 1] <- mean(subset_dummy, na.rm = TRUE)
  }
  
  return(list(
    'devi_all_avg' = devi_all_avg,
    'devi_real_avg' = devi_real_avg,
    'devi_dummy_avg' = devi_dummy_avg,
    'devi_all_cell_type_avg' = devi_list_cell_type_all_list,
    'devi_real_cell_type_avg' = devi_list_cell_type_real_list,
    'devi_dummy_cell_type_avg' = devi_list_cell_type_dummy_list
  ))
}
