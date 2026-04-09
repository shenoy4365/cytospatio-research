get_formula <- function(interactions) {
  formula_string <- ".mpl.Y ~ marks"
  for (interaction in interactions) {
    formula_string <- paste(formula_string, interaction, sep = " + ")
  }
  return(as.formula(formula_string))
}
