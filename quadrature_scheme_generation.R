quadrature_scheme_generation <- function(input_file, output_dir, sr, hr) {

  cat(paste0("Generating quadrature scheme of range ", sr, "\n"))
  
  # Setting the working directory
  script_dir = setwd(getwd())
  
  # Sourcing necessary scripts
  source(paste0(script_dir, '/modeling/mppm.ppp.R'))
  source(paste0(script_dir, '/modeling/mppm.quad.ppp.R'))
  source(paste0(script_dir, '/modeling/bt.frame.ppp.R'))
  source(paste0(script_dir, '/modeling/quadscheme.ppp.R'))
  source(paste0(script_dir, '/modeling/default.dummy.ppp.R'))
  source(paste0(script_dir, '/modeling/default.n.tiling.ppp.R'))
  source(paste0(script_dir, '/modeling/mpl.engine.ppp.R'))
  source(paste0(script_dir, '/modeling/mpl.prepare.ppp.R'))
  source(paste0(script_dir, '/modeling/evalInteraction.ppp.R'))
  source(paste0(script_dir, '/modeling/quadBlockSizes.ppp.R'))
  source(paste0(script_dir, '/modeling/evalInterEngine.ppp.R'))
  source(paste0(script_dir, '/modeling/evaluate.ppp.R'))
  source(paste0(script_dir, '/modeling/evalPairPotential.ppp.R'))
  

  
  cell_data <- read.csv(input_file)
  input_dir <- dirname(input_file)
  input_file_name <- basename(tools::file_path_sans_ext(input_file))
  cell_data[[3]] <- factor(cell_data[[3]])
  W <- owin(c(min(cell_data[[1]]), max(cell_data[[1]])), c(min(cell_data[[2]]), max(cell_data[[2]])))
  p <- ppp(cell_data[[1]], cell_data[[2]], marks=cell_data[[3]], window=W)
  P <- list()
  P[[1]] <- p
  
  interaction <- list()
  for (i in 1:length(P)) {
    num <- length(unique(P[[i]]$marks))
    sr_matrix <- matrix(rep(sr, num*num), nrow = num)
    hr_matrix <- matrix(rep(hr, num*num), nrow = num)
    interaction[[i]] <- MultiStraussHard(iradii = sr_matrix, hradii = hr_matrix)
  }
  H <- hyperframe(Y = P)
  I <- hyperframe(Interaction = interaction)
  Quad_all <- mppm.quad.ppp(Y ~ marks, data=H, interaction = I)
  filename <- file.path(output_dir, paste0(input_file_name, "_quad_SR_", sr, "_HR_", hr, '.Rda'))
  save(Quad_all, file=filename)
  cat(paste0("Quadrature scheme of range ", sr, " generated!\n"))
}