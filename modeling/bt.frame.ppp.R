bt.frame.ppp = function (Q, trend = ~1, interaction = NULL, ..., covariates = NULL, 
          correction = "border", rbord = 0, use.gam = FALSE, allcovar = FALSE) 
{
  prep <- mpl.engine.ppp(Q, trend = trend, interaction = interaction, 
                     ..., covariates = covariates, correction = correction, 
                     rbord = rbord, use.gam = use.gam, allcovar = allcovar, 
                     preponly = TRUE, forcefit = TRUE)
  class(prep) <- c("bt.frame", class(prep))
  return(prep)
}
