mppm.fit.ppp = function (Data, fmla, interaction = Poisson(), ..., iformula = NULL, 
                     random = NULL, weights = NULL, use.gam = FALSE, reltol.pql = 0.001, 
                     gcontrol = list()) 
{

  .mpl.SUBSET <- Data$moadf$.mpl.SUBSET
  .mpl.W <- Data$moadf$.mpl.W
  caseweight <- Data$moadf$caseweight

  fitter <- "glm"
  ctrl <- do.call(glm.control, resolve.defaults(gcontrol, list(maxit = 50)))

  FIT <- glm(formula = fmla, family = quasi(link = "log", variance = "mu"), weights = .mpl.W * caseweight, data = subset(Data$moadf, Data$moadf$.mpl.SUBSET == "TRUE"), control = ctrl)
  vcov_matrix <- vcov(FIT)
  raw_residuals = residuals(FIT, type = 'response')
  FIT <- strip_glm(FIT)
  confidence <- compute_se_and_ci(FIT$coefficients, vcov_matrix)
  
  deviants <- deviance(FIT)
  env <- list2env(Data$moadf, parent = sys.frame(sys.nframe()))
  environment(FIT$terms) <- env
  W <- with(Data$moadf, Data$moadf$.mpl.W * Data$moadf$caseweight)
  SUBSET <- Data$moadf$.mpl.SUBSET
  Z <- (Data$moadf$.mpl.Y != 0)
  maxlogpl <- -(deviants/2 + sum(log(W[Z & SUBSET])) + sum(Z & SUBSET))
  Call = Data$Call
  Info = Data$Info
  Vnamelist = Data$Vnamelist
  Isoffsetlist = Data$Isoffsetlist
  Inter = Data$Inter
  formula = Data$formula
  trend = Data$trend
  iformula = Data$iformula
  random = Data$random
  npat = Data$npat
  result <- list(Call = Call, 
                 Info = Info,
                 Fit = list(fitter = fitter, use.gam = use.gam, fmla = fmla, FIT = FIT,  Vnamelist = Vnamelist, Isoffsetlist = Isoffsetlist), 
                 Inter = Inter,
                 formula = formula, trend = trend, iformula = iformula, 
                 random = random, npat = npat, data = Data$data, Y = Data$Y, maxlogpl = maxlogpl, 
                 datadf = Data$datadf, residuals = raw_residuals, confidence = confidence)
  class(result) <- c("mppm", class(result))
  return(result)
}

checkvars <- function(f, b, extra=NULL, bname=short.deparse(substitute(b))){
  fname <- short.deparse(substitute(f))
  fvars <- variablesinformula(f)
  bvars <- if(is.character(b)) b else names(b)
  bvars <- c(bvars, extra)
  nbg <- !(fvars %in% bvars)
  if(any(nbg)) {
    nn <- sum(nbg)
    stop(paste(ngettext(nn, "Variable", "Variables"),
               commasep(dQuote(fvars[nbg])),
               "in", fname,
               ngettext(nn, "is not one of the", "are not"),
               "names in", bname))
  }
  return(NULL)
}

consistentname <- function(x) {
  xnames <- unlist(lapply(x, getElement, name="name"))
  return(length(unique(xnames)) == 1)
}

firstname <- function(z) { z[[1]]$name }

allpoisson <- function(x) all(sapply(x, is.poisson.interact))

marklevels <- function(x) { levels(marks(x)) }

errorInconsistentRows <- function(what, offending) {
  stop(paste("There are inconsistent",
             what,
             "for the",
             ngettext(length(offending), "variable", "variables"),
             commasep(sQuote(offending)),
             "between different rows of the hyperframe 'data'"),
       call.=FALSE)
}

strip_glm <- function(cm) {
  cm$y = NULL
  cm$model = NULL
  cm$residuals = NULL
  cm$fitted.values = NULL
  cm$effects = NULL
  cm$qr$qr = NULL
  cm$linear.predictors = NULL
  cm$weights = NULL
  cm$prior.weights = NULL
  cm$data = NULL
  
  
  cm$family$variance = NULL
  cm$family$dev.resids = NULL
  cm$family$aic = NULL
  cm$family$validmu = NULL
  cm$family$simulate = NULL
  attr(cm$terms,".Environment") = NULL
  attr(cm$formula,".Environment") = NULL
  
  return(cm)
}

compute_se_and_ci <- function(coef_values, vcov_matrix) {
  se <- sqrt(diag(vcov_matrix))  
  
  ci_lower <- coef_values - 1.96 * se  
  ci_upper <- coef_values + 1.96 * se  
  
  results <- data.frame(
    Coefficient = coef_values,
    SE = se,
    CI95_Lower = ci_lower,
    CI95_Upper = ci_upper
  )
  
  return(results)
}

