mpl.engine.ppp = function (Q, trend = ~1, interaction = NULL, ..., covariates = NULL, 
                           subsetexpr = NULL, clipwin = NULL, covfunargs = list(), 
                           correction = "border", rbord = 0, use.gam = FALSE, gcontrol = list(), 
                           GLM = NULL, GLMfamily = NULL, GLMcontrol = NULL, famille = NULL, 
                           forcefit = FALSE, nd = NULL, eps = eps, allcovar = FALSE, 
                           callstring = "", precomputed = NULL, savecomputed = FALSE, 
                           preponly = FALSE, rename.intercept = TRUE, justQ = FALSE, 
                           weightfactor = NULL) 
{
  GLMname <- if (!missing(GLM)) 
    short.deparse(substitute(GLM))
  else NULL
  if (!is.null(precomputed$Q)) {
    Q <- precomputed$Q
    X <- precomputed$X
    P <- precomputed$U
  }
  else {
    if (verifyclass(Q, "quad", fatal = FALSE)) {
      validate.quad(Q, fatal = TRUE, repair = FALSE, announce = TRUE)
      X <- Q$data
    }
    else if (verifyclass(Q, "ppp", fatal = FALSE)) {
      X <- Q
      Q <- quadscheme.ppp(X, method = "d")
    }
    else stop("First argument Q should be a point pattern or a quadrature scheme")
    P <- union.quad(Q)
  }
  if (!is.null(clipwin)) {
    if (is.data.frame(covariates)) 
      covariates <- covariates[inside.owin(P, w = clipwin), 
                               , drop = FALSE]
    Q <- Q[clipwin]
    X <- X[clipwin]
    P <- P[clipwin]
  }
  if (justQ) 
    return(Q)
  computed <- if (savecomputed) 
    list(X = X, Q = Q, U = P)
  else NULL
  if (!is.null(trend) && !inherits(trend, "formula")) 
    stop(paste("Argument", sQuote("trend"), "must be a formula"))
  if (!is.null(interaction) && !inherits(interaction, "interact")) 
    stop(paste("Argument", sQuote("interaction"), "has incorrect format"))
  check.1.real(rbord, "In ppm")
  explain.ifnot(rbord >= 0, "In ppm")
  if (correction != "border") 
    rbord <- 0
  covfunargs <- as.list(covfunargs)
  if (is.null(trend)) {
    trend <- ~1
    environment(trend) <- parent.frame()
  }
  want.trend <- !identical.formulae(trend, ~1)
  want.inter <- !is.null(interaction) && !is.null(interaction$family)
  spv <- package_version(versionstring.spatstat())
  the.version <- list(major = spv$major, minor = spv$minor, 
                      release = spv$patchlevel, date = "$Date: 2019/02/20 03:34:50 $")
  if (want.inter) {
    if (outdated.interact(interaction)) 
      interaction <- update(interaction)
  }
  if (!want.trend && !want.inter && !forcefit && !allcovar && 
      is.null(subsetexpr)) {
    npts <- npoints(X)
    W <- as.owin(X)
    if (correction == "border" && rbord > 0) {
      npts <- sum(bdist.points(X) >= rbord)
      areaW <- eroded.areas(W, rbord)
    }
    else {
      npts <- npoints(X)
      areaW <- area(W)
    }
    volume <- areaW * markspace.integral(X)
    lambda <- npts/volume
    co <- log(lambda)
    varcov <- matrix(1/npts, 1, 1)
    fisher <- matrix(npts, 1, 1)
    se <- sqrt(1/npts)
    tag <- if (rename.intercept) 
      "log(lambda)"
    else "(Intercept)"
    names(co) <- tag
    dimnames(varcov) <- dimnames(fisher) <- list(tag, tag)
    maxlogpl <- if (npts == 0) 
      0
    else npts * (log(lambda) - 1)
    rslt <- list(method = "mpl", fitter = "exact", projected = FALSE, 
                 coef = co, trend = trend, interaction = NULL, fitin = fii(), 
                 Q = Q, maxlogpl = maxlogpl, satlogpl = NULL, internal = list(computed = computed, 
                                                                              se = se), covariates = mpl.usable(covariates), 
                 covfunargs = covfunargs, subsetexpr = NULL, correction = correction, 
                 rbord = rbord, terms = terms(trend), fisher = fisher, 
                 varcov = varcov, version = the.version, problems = list())
    class(rslt) <- "ppm"
    return(rslt)
  }
  prep <- mpl.prepare.ppp(Q, X, P, trend, interaction, covariates, 
                      want.trend, want.inter, correction, rbord, "quadrature points", 
                      callstring, subsetexpr = subsetexpr, allcovar = allcovar, 
                      precomputed = precomputed, savecomputed = savecomputed, 
                      covfunargs = covfunargs, weightfactor = weightfactor, 
                      ...)
  if (preponly) {
    prep$info <- list(want.trend = want.trend, want.inter = want.inter, 
                      correction = correction, rbord = rbord, interaction = interaction)
    return(prep)
  }
  fmla <- prep$fmla
  glmdata <- prep$glmdata
  problems <- prep$problems
  likelihood.is.zero <- prep$likelihood.is.zero
  is.identifiable <- prep$is.identifiable
  computed <- resolve.defaults(prep$computed, computed)
  IsOffset <- prep$IsOffset
  if (!is.null(prep$covariates)) 
    covariates <- prep$covariates
  if (!is.identifiable) 
    stop(paste("in", callstring, ":", problems$unidentifiable$print), 
         call. = FALSE)
  .mpl.W <- glmdata$.mpl.W
  .mpl.SUBSET <- glmdata$.mpl.SUBSET
  if (is.null(gcontrol)) 
    gcontrol <- list()
  else stopifnot(is.list(gcontrol))
  gcontrol <- if (!is.null(GLMcontrol)) 
    do.call(GLMcontrol, gcontrol)
  else if (want.trend && use.gam) 
    do.call(mgcv::gam.control, gcontrol)
  else do.call(stats::glm.control, gcontrol)
  if (is.null(GLM) && is.null(famille)) {
    if (want.trend && use.gam) {
      FIT <- gam(fmla, family = quasi(link = "log", variance = "mu"), 
                 weights = .mpl.W, data = glmdata, subset = .mpl.SUBSET, 
                 control = gcontrol)
      fittername <- "gam"
    }
    else {
      FIT <- glm(fmla, family = quasi(link = "log", variance = "mu"), 
                 weights = .mpl.W, data = glmdata, subset = .mpl.SUBSET, 
                 control = gcontrol, model = FALSE)
      fittername <- "glm"
    }
  }
  else if (!is.null(GLM)) {
    fam <- GLMfamily %orifnull% quasi(link = "log", variance = "mu")
    FIT <- GLM(fmla, family = fam, weights = .mpl.W, data = glmdata, 
               subset = .mpl.SUBSET, control = gcontrol)
    fittername <- GLMname
  }
  else {
    if (is.function(famille)) 
      famille <- famille()
    stopifnot(inherits(famille, "family"))
    if (want.trend && use.gam) {
      FIT <- gam(fmla, family = famille, weights = .mpl.W, 
                 data = glmdata, subset = .mpl.SUBSET, control = gcontrol)
      fittername <- "experimental"
    }
    else {
      FIT <- glm(fmla, family = famille, weights = .mpl.W, 
                 data = glmdata, subset = .mpl.SUBSET, control = gcontrol, 
                 model = FALSE)
      fittername <- "experimental"
    }
  }
  environment(FIT$terms) <- sys.frame(sys.nframe())
  co <- FIT$coef
  W <- glmdata$.mpl.W
  SUBSET <- glmdata$.mpl.SUBSET
  Z <- is.data.ppp(Q)
  Vnames <- prep$Vnames
  satlogpl <- -(sum(log(W[Z & SUBSET])) + sum(Z & SUBSET))
  maxlogpl <- if (likelihood.is.zero) 
    -Inf
  else (satlogpl - deviance(FIT)/2)
  fitin <- if (want.inter) 
    fii(interaction, co, Vnames, IsOffset)
  else fii()
  unitname(fitin) <- unitname(X)
  rslt <- list(method = "mpl", fitter = fittername, projected = FALSE, 
               coef = co, trend = trend, interaction = if (want.inter) interaction else NULL, 
               fitin = fitin, Q = Q, maxlogpl = maxlogpl, satlogpl = satlogpl, 
               internal = list(glmfit = FIT, glmdata = glmdata, Vnames = Vnames, 
                               IsOffset = IsOffset, fmla = fmla, computed = computed, 
                               vnamebase = prep$vnamebase, vnameprefix = prep$vnameprefix), 
               covariates = mpl.usable(covariates), covfunargs = covfunargs, 
               subsetexpr = subsetexpr, correction = correction, rbord = rbord, 
               terms = terms(trend), version = the.version, problems = problems)
  class(rslt) <- "ppm"
  return(rslt)
}
